import itertools

import numpy as np
import scipy
from scipy.interpolate import RegularGridInterpolator
from scipy.interpolate import interp1d

from astropy.table import Table, Column, join
from astropy.modeling.models import Gaussian1D
from astropy.convolution import convolve

from fermipy.spectrum import *
from haloanalysis.utils import Axis, MapND


def make_prim_model(inj_spectrum, inj_flux, prim_flux, emin, emax):
    """Generate a model for the EBL-absorbed primary spectrum weighted
    according to the injection spectrum given in the first argument.
    Weighting is performed along the injected energy dimension with
    index specified by the eidx parameter."""
    inj_scale = inj_spectrum.flux(emin, emax)/inj_flux
    return prim_flux*inj_scale


def make_casc_model(inj_spectrum, inj_flux, casc_flux, emin, emax, eidx=1):
    """Generate a model for the secondary cascade spectrum weighted
    according to the injection spectrum given in the first argument.
    Weighting is performed along the injected energy dimension with
    index specified by the eidx parameter.

    Parameters
    ----------
    inj_spectrum : `~fermipy.spectrum.Spectrum`
      Model for the injection spectrum.

    inj_flux : `~numpy.ndarray`
       Array of injected flux values in arbitrary units.

    casc_flux : `~numpy.ndarray`
       Array of injected flux values in arbitrary units.

    emin: `~numpy.ndarray`

    emax: `~numpy.ndarray`

    eidx : int
       Index of injected energy axis.

    Returns
    -------
    casc_flux : `~numpy.ndarray`
       Weighted model for the cascade spectrum.
    """

    inj_scale = inj_spectrum.flux(emin, emax)/inj_flux
    return np.sum(casc_flux*inj_scale[..., np.newaxis], axis=eidx)


class HaloLike(object):

    """Class responsible for evaluating the likelihood function for a
    given model of cascade and primary emission:

    ln L = ln L_casc + ln L_prim

    where L_casc is the likelihood for the cascade component and
    L_prim is the likelihood for the primary component.
    """
    def __init__(self, model, fn, sed_casc, sed_prim):
        """
        Parameters
        ----------
        model : `~haloanalysis.model.HaloModelMap`
        """
        
        self._model = model
        self._fn = fn
        self._sed_prim = sed_prim
        self._sed_casc = sed_casc
        self._axis_prim = Axis('eobs',
                               sed_prim.specData.ebins,
                               sed_prim.specData.evals)
        self._axis_casc = sed_casc.axes[1]
        
    def lnl(self, p0, p1):
        """Evaluate the likelihood.

        Parameters
        ----------
        p0 : `~numpy.ndarray`
           Array of IGMF parameters.

        p1 : `~numpy.ndarray`
           Array of spectral parameters.

        """

        prim_flux = self._model.prim_flux(self._fn, p0, p1,
                                          axis_eobs=self._axis_prim)

        casc_flux = self._model.casc_flux(self._fn, p0, p1,
                                          axis_eobs=self._axis_casc)

        casc_r68 = self._model.casc_r68(self._fn, p0, p1,
                                        axis_eobs=self._axis_casc)

        lnl_prim = self._sed_prim(prim_flux)
        lnl_casc = self._sed_casc(casc_flux, casc_r68)
        return lnl_prim + lnl_casc

    def fit(self, p0):

        def fitfn(params):
            p = np.array([10**params[0],params[1],10**params[2]])            
            return self.lnl(p0,p)
        
        #p1 = scipy.optimize.fmin(fn, self._fn.params, disp=False, xtol=1e-3)
        p1 = scipy.optimize.minimize(fitfn, [-9,-1.5,7], method='BFGS', tol=1e-3)

        return p1

        
class HaloModelMap(object):
    """Object representing a library of cascade simulations."""
    
    def __init__(self, axes, casc_flux, casc_r68, prim_flux):
        """
        Parameters
        ----------
        axes : list

        casc_flux : `~haloanalysis.mapnd.MapND`

        casc_r68 : `~haloanalysis.mapnd.MapND`

        prim_flux : `~haloanalysis.mapnd.MapND`

        """
        self._axes = axes
        self._casc_flux = casc_flux
        self._casc_r68 = casc_r68
        self._prim_flux = prim_flux

    @property
    def axes(self):
        return self._axes

    def make_fits_template(self, fn, p0, p1):
        """Make a FITS template for the given spectral model and input
        paramters."""
        pass
    
    def casc_flux(self, fn, p0, p1=None, axis_eobs=None):
        """Calculate the cascade flux given an injection spectrum and a
        sequence of observed energy bins.

        Parameters
        ----------
        fn : `~fermipy.spectrum.Spectrum`
           Model for the injection spectrum.

        p0 : tuple
           Tuple of IGMF parameters.  Parameters can be passed as
           either scalars or numpy arrays.  All non-scalar parameters
           must have the same dimension.

        p1 : `~numpy.ndarray`
           Array of spectral parameters.  If none then the parameters
           of the spectral model will be used.

        axis_eobs : `~haloanalysis.utils.Axis`
           Axis defining the binning in observed energy.

        """

        # Get Injected Flux
        inj_flux = np.squeeze(fn.flux(self.axes[0].lo, self.axes[0].hi, p1))

        p00 = np.array(p0[0],ndmin=1)
        p01 = np.array(p0[1],ndmin=1)
        
        vals = np.meshgrid(self.axes[0].centers,
                           self.axes[1].centers, 
                           indexing='ij', sparse=True)
        
        # Interpolate IGMF Params
        casc_flux = self._casc_flux.interp((vals[0][...,np.newaxis],
                                            vals[1][...,np.newaxis],
                                            p00[np.newaxis,np.newaxis,...],
                                            p01[np.newaxis,np.newaxis,...]))

        # Integrate over injected energy
        inj_flux = inj_flux.reshape(inj_flux.shape + (1,1))
        
        casc_flux = np.sum(inj_flux*casc_flux, axis=0)

        # Remap to binning defined by axis_eobs
        if axis_eobs is not None:

            casc_dfde = casc_flux/self.axes[1].width
            casc_flux = np.interp(axis_eobs.centers, self.axes[1].centers,
                                  casc_dfde)
            casc_flux *= axis_eobs.width

        return np.squeeze(casc_flux)

    def casc_r68(self, fn, p0, p1=None, axis_eobs=None):
        """Calculate the 68% containment radius of the cascade emission given
        an injection spectrum and a sequence of observed energy bins.

        Parameters
        ----------
        fn : `~fermipy.spectrum.Spectrum`
           Model for the injection spectrum.

        p0 : tuple
           Tuple of IGMF parameters.  Parameters can be passed as
           either scalars or numpy arrays.  All non-scalar parameters
           must have the same dimension.

        p1 : `~numpy.ndarray`
           Array of spectral parameters.  If none then the parameters
           of the spectral model will be used.

        axis_eobs : `~haloanalysis.utils.Axis`
           Axis defining the binning in observed energy.

        Returns
        -------
        casc_r68 : `~numpy.ndarray`
           68% containment radius of the cascade emission represented
           as an N x M x K array where N is the number of bins in
           observed energy, M is the number of steps in the IGMF
           parameters, and K is the number of steps in the spectral
           parameters.

        """
        # Get Injected Flux
        inj_eflux = np.squeeze(fn.eflux(self.axes[0].lo, self.axes[0].hi, p1))
        imax = np.argmax(inj_eflux)

        p00 = np.array(p0[0],ndmin=1)
        p01 = np.array(p0[1],ndmin=1)

        vals = np.meshgrid(self.axes[0].centers[imax],
                           self.axes[1].centers, 
                           indexing='ij', sparse=True)
        
        # Interpolate IGMF Params
        casc_r68 = self._casc_r68.interp((vals[0][...,np.newaxis],
                                          vals[1][...,np.newaxis],
                                          p00[np.newaxis,np.newaxis,...],
                                          p01[np.newaxis,np.newaxis,...]))

        # Remap to binning defined by axis_eobs
        if axis_eobs is not None:
            casc_r68 = np.interp(axis_eobs.centers,
                                 self.axes[1].centers,
                                 casc_r68)

        return np.squeeze(casc_r68)

    def prim_flux(self, fn, p0, p1=None, axis_eobs=None):
        """Calculate the primary flux given an injection spectrum and a
        sequence of observed energy bins.

        Parameters
        ----------
        fn : `~fermipy.spectrum.Spectrum`
           Model for the injection spectrum.

        p0 : tuple
           Tuple of IGMF parameters.  Parameters can be passed as
           either scalars or numpy arrays.  All non-scalar parameters
           must have the same dimension.

        p1 : `~numpy.ndarray`
           Array of spectral parameters.  If none then the parameters
           of the spectral model will be used.

        axis_eobs : `~haloanalysis.utils.Axis`
           Tuple with lower and upper bin edges in observed energy.
           If none then the internal energy binning will be used.
        
        """

        inj_flux = np.squeeze(fn.flux(self.axes[0].lo, self.axes[0].hi, p1))
        p00 = np.array(p0[0], ndmin=1)
        p01 = np.array(p0[1], ndmin=1)
        
        x, = np.meshgrid(self.axes[0].centers,
                         indexing='ij', sparse=True)

        # Interpolate IGMF Params
        prim_flux = self._prim_flux.interp((x[:,np.newaxis],
                                            p00[np.newaxis,:],
                                            p01[np.newaxis,:]))
        
        # Rescale to injected flux
        prim_flux = inj_flux[:,np.newaxis]*prim_flux

        # Remap to binning defined by axis_eobs
        if axis_eobs is not None:

            prim_dfde = prim_flux/self.axes[1].width
            prim_flux = np.interp(axis_eobs.centers, self.axes[1].centers,
                                  prim_dfde)
            prim_flux *= axis_eobs.width
        
        return np.squeeze(prim_flux)

#    def get_eflux(self, axes): 
#        # Here are some random constants that we can change
#        B0=1E-16
#        gamma=0.5
#
#        allaxes=[]
#        for axis in axes:
#            allaxes.append(axis.centers)
#        params = np.meshgrid(*allaxes,indexing='ij')
#
#        # Really stupid function which defines the eflux
#        ef=PowerLaw.eval_dfde(10**params[0],[1E-10,-2.0],1E3)*(10**params[1]/B0)**gamma
#        self._eflux=ef

#    def get_th68(self, axes): 
#        # Here are some random constants that we can change
#        B0=1E-16
#        gamma=0.5
#
#        allaxes=[]
#        for axis in axes:
#            allaxes.append(axis.centers)
#        params = np.meshgrid(*allaxes,indexing='ij')
#
#        # Function which defines the containment angle
#        th=(10**params[0]/1E3)**(-0.5)*(10**params[1]/B0)**gamma
#        self._th68=th
        
    def write_fits(self,filename,axes):

        cols  = [Column(name='eflux',dtype='f8',shape=self.eflux.shape,data=self.eflux)]
        cols += [Column(name='th68',dtype='f8',shape=self.th68.shape,data=self.th68)]
        allaxes=[]
        for axis in axes:
            allaxes.append(axis.centers)

        params = np.meshgrid(*allaxes,indexing='ij')
        for i in range(len(axes)):
            cols += [Column(name=axes[i].name,dtype='f8',data=params[i])]
        tab = Table(cols)        
        tab.write(filename,format='fits',overwrite=True)

        return tab

    @staticmethod
    def create_from_fits(filename):

        tab0 = Table.read(filename)
        tab1 = Table.read(filename, 'ENERGIES')
        
        nebin = 44
        model_shape = (9, 9)

        emin = tab1['E_ledge']/1E6
        emax = tab1['E_redge']/1E6
        ectr = np.array(tab1['E_cen']/1E6)
        
        tab_r_68 = tab0['r_68'].reshape(model_shape + (nebin, nebin))
        tab_casc_flux = tab0['casc_flux'].reshape(model_shape + (nebin, nebin))
        tab_prim_flux = tab0['prim_flux'].reshape(model_shape + (nebin,))
        tab_inj_flux = tab0['inj_flux'].reshape(model_shape + (nebin,))
        tab_igmf = np.array(tab0['igmf'].reshape(model_shape))
        tab_lcoh = np.array(tab0['lcoh'].reshape(model_shape))

        log_igmf = np.log10(tab_igmf)
        log_lcoh = np.log10(tab_lcoh)

        data_casc_flux = np.array(tab_casc_flux/tab_inj_flux[..., np.newaxis])
        data_prim_flux = np.array(tab_prim_flux/tab_inj_flux)
        data_casc_r68 = np.array(tab_r_68)

        # Swap axes so that E_inj and E_obs are the first two
        # dimensions
        data_casc_flux = np.rollaxis(data_casc_flux, 2, 0)
        data_casc_flux = np.rollaxis(data_casc_flux, 3, 1)
        data_casc_r68 = np.rollaxis(data_casc_r68, 2, 0)
        data_casc_r68 = np.rollaxis(data_casc_r68, 3, 1)
        data_prim_flux = np.rollaxis(data_prim_flux, 2, 0)

        data_casc_r68[data_casc_r68 < 1E-6] = 1E-6
        
        einj_axis = Axis.create_from_centers('einj', ectr, True)
        eobs_axis = Axis.create_from_centers('eobs', ectr, True)
        lcoh_axis = Axis.create_from_centers('log_lcoh', log_lcoh[:, 0])
        igmf_axis = Axis.create_from_centers('log_igmf', log_igmf[0, :])

        axes = [einj_axis, eobs_axis, lcoh_axis, igmf_axis]

        mapnd_casc_flux = MapND([einj_axis, eobs_axis, lcoh_axis, igmf_axis],
                                data_casc_flux, True)
        mapnd_casc_r68 = MapND([einj_axis, eobs_axis, lcoh_axis, igmf_axis],
                               data_casc_r68, True)
        mapnd_prim_flux = MapND([einj_axis, lcoh_axis, igmf_axis],
                                data_prim_flux, True)

        return HaloModelMap(axes, mapnd_casc_flux, mapnd_casc_r68,
                            mapnd_prim_flux)
        

if __name__ == '__main__':
    
    from astropy.io import fits 
    import stack
    
    def column(matrix, i):
        return[row[i] for row in matrix]

    # Reading in a fits file using the Table class
    hdulist=Table.read("/u/gl/mdwood/ki20/mdwood/fermi/ext_analysis/v9/table_std_all_3fgl_glat050.fits")
    hdu_lnl=Table.read("/u/gl/mdwood/ki20/mdwood/fermi/ext_analysis/v9/table_std_all_3fgl_glat050_lnl.fits")

    joint=join(hdulist,hdu_lnl)
    test_mask=stack.create_mask(joint, {'assoc':['1ES 0229+200']})

    #print joint['name','assoc'][test_mask]

    test_source=joint[test_mask]

    # Creating a toy model from HaloModelMap
    x = np.linspace(3.0,7.0,21)
    eobs = np.vstack((x[:-1],x[1:]))

    E=Axis('E',x)
    Ep=Axis('Ep',convolve(x, [0.2, 0.6, 0.3]))
    B=Axis('B',np.linspace(-13,-16,9))
    axes=[E,B]

    test=HaloModelMap(axes)
    test.get_eflux(axes)
    test.get_th68(axes)
    cols=test.write_fits('test.fits',axes)
    

    # Comparing the toy model to the data
    test_source_eflux=HaloModelCalc(cols, test_source['spectrum_type']) 
    #test_source_eflux.eflux(eobs, p0, p1)

    #print test_source_eflux





