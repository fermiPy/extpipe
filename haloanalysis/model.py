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

    """Class responsible for evaluating the total likelihood
    function:

    ln L = ln L_LAT + ln L_TEV

    Things to be specified:

    * Parameterization of injection spectrum
    * Object representing TeV SED
    
    """
    def __init__(self, model, fn, sed_lat, sed_tev):

        self._model = model
        self._fn = fn
        self._sed_tev = sed_tev
        self._axis_eobs = Axis('eobs',sed_tev.specData.ebins,
                               sed_tev.specData.evals)
        
    def lnl(self, p0, p1):
        """Evaluate the likelihood."""

        prim_flux = self._model.prim_flux(self._fn, p0, p1,
                                          axis_eobs=self._axis_eobs)

        
       
        lnl = self._sed_tev(prim_flux)
        return lnl

    def fit(self, p0):

        def fitfn(params):
            p = np.array([10**params[0],params[1],10**params[2]])            
            return self.lnl(p0,p)

        #print(fn(self._fn.params))
        
        #p1 = scipy.optimize.fmin(fn, self._fn.params, disp=False, xtol=1e-3)
        p1 = scipy.optimize.minimize(fitfn, [-9,-1.5,7], method='BFGS', tol=1e-3)

        return p1
        
        #print(p1)
        #print(fn(p1))
        
class HaloModelCalc(object):

    def __init__(self,mmap,smodel):

        self._mmap = mmap
        self._smodel = smodel


    def eflux(self,eobs,p0,p1):
        """Evaluate the model energy flux at observed energy x for parameters
        p.

        Parameters
        ----------

        eobs : array
           Min/max energy of bin edges in 2 x N.

        p0 : array
           Model parameters of primary spectrum.

        p1 : array
           Model parameters of pair cascade.

        """

        x = np.linspace(3.0,7.0,21)
        etrue = np.vstack((x[:-1],x[1:]))
        
        # primary spectrum start at ~100 GeV approximate as powerlaw
        #think about this implementation... WHAT is in mind for smodel     
        w = smodel.eflux(etrue,p0)/1E-6


        # Model Map
        # Axis 0 : True Energy
        # Axis 1 : Observed Energy
        # Axis 2 : B_IGMF
        # Axis N ...

        # Interpolate model flux onto grid of true, observed energy
        #data = np.meshgrid(etrue, w, indexing='ij',sparse=True)
        #define what exactly we want to do here. 
        #I interpret this is we have the etrue and eflux grid points. We want
        #the eflux if given an eobs grid point. Output eflux, input point eobs
        #mmap = interpolate(self._mmap,etrue,w) 
        #mmap=RegularGridInterpolator((etrue, w), data)
        interp=interp1d(x,w)
        mmap=interp(eobs)

        print interp, mmap
        
        #bin width in energy 
        detrue = (x[1]-x[0])
        deobs =  (x[1]-x[0])

        print detrue, deobs
        
        # Summation over true energy
        eflux = np.sum(w*mmap*deobs*detrue,axis=0)

        return eflux

    def th68(self,x,p):
        """Evaluate the 68% containment angle at observed energy x for
        parameters p.

        Parameters
        ----------

        x : array
           Observed energy.

        p : array
           Model parameters.

        """
        #weight theta by energy flux, not primary spectrum as in w, it's energy flux. 
        pass
        
class HaloModelMap(object):
    """Object representing a library of cascade simulations."""
    
    def __init__(self, axes, casc_flux, casc_r68, prim_flux):
        """
        Parameters
        ----------
        axes : list

        casc_flux : `~haloanalysis.model.mapnd`

        prim_flux : `~haloanalysis.model.mapnd`

        casc_r68 : `~haloanalysis.model.mapnd`

        """
        self._axes = axes
        self._casc_flux = casc_flux
        self._casc_r68 = casc_r68
        self._prim_flux = prim_flux

    @property
    def axes(self):
        return self._axes

    def casc_flux(self, fn, p0, p1=None, axis_eobs=None):
        """Calculate the cascade flux given an injection spectrum and a
        sequence of observed energy bins.

        Parameters
        ----------
        fn : `~fermipy.spectrum.Spectrum`
           Model for the injection spectrum.

        p0 : `~numpy.ndarray`
           Array of IGMF parameters.

        p1 : `~numpy.ndarray`
           Array of spectral parameters.  If none then the parameters
           of the spectral model will be used.

        axis_eobs : `~haloanalysis.model.Axis`
           Axis defining the binning in observed energy.

        """

        # Get Injected Flux
        inj_flux = fn.flux(self.axes[0].lo, self.axes[0].hi, p1)

        x, y = np.meshgrid(self.axes[0].centers,
                           self.axes[1].centers,
                           indexing='ij', sparse=True)

        # Interpolate IGMF Params
        casc_flux = self._casc_flux.interp((x, y, p0[0], p0[1]))

        # Integrate over injected energy
        casc_flux = np.sum(inj_flux[:, np.newaxis]*casc_flux, axis=0)

        # Remap to binning defined by axis_eobs
        if axis_eobs is not None:

            casc_dfde = casc_flux/self.axes[1].width
            casc_flux = np.interp(axis_eobs.centers, self.axes[1].centers,
                                  casc_dfde)
            casc_flux *= axis_eobs.width

        return casc_flux

    def prim_flux(self, fn, p0, p1=None, axis_eobs=None):
        """Calculate the primary flux given an injection spectrum and a
        sequence of observed energy bins.

        Parameters
        ----------
        fn : `~fermipy.spectrum.Spectrum`
           Model for the injection spectrum.

        p0 : `~numpy.ndarray`
           Array of IGMF parameters.

        p1 : `~numpy.ndarray`
           Array of spectral parameters.  If none then the parameters
           of the spectral model will be used.

        axis_eobs : tuple
           Tuple with lower and upper bin edges in observed energy.
           If none then the internal energy binning will be used.
        
        """
        
        inj_flux = fn.flux(self.axes[0].lo,self.axes[0].hi,p1)

        x, = np.meshgrid(self.axes[0].centers,
                         indexing='ij',sparse=True)

        # Interpolate IGMF Params
        prim_flux = self._prim_flux.interp((x,p0[0],p0[1]))
        
        # Rescale to injected flux
        prim_flux = inj_flux*prim_flux

        # Remap to binning defined by axis_eobs
        if axis_eobs is not None:

            prim_dfde = prim_flux/self.axes[1].width
            prim_flux = np.interp(axis_eobs.centers, self.axes[1].centers,
                                  prim_dfde)
            prim_flux *= axis_eobs.width
        
        return prim_flux
        

    def get_eflux(self, axes): 
        # Here are some random constants that we can change
        B0=1E-16
        gamma=0.5

        allaxes=[]
        for axis in axes:
            allaxes.append(axis.centers)
        params = np.meshgrid(*allaxes,indexing='ij')

        # Really stupid function which defines the eflux
        ef=PowerLaw.eval_dfde(10**params[0],[1E-10,-2.0],1E3)*(10**params[1]/B0)**gamma
        self._eflux=ef

    def get_th68(self, axes): 
        # Here are some random constants that we can change
        B0=1E-16
        gamma=0.5

        allaxes=[]
        for axis in axes:
            allaxes.append(axis.centers)
        params = np.meshgrid(*allaxes,indexing='ij')

        # Function which defines the containment angle
        th=(10**params[0]/1E3)**(-0.5)*(10**params[1]/B0)**gamma
        self._th68=th
        
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
        ectr = tab1['E_cen']/1E6
        
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
                                data_prim_flux)

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





