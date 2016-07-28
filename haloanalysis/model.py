import itertools
from collections import OrderedDict

import numpy as np
import scipy
from scipy.interpolate import RegularGridInterpolator
from scipy.interpolate import interp1d
from scipy.optimize import brentq

from astropy.table import Table, Column, join
from astropy.modeling.models import Gaussian1D
from astropy.convolution import convolve

from fermipy.spectrum import *
from haloanalysis.utils import Axis, MapND
from haloanalysis.sed import *

def scan_igmf_likelihood(tab, modelfile, sedtab, outputfile, nstep):
    """This function uses the primary and cascade SEDs to scan the
    likelihood space in the IGMF parameter (B,L_coh).

    Parameters
    ----------
    tab : `astropy.table.Table`

    modelfile : str
       FITS file containing the pre-computed IGMF models.

    sedtab : `astropy.table.Table`
    """
    
    fn = PLExpCutoff([1E-11,-1.5,1E6],scale=1E3)

    lcoh_scan = np.linspace(-4,4,nstep)
    igmf_scan = np.linspace(-20,-12,nstep)    
    bpars = np.meshgrid(lcoh_scan, igmf_scan)

    #sed_prim = CastroData.create_from_flux_points(sedfile)
    sed_prim = SED.create_from_row(sedtab)
    sed_casc = HaloSED.create_from_fits(tab)
    hmm = CascModel.create_from_fits(modelfile)
    hl = CascLike(hmm, fn, sed_casc, sed_prim)

    model_lnl = np.zeros(bpars[0].shape)*np.nan
    p1 = fn.params

    for idx, x in np.ndenumerate(bpars[0]):

        p0 = [bpars[0][idx], bpars[1][idx]]
        lnl, p1 = hl.fit(p0,p1,method='SLSQP')
        model_lnl[idx] = lnl
        print(idx, lnl, p0, p1)
        #print bpars, lnl, p1

    cols_dict = OrderedDict()
    cols_dict['name'] = dict(dtype='S32', format='%s', description='name')
    cols_dict['assoc'] = dict(dtype='S32', format='%s', description='assoc')
    cols_dict['igmf_scan_dlnl'] = dict(dtype='f8', format='%.3f',
                                       shape=model_lnl.shape)
    
    tab_scan = Table([Column(name=k, **v) for k, v in cols_dict.items()])
    row_dict = {}
    row_dict['name'] = tab['name']
    row_dict['assoc'] = tab['assoc']
    row_dict['igmf_scan_dlnl'] = model_lnl
    tab_scan.add_row([row_dict[k] for k in cols_dict.keys()])

    cols_dict = OrderedDict()
    cols_dict['lcoh'] = dict(dtype='f8', format='%.3f', shape=bpars[0].shape,
                             data=bpars[0])
    cols_dict['igmf'] = dict(dtype='f8', format='%.3f', shape=bpars[1].shape,
                             data=bpars[1])
    tab_scan_grid = Table([Column(name=k, **v) for k, v in cols_dict.items()])

    return tab_scan, tab_scan_grid

    
    hdulist = fits.HDUList()
    hdulist.append(fits.table_to_hdu(tab_scan))
    hdulist.append(fits.table_to_hdu(tab_scan_grid))

    hdulist[1].name = 'SCAN_DATA'
    hdulist[2].name = 'SCAN_GRID'
    
    hdulist.writeto(outputfile,clobber=True)


def find_eflux_peak(fn, p1, inj_eflux, axis):

    imax = np.argmax(inj_eflux)
    emin = axis.lo[max(imax-1,0)]
    emax = axis.hi[min(imax+1,axis.nbin-1)]
    e2dfde_deriv_emin = np.squeeze(fn.e2dfde_deriv(10**emin,p1))
    e2dfde_deriv_emax = np.squeeze(fn.e2dfde_deriv(10**emax,p1))

    x = np.linspace(emin,emax,11)
    if (np.sign(e2dfde_deriv_emin) != np.sign(e2dfde_deriv_emax)):
        if e2dfde_deriv_emin > 0:                
            epeak = np.interp(0.,np.squeeze(fn.e2dfde_deriv(10**x,p1))[::-1],x[::-1])
        else:
            epeak = np.interp(0.,np.squeeze(fn.e2dfde_deriv(10**x,p1)),x)
    elif e2dfde_deriv_emax > 0:
        epeak = emax
    elif e2dfde_deriv_emin < 0:
        epeak = emin
    else:
        raise Exception('error finding peak')

    return epeak


def interp_flux(flux, axis0, axis1):

    dfde = flux/axis0.width[:, np.newaxis]
    flux_interp = np.zeros((len(axis1.centers),
                            flux.shape[1]))

    for i in range(flux.shape[1]):
            
        flux_interp[:, i] = np.interp(axis1.centers,
                                      axis0.centers,
                                      dfde[:, i])
        flux_interp[:, i] *= axis1.width

    return flux_interp


def interp_r68(r68, axis0, axis1):

    r68_interp = np.zeros((len(axis1.centers),
                           r68.shape[1]))

    for i in range(r68.shape[1]):

        r68_interp[:, i] = np.interp(axis1.centers,
                                     axis0.centers,
                                     r68[:, i])

    return r68_interp


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


class CascLike(object):

    """Class responsible for evaluating the likelihood function for a
    given model of cascade and primary emission:

    ln L = ln L_casc + ln L_prim

    where L_casc is the likelihood for the cascade component and
    L_prim is the likelihood for the primary component.
    """
    def __init__(self, model, fn, sed_casc, sed_prim):
        """
        Constructor.

        Parameters
        ----------
        model : `~haloanalysis.model.CascModel`
           Cascade model object.  This contains precomputed tables for
           the cascade flux and angular size.

        fn : `~fermipy.spectrum.SpectralFunction`

        sed_casc : `~haloanalysis.sed.HaloSED`
           SED for the cascade component.

        sed_prim : `~haloanalysis.sed.SED`
           SED for the primary component.
        
        """
        
        self._model = model
        self._fn = fn
        self._sed_prim = sed_prim
        self._sed_casc = sed_casc
        self._axis_prim = Axis('eobs',np.log10(sed_prim.ebins),
                               np.log10(sed_prim.ectr))
        
#        self._axis_prim = Axis('eobs',
#                               np.log10(sed_prim.specData.ebins),
#                               np.log10(sed_prim.specData.evals))
        self._axis_casc = sed_casc.axes[1]


    def lnl_prim(self, p0, p1):
        """Evaluate the likelihood for the primary spectrum.

        Parameters
        ----------
        p0 : list
           List of IGMF parameters.  Parameters can be passed as
           either scalars or numpy arrays.  All non-scalar parameters
           must have the same shape.

        p1 : `~numpy.ndarray`
           Array of spectral parameters.

        """
        prim_flux = self._model.prim_flux(self._fn, p0, p1,
                                          axis_eobs=self._axis_prim)        
        return self._sed_prim(prim_flux)

    def lnl(self, p0, p1, casc_scale=1.0):
        """Evaluate the total likelihood.

        Parameters
        ----------
        p0 : list
           List of IGMF parameters.  Parameters can be passed as
           either scalars or numpy arrays.  All non-scalar parameters
           must have the same shape.

        p1 : `~numpy.ndarray`
           Array of spectral parameters.

        casc_scale : float
           Scaling factor to be applied to the cascade flux.

        """
        casc_flux = self._model.casc_flux(self._fn, p0, p1,
                                          axis_eobs=self._axis_casc)

        casc_flux *= casc_scale
        
        casc_r68 = self._model.casc_r68(self._fn, p0, p1,
                                        axis_eobs=self._axis_casc)

        lnl_prim = self.lnl_prim(p0, p1)
        lnl_casc = self._sed_casc(casc_flux, casc_r68)
        return lnl_prim + lnl_casc

    def fit(self, p0, p1=None, method='SLSQP', casc_scale=1.0):
        """Perform a fit of the parameters of the primary spectrum while
        holding the IGMF parameters fixed."""
        
        def fitfn(params):
            p = self._fn.log_to_params(params)
            v = self.lnl(p0,p, casc_scale=casc_scale)
            #print params, v
            return v
        
        if p1 is None:
            p1 = self._fn.log_params
        else:
            p1 = self._fn.params_to_log(p1)

        o = scipy.optimize.minimize(fitfn, p1,
                                    bounds=[(-14., None),
                                            (-3.0, -1.5),
                                            (5., 9.)],
                                    method=method, tol=1e-4)

        p1 = self._fn.log_to_params(o['x'])
        return o['fun'], p1

        
class CascModel(object):
    """Object representing a library of cascade simulations."""
    
    def __init__(self, axes, casc_flux, casc_r68, prim_flux):
        """
        Parameters
        ----------
        axes : list

        casc_flux : `~haloanalysis.utils.MapND`

        casc_r68 : `~haloanalysis.utils.MapND`

        prim_flux : `~haloanalysis.utils.MapND`

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

        p0 : list
           List of IGMF parameters.  Parameters can be passed as
           either scalars or numpy arrays.  All non-scalar parameters
           must have the same shape.

        p1 : `~numpy.ndarray`
           Array of spectral parameters.  If none then the parameters
           of the spectral model will be used.

        axis_eobs : `~haloanalysis.utils.Axis`
           Axis defining the binning in observed energy.

        Returns
        -------
        casc_flux : `~numpy.ndarray`
           Flux of the cascade emission represented as an N x M array
           where N is the number of bins in observed energy and M is
           the number of steps in the IGMF parameters.

        """

        # Get Injected Flux
        inj_flux = np.squeeze(fn.flux(10**self.axes[0].lo,
                                      10**self.axes[0].hi, p1))

        p00 = np.array(p0[0], ndmin=1)
        p01 = np.array(p0[1], ndmin=1)
        vals = np.meshgrid(self.axes[0].centers,
                           self.axes[1].centers,
                           indexing='ij', sparse=True)
        
        # Interpolate IGMF Params
        casc_flux = self._casc_flux.interp((vals[0][..., np.newaxis],
                                            vals[1][..., np.newaxis],
                                            p00[np.newaxis, np.newaxis, ...],
                                            p01[np.newaxis, np.newaxis, ...]))

        # Integrate over injected energy
        inj_flux = inj_flux.reshape(inj_flux.shape + (1, 1))
        casc_flux = np.sum(inj_flux*casc_flux, axis=0)

        # Remap to binning defined by axis_eobs
        if axis_eobs is not None:
            casc_flux = interp_flux(casc_flux, self.axes[1],
                                    axis_eobs)

        return np.squeeze(casc_flux)

    def casc_r68(self, fn, p0, p1=None, axis_eobs=None):
        """Calculate the 68% containment radius of the cascade emission given
        an injection spectrum and a sequence of observed energy bins.  

        Parameters
        ----------
        fn : `~fermipy.spectrum.Spectrum`
           Model for the injection spectrum.

        p0 : list
           List of IGMF parameters.  Parameters can be passed as
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
           as an N x M array where N is the number of bins in
           observed energy and M is the number of steps in the IGMF
           parameters.

        """
        # Get Injected Flux
        inj_eflux = np.squeeze(fn.eflux(10**self.axes[0].lo,
                                        10**self.axes[0].hi, p1))

        inj_e2dfde = np.squeeze(fn.e2dfde(10**self.axes[0].centers,
                                          p1))
        
        imax = np.argmax(inj_eflux)
        emin = self.axes[0].lo[max(imax-1,0)]
        emax = self.axes[0].hi[min(imax+1,self.axes[0].nbin-1)]
        e2dfde_deriv_emin = np.squeeze(fn.e2dfde_deriv(10**emin,p1))
        e2dfde_deriv_emax = np.squeeze(fn.e2dfde_deriv(10**emax,p1))

        if isinstance(fn,PLExpCutoff):
            if p1 is not None:            
                epeak = np.log10(p1[2])
            else:
                epeak = np.log10(fn.params[2])
        else:
            epeak = find_eflux_peak(fn,p1,inj_eflux,self.axes[0])            

        epeak = max(epeak,5.0)
        
#        print e2dfde_deriv_emin, e2dfde_deriv_emax
#        print emin,emax, epeak, epeak2, epeak3
            
 #       import matplotlib.pyplot as plt
 #       plt.figure()

 #       plt.plot([emin,emax],[e2dfde_deriv_emin,e2dfde_deriv_emax])

 #       plt.axvline(epeak,color='k')
 #       plt.axvline(epeak2,color='b')
 #       plt.axvline(epeak3,color='r')
        
#        print epeak
            
        p00 = np.array(p0[0], ndmin=1)
        p01 = np.array(p0[1], ndmin=1)

        vals = np.meshgrid(epeak,
                           self.axes[1].centers,
                           indexing='ij', sparse=True)
        
        # Interpolate IGMF Params
        casc_r68 = self._casc_r68.interp((vals[0][..., np.newaxis],
                                          vals[1][..., np.newaxis],
                                          p00[np.newaxis, np.newaxis, ...],
                                          p01[np.newaxis, np.newaxis, ...]))

        casc_r68 = np.squeeze(casc_r68, axis=0)
        
        # Remap to binning defined by axis_eobs
        if axis_eobs is not None:
            casc_r68 = interp_r68(casc_r68, self.axes[1], axis_eobs)
            #np.interp(axis_eobs.centers,
            #                     self.axes[1].centers,
            #                     casc_r68)

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

        inj_flux = np.squeeze(fn.flux(10**self.axes[0].lo,
                                      10**self.axes[0].hi, p1))
        p00 = np.array(p0[0], ndmin=1)
        p01 = np.array(p0[1], ndmin=1)

        x, = np.meshgrid(self.axes[0].centers,
                         indexing='ij', sparse=True)

        # Interpolate IGMF Params
        prim_flux = self._prim_flux.interp((x[:, np.newaxis],
                                            p00[np.newaxis, :],
                                            p01[np.newaxis, :]))

        # Rescale to injected flux
        prim_flux = inj_flux[:, np.newaxis]*prim_flux

        # Remap to binning defined by axis_eobs
        if axis_eobs is not None:
            prim_flux = interp_flux(prim_flux, self.axes[1], axis_eobs)
        
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
        
        einj_axis = Axis.create_from_centers('einj', np.log10(ectr))
        eobs_axis = Axis.create_from_centers('eobs', np.log10(ectr))
        lcoh_axis = Axis.create_from_centers('log_lcoh', log_lcoh[:, 0])
        igmf_axis = Axis.create_from_centers('log_igmf', log_igmf[0, :])

        axes = [einj_axis, eobs_axis, lcoh_axis, igmf_axis]

        mapnd_casc_flux = MapND([einj_axis, eobs_axis, lcoh_axis, igmf_axis],
                                data_casc_flux, True)
        mapnd_casc_r68 = MapND([einj_axis, eobs_axis, lcoh_axis, igmf_axis],
                               data_casc_r68, True)
        mapnd_prim_flux = MapND([einj_axis, lcoh_axis, igmf_axis],
                                data_prim_flux, True)

        return CascModel(axes, mapnd_casc_flux, mapnd_casc_r68,
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

    # Creating a toy model from CascModel
    x = np.linspace(3.0,7.0,21)
    eobs = np.vstack((x[:-1],x[1:]))

    E=Axis('E',x)
    Ep=Axis('Ep',convolve(x, [0.2, 0.6, 0.3]))
    B=Axis('B',np.linspace(-13,-16,9))
    axes=[E,B]

    test=CascModel(axes)
    test.get_eflux(axes)
    test.get_th68(axes)
    cols=test.write_fits('test.fits',axes)
    

    # Comparing the toy model to the data
    test_source_eflux=HaloModelCalc(cols, test_source['spectrum_type']) 
    #test_source_eflux.eflux(eobs, p0, p1)

    #print test_source_eflux





