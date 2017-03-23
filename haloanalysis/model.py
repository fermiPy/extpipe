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
from fermipy import spectrum 
from haloanalysis.utils import Axis, MapND
from haloanalysis.sed import *

from ebltable.tau_from_model import OptDepth
class LogParabolaExpCutoff(spectrum.SpectralFunction):
    """Class that evaluates a function with the parameterization:

    F(x) = p_0 * (x/x_s)^(p_1 - p_2*log(x/x_s) ) * exp(- x/p_3)

    where x_s is a scale parameter.  The `params` array should be
    defined with:
    * params[0] : Prefactor (p_0)
    * params[1] : Index (p_1)
    * params[2] : Curvature (p_2)
    * params[3] : Exponential cut-off energy (p_3)
    """
    def __init__(self, params=None, scale=1.0, extra_params=None):
	params = (params if params is not None else
	np.array([5e-13, -2.0, 0.0, 1E4]))
	super(LogParabolaExpCutoff, self).__init__(params, scale)
	return

    @staticmethod
    def nparam():
	return 4
									    
    @staticmethod
    def params_to_log(params):
	return [np.log10(params[0]),
		params[1], params[2],
		np.log10(params[3])]

    @staticmethod
    def log_to_params(params):
	return [10**params[0],
	    params[1],params[2],
	    10**params[3]]

    @staticmethod
    def _eval_dnde(x, params, scale=1.0, extra_params=None):
	return (params[0] * (x / scale) **
	    (params[1] - params[2] * np.log(x / scale))) * np.exp(-x / params[3])

    @staticmethod
    def _eval_dfde(x, params, scale=1.0, extra_params=None):
	return (params[0] * (x / scale) **
	    (params[1] - params[2] * np.log(x / scale))) * np.exp(-x / params[3])

def scan_igmf_likelihood(casc_model, rows_sed_tev, rows_sed_gev,
                         rows_casc, tab_pars, tab_ebounds,
                         nstep, casc_scale=1.0, casc_r68_scale=1.0,
			 p0 = [1E-13,-1.5,1E7], scale = 1E3, 
			 fint = PLExpCutoff):
    """This function uses the primary and cascade SEDs to scan the
    likelihood space in the IGMF parameter (B,L_coh).

    Parameters
    ----------
    casc_model : `haloanalysis.model.CascModel`
    		cascade model table

    rows_sed_tev : `astropy.table.Table`
    		row of table with all TeV sources

    rows_sed_gev : `astropy.table.Table`
    		row of table with SEDs from likelihood scan

    rows_casc : `astropy.table.Table`
    		row of table with likelihood scan of extended sources
		(in latest catalog version, this is also included in rows_sed_gev)

    tab_pars : `astropy.table.Table`
    		table with likelihood scan parameters of catalog

    tab_ebounds : `astropy.table.Table`
    		table with energy bounds of likelihood scan
    
    nstep : int
    	Number of tested IGMF strengths and coherence lengths

    {options}

    calc_scale : float (default : 1.)

    calc_r68_scale : float (default : 1.)

    p0 : list (default : [1E-13,-1.5,1E7])
    	initial guess for intrinisc spectrum fit parameters
    
    scale : float
        pivot energy of initial spectrum

    fint : `~fermipy.spectrum.Spectrum` or None
        fermipy spectrum function pointer of initial spectrum. 
	If None, use power law with exponential cut- off.
    """
    
    fn = fint(p0,scale=scale)

    # get the scan width from bin 
    # centers of the right axes:
    if casc_model.axes[-2].name == 'log_lcoh':
	lcoh_scan = np.linspace(casc_model.axes[-2].centers[0],
				casc_model.axes[-2].centers[-1],
				nstep)    
    else:
	lcoh_scan = np.linspace(-4,4,nstep)

    if casc_model.axes[-1].name == 'log_igmf':
	igmf_scan = np.linspace(casc_model.axes[-1].centers[0],
				casc_model.axes[-1].centers[-1],
				nstep)    
    else:
	igmf_scan = np.linspace(-20,-12,nstep)    

    bpars = np.meshgrid(lcoh_scan, igmf_scan, indexing='ij')

    sed_prim0 = SED.create_from_row(rows_sed_tev)
    sed_prim1 = SED.create_from_row2(rows_sed_gev, tab_ebounds)
    sed_casc = HaloSED.create_from_fits(rows_casc[0], tab_ebounds, tab_pars)
    hl = CascLike(casc_model, fn, sed_casc, [sed_prim0, sed_prim1])

    sed_prim0_nbin = 25
    
    halo_flux_ul = np.zeros(bpars[0].shape + tab_ebounds['e_ref'].shape)*np.nan

    try:
	fn.nparam()
    except AttributeError: # assume the PLExpCutoff
	def nparam(): return 3
	fn.nparam = nparam

    model_loglike_comp = np.zeros(bpars[0].shape + (3,))*np.nan
    model_dloglike = np.zeros(bpars[0].shape)*np.nan
    model_fit_pars = np.zeros(bpars[0].shape + (fn.nparam(),))*np.nan
    model_prim0_emin = np.zeros((1,1,sed_prim0_nbin))*np.nan
    model_prim0_ectr = np.zeros((1,1,sed_prim0_nbin))*np.nan
    model_prim0_emax = np.zeros((1,1,sed_prim0_nbin))*np.nan
    model_prim0_flux = np.zeros(bpars[0].shape +
                               (sed_prim0_nbin,))*np.nan
    model_prim1_flux = np.zeros(bpars[0].shape +
                               (sed_prim1.axis.nbin,))*np.nan
    model_casc_flux = np.zeros(bpars[0].shape +
                               (sed_casc.axes[1].nbin,))*np.nan
    model_casc_r68 = np.zeros(bpars[0].shape +
                              (sed_casc.axes[1].nbin,))*np.nan
    
    p1 = fn.params

    redshift = rows_sed_tev[0]['REDSHIFT']
    
    # fit over null hypothesis
    if fn.nparam() == 3:
	bounds = [(-14., -8.0),(-3.0, -1.5),(5., 9.)]

    elif fn.nparam() == 4:
	bounds = [(-14., -8.0),(-3.0, -1.5),(0.0, 5.0),(5., 9.)]

    null_fit = hl.fit([redshift,0.0,-12.0],fn.params,method='SLSQP',
                      casc_scale=1E-10, bounds = bounds)
    null_fit = hl.fit([redshift,0.0,-12.0],null_fit[1],method='SLSQP',
                      casc_scale=1E-10, bounds = bounds)

    null_params = null_fit[1]
    
    flux = np.logspace(-14,-9,200)
    ff, ee = np.meshgrid(flux, 10.**sed_casc.axes[1].centers)

    #full_loge = Axis.create_from_centers('eobs',
		    #np.concatenate((sed_casc.axes[1].centers,
		    #sed_prim0.axis.centers)))

    for idx, x in np.ndenumerate(bpars[0]):

        p0 = [redshift,bpars[0][idx], bpars[1][idx]]

	# fit over all B fields. Use best parameters of null fit 
	# as initial guesses
        lnl, p1, o = hl.fit(p0,null_params,method='SLSQP',
                            casc_scale=casc_scale,
                            casc_r68_scale=casc_r68_scale, bounds = bounds)
        lnl, p1, o = hl.fit(p0,p1,method='SLSQP',
                            casc_scale=casc_scale,
                            casc_r68_scale=casc_r68_scale, bounds = bounds)
        dloglike = -(lnl - null_fit[0])
        
	# get the results and theoretical fluxes
        model_dloglike[idx] = dloglike
	model_loglike_comp[idx] = hl._nll_prim + [hl._lnl_casc]
        model_fit_pars[idx][:fn.nparam()] = p1

        sed_prim0_flux = casc_model.prim_flux(fn,p0,p1,
                                              axis_eobs=sed_prim0.axis)
        
        model_prim0_flux[idx][:len(sed_prim0_flux)] = sed_prim0_flux

        model_prim1_flux[idx] = casc_model.prim_flux(fn,p0,p1,
                                                     axis_eobs=sed_prim1.axis)
        model_casc_flux[idx] = casc_model.casc_flux(fn,p0,p1,
                                                    axis_eobs=sed_casc.axes[1])
                                                    #axis_eobs=full_loge)
        model_casc_r68[idx] = casc_model.casc_r68(fn,p0,p1,
                                                  axis_eobs=sed_casc.axes[1])
                                                  #axis_eobs=full_loge)
	# calculate the flux upper limit of the halo for a 
	# halo with same r68 as best fit cascade
	ff, rr = np.meshgrid(flux,model_casc_r68[idx])
	# get a log l grid over large flux array and 
	#r68 array for each energy bin
	lnl = sed_casc.nll_bin(ff,rr)
	# find the bin closet to 4. which corresponds to 2sigma ul
	# use 3.84 for 95% ul
	ul_id = np.argmin(np.abs(2. * (lnl - np.min(lnl, axis = 0)) - 4.), axis = 1 )
	m = np.min(np.abs(2. * (lnl - np.min(lnl, axis = 0)) - 4.), axis = 1) > 0.5
	# get the flux 
	halo_flux_ul[idx] = np.diag(ff[:,ul_id])
	halo_flux_ul[idx][m] = np.max(ff, axis = 1)[m]
        print(idx, dloglike, p0, p1)

        #print bpars, lnl, p1

    model_prim0_emin[0,0][:sed_prim0.axis.nbin] = 10**sed_prim0.axis.lo
    model_prim0_ectr[0,0][:sed_prim0.axis.nbin] = 10**sed_prim0.axis.centers
    model_prim0_emax[0,0][:sed_prim0.axis.nbin] = 10**sed_prim0.axis.hi
        
    cols_dict = OrderedDict()
    cols_dict['name'] = dict(dtype='S32', format='%s', description='name')
    cols_dict['SOURCE_FULL'] = dict(dtype='S32', format='%s', description='name')
    cols_dict['assoc'] = dict(dtype='S32', format='%s', description='assoc')
    cols_dict['redshift'] = dict(dtype='S32', format='%s', description='redshift')
    cols_dict['dloglike'] = dict(dtype='f8', format='%.3f',
                                       shape=model_dloglike.shape)
    cols_dict['loglike_comp'] = dict(dtype='f8', format='%.3f',
                                       shape=model_loglike_comp.shape)
    cols_dict['prim_tev_emin'] = dict(dtype='f8', format='%.3f',
                                      shape=model_prim0_emin.shape)
    cols_dict['prim_tev_ectr'] = dict(dtype='f8', format='%.3f',
                                      shape=model_prim0_ectr.shape)
    cols_dict['prim_tev_emax'] = dict(dtype='f8', format='%.3f',
                                      shape=model_prim0_emax.shape)
    cols_dict['prim_tev_flux'] = dict(dtype='f8', format='%.3f',
                                      shape=model_prim0_flux.shape)
    cols_dict['prim_flux'] = dict(dtype='f8', format='%.3f',
                                  shape=model_prim1_flux.shape)
    cols_dict['casc_flux'] = dict(dtype='f8', format='%.3f',
                                  shape=model_casc_flux.shape)
    cols_dict['casc_r68'] = dict(dtype='f8', format='%.3f',
                                 shape=model_casc_r68.shape)
    cols_dict['src_fit_pars'] = dict(dtype='f8', 
                                     shape=model_dloglike.shape + (fn.nparam(),))
    cols_dict['lcoh'] = dict(dtype='f8', format='%.3f', shape=(nstep,1))
    cols_dict['igmf'] = dict(dtype='f8', format='%.3f', shape=(1,nstep))
    cols_dict['casc_flux'] = dict(dtype='f8', format='%.3f',
                                  shape=(nstep,nstep) + (sed_casc.axes[1].nbin,))
    cols_dict['halo_flux_ul'] = dict(dtype='f8', format='%.3f', 
				    shape=halo_flux_ul.shape)

    
    row = rows_casc[0]
    
    tab_scan = Table([Column(name=k, **v) for k, v in cols_dict.items()])
    row_dict = {}
    row_dict['name'] = row['name']
    row_dict['SOURCE_FULL'] = rows_sed_tev[0]['SOURCE_FULL']
    row_dict['assoc'] = row['assoc']
    row_dict['redshift'] = rows_sed_tev[0]['REDSHIFT']
    row_dict['dloglike'] = model_dloglike
    row_dict['loglike_comp'] = model_loglike_comp
    row_dict['prim_tev_emin'] = model_prim0_emin
    row_dict['prim_tev_ectr'] = model_prim0_ectr
    row_dict['prim_tev_emax'] = model_prim0_emax
    row_dict['prim_tev_flux'] = model_prim0_flux
    row_dict['prim_flux'] = model_prim1_flux
    row_dict['casc_flux'] = model_casc_flux
    row_dict['casc_r68'] = model_casc_r68
    row_dict['src_fit_pars'] = model_fit_pars
    row_dict['lcoh'] = lcoh_scan[:,np.newaxis]
    row_dict['igmf'] = igmf_scan[np.newaxis,:]
    row_dict['halo_flux_ul'] = halo_flux_ul
    tab_scan.add_row([row_dict[k] for k in cols_dict.keys()])

    return tab_scan
    
    #cols_dict = OrderedDict()
    #cols_dict['lcoh'] = dict(dtype='f8', format='%.3f', shape=bpars[0].shape,
    #                         data=bpars[0])
    #cols_dict['igmf'] = dict(dtype='f8', format='%.3f', shape=bpars[1].shape,
    #                         data=bpars[1])
    #tab_scan_grid = Table([Column(name=k, **v) for k, v in cols_dict.items()])
    #return tab_scan, tab_scan_grid

    
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

def expand_array(v0, v1):

    shape = list(v0.shape)
    for i in range(1,len(shape)):
        shape[i] = 1
    
    return v1.reshape(shape)


def interp_flux(flux, axis0, axis1):

    width0 = expand_array(flux,axis0.width)    
    dfde = flux/width0
    flux_interp = np.zeros((len(axis1.centers),) +
                            flux.shape[1:])
    width1 = expand_array(flux_interp,axis1.width)
    
    for idx, t in np.ndenumerate(flux[0,...]):

        sidx = [slice(None)]
        for i in idx:
            sidx += [np.index_exp[i]]

        v = np.interp(axis1.centers,
                      axis0.centers,
                      np.squeeze(dfde[sidx]))
            
        flux_interp[sidx] = expand_array(flux_interp[sidx],v)
        flux_interp[sidx] *= expand_array(flux_interp[sidx],axis1.width)

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
        self._casc_r68 = None
        
#        self._axis_prim = Axis('eobs',
#                               np.log10(sed_prim.specData.ebins),
#                               np.log10(sed_prim.specData.evals))
        self._axis_casc = sed_casc.axes[1]

    @property
    def model(self):
        return self._model
        
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

        nll = []
        for sed in self._sed_prim:        
            prim_flux = self._model.prim_flux(self._fn, p0, p1,
                                              axis_eobs=sed.axis)
            nll += [sed.nll(prim_flux)]

	# save the likelihood for the individual seds
	self._nll_prim = nll

        return reduce(lambda x,y: np.add(x,y),nll)

    def lnl_casc(self, p0, p1, casc_scale=1.0, casc_r68_scale=1.0, cache=True):
        """Evaluate the likelihood for the cascade component.

        """
        
        casc_flux = self._model.casc_flux(self._fn, p0, p1,
                                          axis_eobs=self._axis_casc)

        casc_flux *= casc_scale

        if cache and self._casc_r68 is not None:
            casc_r68 = self._casc_r68
        else:
            casc_r68 = self._model.casc_r68(self._fn, p0, p1,
                                            axis_eobs=self._axis_casc)
            self._casc_r68 = casc_r68

        casc_r68 *= casc_r68_scale
        return self._sed_casc(casc_flux, casc_r68)
            
    def lnl(self, p0, p1, casc_scale=1.0, casc_r68_scale=1.0, cache=True):
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
        
        lnl_prim = self.lnl_prim(p0, p1)
        lnl_casc = self.lnl_casc(p0, p1, casc_scale, casc_r68_scale, cache)

	# save likelihood of cascade width
	self._lnl_casc = lnl_casc

        return lnl_prim + lnl_casc

    def fit(self, p0, p1=None, 
	    method='SLSQP', 
	    casc_scale=1.0, casc_r68_scale=1.0, 
	    cache=True,
	    bounds = [(-14., -8.0),(-3.0, -1.5),(5., 9.)],
	    tol = 1e-4):
        """Perform a fit of the parameters of the primary spectrum while
        holding the IGMF parameters fixed."""

        self._casc_r68 = None
        self.model.set_pars(p0)
        
        def fitfn(params):
            p = self._fn.log_to_params(params)
            v = self.lnl(p0, p, casc_scale=casc_scale,
                         casc_r68_scale=casc_r68_scale,
                         cache=cache)
            #print params, v
            return v
        
        if p1 is None:
            p1 = self._fn.log_params
        else:
            p1 = self._fn.params_to_log(p1)

        o = scipy.optimize.minimize(fitfn, p1,
                                    bounds=bounds,
                                    method=method, tol=tol)

        p1 = self._fn.log_to_params(o['x'])
        self.model.clear_pars()
        return o['fun'], p1, o

        
class CascModel(object):
    """Object representing a library of cascade simulations."""
    
    def __init__(self, axes, casc_flux, prim_flux):
        """
        Parameters
        ----------
        axes : list

        casc_flux : `~haloanalysis.utils.MapND`

        casc_r68 : `~haloanalysis.utils.MapND`

        prim_flux : `~haloanalysis.utils.MapND`

	{options}

	eblmodel : str or None
		If None, use simulation for prim. flux.
		If str, init EBL model with this name and 
		prim. flux will be analytically calculated through
		inj_flux * exp(-tau)
		where tau is the optical depth given by that EBL model
        """
        self._axes = axes
        self._casc_theta_flux = casc_flux
        self._casc_flux = casc_flux.marginalize([2])
        self._prim_flux = prim_flux

        self._casc_flux_cache = self._casc_flux
        self._casc_theta_flux_cache = self._casc_theta_flux
	self._eblmodel = None

        
    @property
    def axes(self):
        return self._axes

    def set_eblmodel(self, eblmodel = 'dominguez'):
	"""Initialize an EBL model to compute the opacity"""
	self._eblmodel = eblmodel
	if not self._eblmodel == None:
	    self._tau = OptDepth.readmodel(eblmodel)

    def set_pars(self,p0):
        """Set the fixed IGMF and source parameters."""
        self._casc_flux_cache = self._casc_flux.slice(np.arange(2,2+len(p0)),p0) 
        self._casc_theta_flux_cache = self._casc_theta_flux.slice(np.arange(3,3+len(p0)),p0)

    def clear_pars(self):
        self._casc_flux_cache = self._casc_flux
        self._casc_theta_flux_cache = self._casc_theta_flux
        
    def make_fits_template(self, fn, p0, p1):
        """Make a FITS template for the given spectral model and input
        paramters."""
        pass
    
    def casc_flux(self, fn, p0, p1=None, axis_eobs=None, sum_theta=True, squeeze=True):
        """Calculate the cascade flux given an injection spectrum and a
        sequence of observed energy bins.

        Parameters
        ----------
        fn : `~fermipy.spectrum.Spectrum`
           Model for the injection spectrum.

        p0 : list        
           List of IGMF and source parameters.  Parameters can be
           passed as scalars or numpy arrays.  All non-scalar
           parameters must have the same shape.

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

        # Interpolate IGMF Params
        if sum_theta is True:
            casc_flux_map = self._casc_flux_cache
            vals = np.meshgrid(self.axes[0].centers,
                               self.axes[1].centers,
                               indexing='ij', sparse=True)

            for i, v in enumerate(vals):
                vals[i] = vals[i][...,np.newaxis]
            
            for i, ax in enumerate(casc_flux_map.axes[2:]):
                p = np.array(p0[i], ndmin=1)
                vals += [p[np.newaxis, np.newaxis, ...]]
            
            casc_flux = casc_flux_map.interp(tuple(vals))
            
        else:
            casc_flux_map = self._casc_theta_flux_cache
            vals = np.meshgrid(self.axes[0].centers,
                               self.axes[1].centers,
                               self.axes[2].centers,
                               indexing='ij', sparse=True)

            for i, v in enumerate(vals):
                vals[i] = vals[i][...,np.newaxis]
            
            for i, ax in enumerate(casc_flux_map.axes[3:]):
                p = np.array(p0[i], ndmin=1)
                vals += [p[np.newaxis, np.newaxis, ...]]

            casc_flux = casc_flux_map.interp(tuple(vals))

            
        # Integrate over injected energy
        inj_flux = expand_array(casc_flux,inj_flux)
        casc_flux = np.sum(inj_flux*casc_flux, axis=0)

        # Remap to binning defined by axis_eobs
        if axis_eobs is not None:
            casc_flux = interp_flux(casc_flux, self.axes[1],
                                    axis_eobs)

        if squeeze:
            casc_flux = np.squeeze(casc_flux)

        return casc_flux

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
        cf = self.casc_flux(fn,p0,p1,axis_eobs,False,False)

        cf = np.cumsum(cf,axis=1)
        cf /= cf[:,-1,...][:,np.newaxis,...]

        shape = (cf.shape[0],) + (1,) + cf.shape[2:]
        
        cf = np.concatenate((np.zeros(shape),cf),axis=1)        
        r68 = np.zeros(shape)
        
        for i, t in np.ndenumerate(cf[:,0,...]):
            idx = (np.index_exp[i[0]],slice(None),np.index_exp[i[1]])
            y = np.interp(0.68,np.squeeze(cf[idx]),self.axes[2].edges)
            r68[idx] = 10**y
            
        return np.squeeze(r68)
        
    def casc_r68_old(self, fn, p0, p1=None, axis_eobs=None):
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

    def prim_flux(self, fn, p0, p1=None, axis_eobs=None, use_analytic = True):
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

	{options}
	use_analytic : bool
	    if True and an eblmodel is specified, use the analytic function 
	    to calculate the prim flux
        """

        inj_flux = np.squeeze(fn.flux(10**self.axes[0].lo,
                                      10**self.axes[0].hi, p1))

	# centers of log energy bins in self.axes 
        vals = np.meshgrid(self.axes[0].centers,
                           indexing='ij', sparse=True)

	# use analytic formular for primary flux
	if not self._eblmodel == None and use_analytic:
	    if not len(p0) == 3:
		raise ValueError("If you use analytic description for prim flux, redshift required")
	    else:
		if not np.isscalar(p0[0]):
		    p0[0] = np.array(p0[0])
		prim_flux = inj_flux * \
				np.exp(-1. * self._tau.opt_depth(p0[0], 
					    10.**(vals[0] - 6.)))
		prim_flux = prim_flux.T
	# else, get it from the simulation
	else:
	    for i, v in enumerate(vals):
		vals[i] = vals[i][...,np.newaxis]
            
	    for i, ax in enumerate(self._prim_flux.axes[1:]):
		p = np.array(p0[i], ndmin=1)
		vals += [p[np.newaxis, ...]]
		
	    prim_flux = self._prim_flux.interp(tuple(vals))
        
	    # Rescale to injected flux
	    inj_flux = expand_array(prim_flux,inj_flux)
	    prim_flux = inj_flux*prim_flux        
        
        # Remap to binning defined by axis_eobs
        if axis_eobs is not None and self._eblmodel is not None and use_analytic:
            prim_flux = interp_flux(prim_flux, self.axes[1], axis_eobs)
        
        return np.squeeze(prim_flux)
        
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
    def create_from_fits(filename, eblmodel = None):

        tab0 = Table.read(filename)
        tab1 = Table.read(filename, 'ENERGIES')
        tab2 = Table.read(filename, 'THETA')
        
        nebin = len(tab1)
        nthbin = len(tab2)
        if 'z' in tab0.columns:        
            model_shape = (len(tab0), np.unique(tab0['lcoh'][0]).shape[0], np.unique(tab0['igmf'][0]).shape[0])
            axis_offset = 1
        else:
            model_shape = (np.unique(tab0['lcoh'][0]).shape[0], np.unique(tab0['igmf'][0]).shape[0])
            axis_offset = 0
	print model_shape, 'here'
            
        emin = np.array(tab1['E_ledge']/1E6)
        emax = np.array(tab1['E_redge']/1E6)
        ectr = np.array(tab1['E_cen']/1E6)

        thmin = np.array(10**tab2['th_ledge'])
        thmax = np.array(10**tab2['th_redge'])
        thctr = np.array(tab2['th_cen'])
        thedges = np.insert(tab2['th_redge'],0,tab2['th_ledge'][0])
        
        domega = np.pi*(thmax**2-thmin**2)*(np.pi/180.)**2
        
        tab_casc_flux = tab0['casc_flux'].reshape(model_shape + (nebin, nebin, nthbin))
        tab_prim_flux = tab0['prim_flux'].reshape(model_shape + (nebin,))
        tab_inj_flux = tab0['inj_flux'].reshape(model_shape + (nebin,))
        tab_igmf = np.unique(np.array(tab0['igmf'].reshape(model_shape)))
        tab_lcoh = np.unique(np.array(tab0['lcoh'].reshape(model_shape)))
            
        log_igmf = np.log10(tab_igmf)
        log_lcoh = np.log10(tab_lcoh)

        data_casc_flux = np.array(tab_casc_flux/tab_inj_flux[..., np.newaxis, np.newaxis])
        data_prim_flux = np.array(tab_prim_flux/tab_inj_flux)

        data_prim_flux_avg = (np.apply_over_axes(np.sum,data_prim_flux,
                                                 axes=[0+axis_offset,1+axis_offset])/
                              (model_shape[0+axis_offset]*model_shape[1+axis_offset]))
        data_prim_flux_avg = data_prim_flux_avg*np.ones(data_prim_flux.shape)
        data_prim_flux = data_prim_flux_avg
        
        # Swap axes so that E_inj and E_obs are the first two
        # dimensions
        data_casc_flux = np.rollaxis(data_casc_flux, 2+axis_offset, 0)
        data_casc_flux = np.rollaxis(data_casc_flux, 3+axis_offset, 1)
        data_casc_flux = np.rollaxis(data_casc_flux, 4+axis_offset, 2)
        data_prim_flux = np.rollaxis(data_prim_flux, 2+axis_offset, 0)

        axis_einj = Axis.create_from_centers('einj', np.log10(ectr))
        axis_eobs = Axis.create_from_centers('eobs', np.log10(ectr))
        axis_theta = Axis('th', thedges)
        axis_lcoh = Axis.create_from_centers('log_lcoh', log_lcoh)
        axis_igmf = Axis.create_from_centers('log_igmf', log_igmf)

        model_axes = [axis_lcoh, axis_igmf]
        axes = [axis_einj, axis_eobs, axis_theta]

        if 'z' in tab0.columns:
	    # too many z values
	    if not tab0['z'].shape[1] == np.prod(model_shape[1:]):
		tabz = tab0['z'][:,:np.prod(model_shape[1:])]
	    else:
		tabz = tab0['z']

            tab_z = np.unique(np.array(tabz.reshape(model_shape)))
            axis_z = Axis.create_from_centers('z', tab_z)            
            model_axes = [axis_z] + model_axes
        
        mapnd_casc_flux = MapND([axis_einj, axis_eobs, axis_theta] + model_axes,
                                data_casc_flux, True)
        mapnd_prim_flux = MapND([axis_einj] + model_axes,
                                data_prim_flux, True)

        return CascModel(axes + model_axes, mapnd_casc_flux, mapnd_prim_flux)
        

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





