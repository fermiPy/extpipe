import numpy as np

import astropy.units as u

from haloanalysis.utils import Axis, MapND


def create_lnl_from_errors(norm, norm_err):

    stephi = np.linspace(0, 1, 11)
    steplo = -np.linspace(0, 1, 11)[1:][::-1]

    loscale = 3 * norm_err
    hiscale = 3 * norm_err
    loscale[loscale > norm] = norm[loscale > norm]
    
    norm_vals_hi = norm[:, np.newaxis] + stephi[np.newaxis, :] * hiscale[:, np.newaxis]
    norm_vals_lo = norm[:, np.newaxis] + steplo[np.newaxis, :] * loscale[:, np.newaxis]

    norm_vals = np.hstack((norm_vals_lo, norm_vals_hi))
    nll_vals = 0.5 * (norm_vals - norm[:, np.newaxis]) ** 2 / norm_err[:, np.newaxis] ** 2
    #norm_vals *= flux[:, np.newaxis] / norm[:, np.newaxis]

    return norm_vals, nll_vals


class SED(object):

    def __init__(self, lnl, emin, ectr, emax, flux, flux_err):
        
        self._emin = emin
        self._emax = emax
        self._ectr = ectr
        self._ebins = np.insert(emax,0,emin[0])
        self._flux = flux
        self._flux_err = flux_err
        self._lnl = lnl

    @property
    def flux(self):
        return self._flux

    @property
    def flux_err(self):
        return self._flux_err
        
    @property
    def ectr(self):
        return self._ectr

    @property
    def ebins(self):
        return self._ebins
        
    @staticmethod
    def create_from_row(row):

        row['E_REF'].unit = 'TeV'
        row['NORM'].unit = 'ph / (m2 TeV s)'
        row['NORM_ERRP'].unit = 'ph / (m2 TeV s)'
        row['NORM_ERRN'].unit = 'ph / (m2 TeV s)'
        
        dfde = row['NORM']
        dfde_unit = u.ph / (u.MeV * u.cm ** 2 * u.s)
        loge = np.log10(np.array(row['E_REF'].to(u.MeV)[0]))

        norm = np.array(row['NORM'].to(dfde_unit)[0])
        norm_errp = np.array(row['NORM_ERRP'].to(dfde_unit)[0])
        norm_errn = np.array(row['NORM_ERRN'].to(dfde_unit)[0])
        norm_err = 0.5 * (norm_errp + norm_errn)

        m = np.isfinite(loge)
        loge = loge[m]
        norm = norm[m]
        norm_err = norm_err[m]

        dloge = loge[1:] - loge[:-1]
        dloge = np.insert(dloge, 0, dloge[0])
        emin = 10 ** (loge - dloge * 0.5)
        emax = 10 ** (loge + dloge * 0.5)
        ectr = 10 ** loge
        deltae = emax - emin
        flux = norm * deltae
        flux_err = norm_err * deltae
        
        # Build Sequence of 1D likelihoods
        
        norm_vals, nll_vals = create_lnl_from_errors(flux, flux_err)
        lnl_maps = []
        
        for i, x in enumerate(loge):
            #print i, x, norm_vals[i]
            axis = Axis.create_from_centers('flux',norm_vals[i])
            lnl_maps += [MapND([axis], nll_vals[i])]

        return SED(lnl_maps, emin, ectr, emax, flux, flux_err)
            
    def __call__(self, flux):
        return self.nll(flux)

    def nll(self, flux):

        nll_sum = np.zeros(flux.shape)
        for i, t in enumerate(flux):
            args = (np.array(t,ndmin=1),)
            nll_sum += self._lnl[i].interp(args)
            
        return np.sum(nll_sum,axis=0)


class HaloSED(object):
    """Representation of a halo SED containing log-likelihood values as a
    function of energy and angular size.

    """

    def __init__(self, sed):

        self._sed = sed

    @property
    def axes(self):
        return self._sed.axes
        
    @staticmethod
    def create_from_fits(row):

        # These are hard-coded for now until we can figure out how to
        # extract them from the file
        axis0 = Axis.create_from_centers('width',np.linspace(-1.125, 1.125, 13))
        axis1 = Axis('eobs',np.linspace(3,5.5,21))
        axis2 = Axis.create_from_centers('eflux',np.linspace(-9, -5, 41))

        sed = MapND([axis0,axis1,axis2],-np.array(row['fit_halo_sed_scan_dlnl']))

        return HaloSED(sed)
        
    def __call__(self, flux, r68):
        return self.nll(flux, r68)

    def nll(self, flux, width):
        """Evaluate the negative log-likelihood function."""
        
        ectr = self._sed.axes[1].centers

        if flux.ndim > 1:
            ectr = np.expand_dims(ectr,-1)

        log_eflux = np.log10(flux*10**ectr)
        log_width = np.log10(width)

        #log_eflux[log_eflux < self.axes[2].lo[0]]
        #log_width[log_width < -1.0] = -1.0
        
        args = (log_width, ectr, log_eflux)
        return np.sum(self._sed.interp(args),axis=0)
        
