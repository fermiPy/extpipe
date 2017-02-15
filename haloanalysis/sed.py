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

    def __init__(self, lnl, emin, ectr, emax, flux, flux_err, flux_ul95):
        
        self._emin = emin
        self._emax = emax
        self._ewidth = emax-emin
        self._ectr = ectr
        self._ebins = np.insert(emax,0,emin[0])
        self._flux = flux
        self._flux_err = flux_err
        self._flux_ul95 = flux_ul95
        self._lnl = lnl
        self._axis = Axis('eobs',np.log10(self._ebins),
                          np.log10(self._ectr))

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
    def ewidth(self):
        return self._ewidth
    
    @property
    def ebins(self):
        return self._ebins

    @property
    def axis(self):
        return self._axis
    
    @staticmethod
    def create_from_row2(row, tab_ebounds):

        ref_flux = np.array(row['ref_flux'][0])
        norm = np.array(row['norm'][0])
        norm_err = np.array(row['norm_err'][0])
        norm_ul = np.array(row['norm_ul'][0])
                
        flux = norm*ref_flux
        flux_err = norm_err*ref_flux
        flux_ul = norm_ul*ref_flux

        ectr = np.array(tab_ebounds['e_ref'])
        emin = np.array(tab_ebounds['e_min'])
        emax = np.array(tab_ebounds['e_max'])

        norm_vals = np.array(row['norm_scan'][0])*ref_flux[:,np.newaxis]
        nll_vals = -np.array(row['dloglike_scan'][0])
        lnl_maps = []
        
        for i, x in enumerate(ectr):
            axis = Axis.create_from_centers('flux',norm_vals[i])
            lnl_maps += [MapND([axis], nll_vals[i])]

        return SED(lnl_maps, emin, ectr, emax, flux, flux_err, flux_ul)
        
        
    @staticmethod
    def create_from_row(row):

        row['e_ref'].unit = 'TeV'
        row['norm'].unit = 'ph / (m2 TeV s)'
        row['norm_errp'].unit = 'ph / (m2 TeV s)'
        row['norm_errn'].unit = 'ph / (m2 TeV s)'
        
        dfde = row['norm']
        dfde_unit = u.ph / (u.MeV * u.cm ** 2 * u.s)
        loge = np.log10(np.array(row['e_ref'].to(u.MeV)[0]))

        norm = np.array(row['norm'].to(dfde_unit)[0])
        norm_errp = np.array(row['norm_errp'].to(dfde_unit)[0])
        norm_errn = np.array(row['norm_errn'].to(dfde_unit)[0])
        norm_err = 0.5 * (norm_errp + norm_errn)
        norm_ul = 2.0*norm_errp

        m = np.isfinite(loge)
        loge = loge[m]
        norm = norm[m]
        norm_err = norm_err[m]
        norm_ul = norm_ul[m]

        dloge = loge[1:] - loge[:-1]
        dloge = np.insert(dloge, 0, dloge[0])
        emin = 10 ** (loge - dloge * 0.5)
        emax = 10 ** (loge + dloge * 0.5)
        ectr = 10 ** loge
        deltae = emax - emin
        flux = norm * deltae
        flux_err = norm_err * deltae
        flux_ul = norm_ul * deltae
        
        # Build Sequence of 1D likelihoods
        
        norm_vals, nll_vals = create_lnl_from_errors(flux, flux_err)
        lnl_maps = []
        
        for i, x in enumerate(loge):
            #print i, x, norm_vals[i]
            axis = Axis.create_from_centers('flux',norm_vals[i])
            lnl_maps += [MapND([axis], nll_vals[i])]

        return SED(lnl_maps, emin, ectr, emax, flux, flux_err, flux_ul)
            
    def __call__(self, flux):
        return self.nll(flux)

    def nll(self, flux):

        nll_sum = np.zeros(flux.shape)
        for i, t in enumerate(flux):
            args = (np.array(t,ndmin=1),)
            nll_sum += self._lnl[i].interp(args)
            
        return np.sum(nll_sum,axis=0)

    def plot(self, ax=None, **kwargs):

        import matplotlib.pyplot as plt
        
        ax = ax if ax is not None else plt.gca()

        kwargs.setdefault('linestyle','None')
        kwargs.setdefault('marker','o')

        efct = self._ectr**2/self._ewidth*1.602177e-06

        m = self._flux_err > self._flux
        
        ax.errorbar(self._ectr[~m]/1E6,
                    (efct*self._flux)[~m],
                    (efct*self._flux_err)[~m],
                    **kwargs)
        if np.any(m):

            import copy
            kw = copy.deepcopy(kwargs)
            kw.pop('label')
            
            ax.errorbar(self._ectr[m]/1E6,
                        (efct*self._flux_ul95)[m],
                        0.5*(efct*self._flux_ul95)[m],
                        uplims=True,
                        **kw)
        


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
    def create_from_fits(row, tab_ebounds, tab_pars):

        width = np.log10(np.array(tab_pars['halo_scan_width'][0]))
        index = np.log10(tab_pars['halo_scan_index'][0])
        eobs = np.array(tab_ebounds['E_MAX'])
        eobs = np.insert(eobs,0,tab_ebounds['E_MIN'][0])
        eobs = np.log10(eobs)
        eflux = np.log10(np.array(tab_pars['halo_scan_eflux'][0]))

        axis0 = Axis.create_from_centers('width',width)
        axis1 = Axis('eobs',eobs)
        axis2 = Axis.create_from_centers('eflux',eflux)

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
        
        log_eflux[log_eflux < self.axes[2].edges[0] ] = self.axes[2].edges[0]
        log_width[log_width < self.axes[0].edges[0] ] = self.axes[0].edges[0]
        
        args = (log_width, ectr, log_eflux)

#        print 'log_eflux: ', log_eflux
#        print 'log_width: ', log_width
#        print self._sed.interp(args)
        
        return np.sum(self._sed.interp(args),axis=0)
        
