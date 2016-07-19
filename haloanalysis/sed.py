import numpy as np

from haloanalysis.utils import Axis, MapND


class HaloSED(object):
    """Representation of a halo SED containing log-likelihood values as a
    function of energy and angular size.

    """

    def __init__(self, sed):

        self._sed = sed

    @staticmethod
    def create_from_fits(row):

        axis0 = Axis('eobs',np.linspace(3,5.5,21))
        axis1 = Axis.create_from_centers('width',np.linspace(-1.125, 1.125, 13))
        axis2 = Axis.create_from_centers('eflux',np.linspace(-9, -5, 41))

        sed = MapND([axis1,axis0,axis2],-np.array(row['fit_halo_sed_scan_dlnl']))

        return HaloSED(sed)
        
    def __call__(self, flux, th68):
        return self.nll(flux, th68)

    def nll(self, flux, width):
        """Evaluate the negative log-likelihood function."""

        ectr = 10**self._sed.axes[1].centers        
        log_eflux = np.log10(flux*ectr)
        log_width = np.log10(width)
        args = (log_width, self._sed.axes[1].centers, log_eflux)
        
        return np.sum(self._sed.interp(args))
        
