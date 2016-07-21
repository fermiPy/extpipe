import numpy as np

from haloanalysis.utils import Axis, MapND


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
        log_eflux = np.log10(flux*ectr)
        log_width = np.log10(width)


        #log_eflux[log_eflux < self.axes[2].lo[0]]
        #log_width[log_width < -1.0] = -1.0
        
        args = (log_width, self._sed.axes[1].centers, log_eflux)
        
        return np.sum(self._sed.interp(args))
        
