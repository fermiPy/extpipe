import numpy as np
from fermipy.skymap import Map
from fermipy.wcs_utils import create_wcs

from astropy.io import fits
from astropy.table import Table,Column

m = Map.create_from_fits('../../diffuse/gll_iem_v06.fits')

h = fits.open('../../diffuse/gll_iem_v06.fits')

nbin = h[0].data.shape[0]
data = np.zeros((nbin+5,1441,2880),dtype=np.float32)

loge_ctr = np.log(m._ectr)
dloge = loge_ctr[-1]-loge_ctr[-2]

loge_ctr = np.concatenate((loge_ctr,
                           loge_ctr[-1] + np.linspace(1,5,5)*dloge))

logm0 = np.log(h[0].data[-1,...])
logm1 = np.log(h[0].data[-2,...])
dydx = (logm0-logm1)/dloge

data[:nbin,...] = h[0].data
for i in range(5):
    data[nbin+i,...] = np.exp(logm0 + dydx*dloge*float(i+1))

tab = Table([Column(name='Energy',data=np.exp(loge_ctr))])    
h[1] = fits.table_to_hdu(tab)
h[1].name = 'ENERGIES'
h[0].data = data
h.writeto('test.fits',overwrite=True)
