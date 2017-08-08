import numpy as np
import os
import sys
from fermipy.skymap import Map
from fermipy.wcs_utils import create_wcs

from astropy.io import fits
from astropy.table import Table,Column

outfile = os.path.splitext(sys.argv[1])[0] + '_ext.txt'

iso = np.loadtxt(sys.argv[1],unpack=True)

x = np.log(iso[0])
y = np.log(iso[1])
yerr = iso[2]
dx = (x[-1] - x[-2])
dydx = (y[-1] - y[-2])/dx

nbin = iso.shape[1]

shape = list(iso.shape)
shape[1] += 4
iso_new = np.zeros(shape)

iso_new[:,:nbin] = iso
iso_new[0,nbin:] = np.exp(x[-1] + dx*np.linspace(1,4,4))
iso_new[1,nbin:] = np.exp(y[-1] + dydx*np.linspace(1,4,4)*dx)


np.savetxt(outfile,iso_new.T)
