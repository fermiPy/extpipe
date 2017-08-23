import yaml
import sys
import numpy as np
from astropy.table import Table
from astropy.coordinates import SkyCoord
from fermipy.catalog import *
from fermipy.utils import *


def get_coord(name,tab):

    row = tab[tab['Source_Name'] == name]
    return SkyCoord(float(row['RAJ2000']), float(row['DEJ2000']),unit='deg')

def avg_coords(coords):
    xyz = np.zeros(3)
    for t in coords:
        xyz += t.cartesian.xyz

    xyz /= np.sum(xyz**2)**0.5
        
    c = SkyCoord(xyz[0], xyz[1], xyz[2],representation='cartesian')
    c.representation='spherical'
    return c



tab = Table.read(sys.argv[1])
src_names = []
    
m = np.abs(tab['glat']) < 0.
#m |= (tab['fit_ext_gauss_ts_ext'] > 9.0)
#m |= (tab['fit_ext_disk_ts_ext'] > 9.0)
m |= (tab['fit_halo_ts'] > 16.0)
#m |= (tab['ts'] > 20000.0)


for row in tab[m]:
    src_names += [row['codename']]

src_names = sorted(list(set(src_names)))
        
o = {}

for name in src_names:
    #coords = [get_coord(t,cat.table) for t in names]
    #c0 = avg_coords(coords)

    print(name)
    #print(create_source_name(c0))
    
    names = [name]
    
    row = tab[tab['codename'] == names[0].lower().replace(' ','_')]    
    c0 = SkyCoord(row['ra'],row['dec'],unit='deg')    
    name = create_source_name(c0).replace('PS','FHES') + 'e'
    #print(c0.ra.deg,c0.dec.deg)
    
    #print(names[0])
    #print(row['codename'])
    
    src = {'name' : name,
           'ra' : float(c0.ra.deg), 'dec' : float(c0.dec.deg),
           'SpectrumType' : 'PowerLaw', 'SpatialModel' : 'RadialGaussian',
           'SpatialWidth' : float(row['fit_halo_r68']),
           'Index' : float(row['fit_halo_index'])}
    
    o[name.lower().replace(' ','_')] = {'selection' : {'target' : name},
               'model' : {'sources' : [src]} }


yaml.dump(o,open('out.yaml','w'))
