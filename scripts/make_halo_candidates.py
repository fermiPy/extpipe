import yaml
import sys
import numpy as np
from astropy.table import Table

tab = Table.read(sys.argv[1])
src_names = []
    
m = np.abs(tab['glat']) < 0.
m |= (tab['fit_ext_gauss_ts_ext'] > 9.0)
m |= (tab['fit_ext_disk_ts_ext'] > 9.0)
m |= (tab['fit_halo_ts'] > 9.0)
#m |= (tab['ts'] > 20000.0)

for row in tab[m]:
    src_names += [row['codename']]

src_names = sorted(list(set(src_names)))
        
o = {}

for key in src_names:

    #key = name.lower().replace(' ','_')

    if 'fhes' in key or 'lmc' in key:
        continue

    if '3fgl' in key:
        name = key.replace('3fgl_j','3FGL J')
    if '3fhl' in key:
        name = key.replace('3fhl_j','3FHL J')
    
    o[str(key)] = {'selection' : {'target' : str(name)} }
    


yaml.dump(o,open('list.yaml','w'))
