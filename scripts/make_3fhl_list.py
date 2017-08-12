import yaml
import sys
from astropy.table import Table
from astropy.coordinates import SkyCoord
from fermipy.utils import *



tab0 = Table.read('/u/gl/mdwood/fermi/catalogs/gll_psch_v11.fit')
cat_3fgl = Table.read('/u/gl/mdwood/fermi/catalogs/gll_psc_v16_ext.fit')

skydir_3fgl = SkyCoord(cat_3fgl['RAJ2000'],cat_3fgl['DEJ2000'],unit='deg')

m0 = np.zeros(len(tab0),dtype=bool)
for i, row in enumerate(tab0):
    if '3FGL' in row['ASSOC_GAM']:
        m0[i] = True

m1 = np.abs(tab0['GLAT']) > 5.0

tab0 = tab0[m1 & ~m0]

#3FHL J0531.8-6639e

config = {
    'selection' : {},
    'model' :
        { 'sources' : [{'SpectrumType' : 'PowerLaw', 'Index' : 2.0, 'Prefactor' : 1E-13, 'Scale' : 1E3}]

        },
    }

configs = {}
for row in tab0:    

    skydir = SkyCoord(row['RAJ2000'],row['DEJ2000'],unit='deg')
    sep = skydir.separation(skydir_3fgl).deg
    min_sep = np.min(sep)
    
    idx = np.argmin(sep)
    row_3fgl = cat_3fgl[idx]
    sigma95 = np.sqrt(np.array(row_3fgl['Conf_95_SemiMajor'])*np.array(row_3fgl['Conf_95_SemiMinor']))
    
    name = row['Source_Name'].strip()

    if min_sep/sigma95 < 2.0:
        print('skipping',name,min_sep,sigma95,min_sep/sigma95)    
        continue
    
    c = copy.deepcopy(config)
    c['selection']['target'] = name
    c['model']['sources'][0]['name'] = name
    c['model']['sources'][0]['ra'] = row['RAJ2000']
    c['model']['sources'][0]['dec'] = row['DEJ2000']
    c['model']['sources'][0]['Index'] = min(3.0,max(1.5,row['Spectral_Index']))

    configs[name.lower().replace(' ','_')] = c
    
#    srcs_dict[s.name] = s
    
print 'found ', len(configs), ' sources'

#src_names = []
#for k,s in sorted(srcs_dict.items()):
#    src_names += [k]

yaml.dump(tolist(configs),open(sys.argv[1],'w'))

