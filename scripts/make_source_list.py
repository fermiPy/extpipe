import yaml
import sys
from fermipy.utils import *

from fermipy.roi_model import ROIModel

roi = ROIModel({'catalogs' : ['/u/gl/mdwood/fermi/catalogs/gll_psc_v16_ext.fit','lmc_psc_v0.fit'],
                'extdir' : '/u/gl/mdwood/fermi/catalogs/Extended_archive_v16'})

#cuts0 = [('TS_value',10000.0,None),('GLAT',5.,90.)]
#cuts1 = [('TS_value',10000.0,None),('GLAT',-90.,-5.)]

cuts0 = [('GLAT',5.,90.)]
cuts1 = [('GLAT',-90.,-5.)]

srcs = roi.get_sources(cuts=cuts0) + roi.get_sources(cuts=cuts1)

config = {
    'selection' : {},
    }

configs = {}
for s in srcs:

    if s.extended:
        print(s['glat'],s.name)
#        continue

    c = copy.deepcopy(config)
    c['selection']['target'] = s.name
    configs[s.name.lower().replace(' ','_')] = c
    
#    srcs_dict[s.name] = s
    
print 'found ', len(configs), ' sources'

#src_names = []
#for k,s in sorted(srcs_dict.items()):
#    src_names += [k]

yaml.dump(tolist(configs),open(sys.argv[1],'w'))

