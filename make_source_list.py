import yaml
import sys
from fermipy.utils import *

from fermipy.roi_model import ROIModel

roi = ROIModel({'catalogs' : ['3FGL']})

#cuts0 = [('TS_value',10000.0,None),('GLAT',5.,90.)]
#cuts1 = [('TS_value',10000.0,None),('GLAT',-90.,-5.)]

cuts0 = [('GLAT',5.,90.)]
cuts1 = [('GLAT',-90.,-5.)]

srcs = roi.get_sources(cuts0) + roi.get_sources(cuts1)

config = {
    'selection' : {},
    }

configs = {}
for s in srcs:

    if s.extended:
        continue

    c = copy.deepcopy(config)
    c['selection']['target'] = s.name
    configs[s.name.lower().replace(' ','_')] = c
    
#    srcs_dict[s.name] = s
    
print 'found ', len(configs), ' sources'

#src_names = []
#for k,s in sorted(srcs_dict.items()):
#    src_names += [k]

yaml.dump(tolist(configs),open(sys.argv[1],'w'))

