import yaml
from fermipy.utils import *

from fermipy.roi_model import ROIModel




roi = ROIModel({'catalogs' : ['3FGL']})

#cuts0 = [('TS_value',10000.0,None),('GLAT',5.,90.)]
#cuts1 = [('TS_value',10000.0,None),('GLAT',-90.,-5.)]

cuts0 = [('GLAT',5.,90.)]
cuts1 = [('GLAT',-90.,-5.)]

srcs = roi.get_sources(cuts0) + roi.get_sources(cuts1)

srcs_dict = {}

for s in srcs:

    if s.extended:
        continue
    srcs_dict[s.name] = s
    
print 'found ', len(srcs_dict), ' sources'

src_names = []

for k,s in sorted(srcs_dict.items()):
    src_names += [k]

yaml.dump(tolist(src_names),open(sys.argv[1],'w'))

