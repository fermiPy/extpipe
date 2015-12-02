import yaml
import numpy as np
from fermipy.utils import *

from fermipy.roi_model import ROIModel

import sys
import time, os, stat







config = {}
config['model'] = dict(
    src_roiwidth = 10.0,
    galdiff  = '/nfs/slac/g/ki/ki20/cta/mdwood/fermi/diffuse/v5r0/gll_iem_v06.fits',
    isodiff  = '/nfs/slac/g/ki/ki20/cta/mdwood/fermi/diffuse/v5r0/iso_P8R2_SOURCE_V6_v06.txt',
    catalogs = ['gll_psc_v14.fit'])

config['run'] = {}

roi = ROIModel(config['model'])
roi.load()


o = {'ts_ext' : [], 'ext_ul95' : [], 'ext' : [], 'name' : [], 'assoc' : []}

for d in sys.argv[1:]:

    
#    results = yaml.load(open(os.path.join(d,'fit1.yaml')),Loader=yaml.CLoader)

    results_file = os.path.join(d,'fit1.npy')

    if not os.path.isfile(results_file):
        print 'Skipping ', d
        continue
    
    results = np.load(results_file)
    results = np.array(results,ndmin=1)[0]
    
    srcname = results['config']['run']['sed'][0]

    src = roi.get_source_by_name(srcname)


    ext = results['sources'][srcname]['extension']
    
    o['ts_ext'] += [ext['ts_ext']]
    o['ext_ul95'] += [ext['ext_ul95']]
    o['ext'] += [ext['ext']]
    o['name'] += [str(srcname)]

    if src['ASSOC1']:    
        o['assoc'] += [str(src['ASSOC1'])]
    else:
        o['assoc'] += ['']
        
    print '%30s %30s %10.2f %10.2f %10.2f'%(srcname.lower().replace(' ','_'),
                                            src['ASSOC1'],
                                            src['TS_value'],
                                            ext['ts_ext'],
                                            ext['ext'])
    
#    print d, srcname                       
#print o
yaml.dump(o,open('results.yaml','w'))

sys.exit(0)

cuts0 = [('TS_value',100.0,None),('GLAT',5.,90.)]
cuts1 = [('TS_value',100.0,None),('GLAT',-90.,-5.)]



srcs = roi.get_sources(cuts0) + roi.get_sources(cuts1)

srcs = [roi.get_source_by_name('vela')] + srcs
srcs = [roi.get_source_by_name('geminga')] + srcs
srcs = [roi.get_source_by_name('crab')] + srcs

srcs_dict = {}

for s in srcs:    
    srcs_dict[s.name] = s
    
print 'found ', len(srcs), ' sources'

src_names = []

for k,s in sorted(srcs_dict.items()):
    src_names += [k]

yaml.dump(src_names,open(sys.argv[1],'w'))

