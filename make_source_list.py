import yaml
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

cuts0 = [('TS_value',10000.0,None),('GLAT',5.,90.)]
cuts1 = [('TS_value',10000.0,None),('GLAT',-90.,-5.)]



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

