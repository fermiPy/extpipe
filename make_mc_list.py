import sys

import itertools
import yaml
from fermipy.utils import *

min_prefactor = 1E-14

testsource = {'name' : 'testsource', 'Index' : 2.0, 'Prefactor' : min_prefactor,
              'Scale' : 3162., 'glon' : 90.0, 'glat' : 50.0, 'SpatialModel' : 'PSFSource'}

config = {
    'selection' : {'glon' : 90.0, 'glat' : 50.0},
    'model' : { 'sources' : [testsource], 'catalogs' : [] },
    'mc' : {'seed' : 0 }
    }


srcs_dict = {}
index = [1.5,2.0,2.5,3.0]
amp_scale = np.linspace(0,2,5)
seeds = np.arange(0,10)
configs = {}

for t, amp, s in itertools.product(index,amp_scale,seeds):

    name = 'powerlaw_testext_%3.1f_%4.2f_%02i'%(t,amp,s)    
    c = copy.deepcopy(config)
    c['model']['sources'][0]['Index'] = t
    c['model']['sources'][0]['Prefactor'] = min_prefactor*10**amp
    c['mc']['seed'] = s
    configs[name] = c
    
yaml.dump(tolist(configs),open(sys.argv[1],'w'))

