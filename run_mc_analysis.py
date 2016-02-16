import os
import sys
import copy
from fermipy.gtanalysis import GTAnalysis
import numpy as np
import itertools
import argparse
import yaml

usage = "usage: %(prog)s [config file]"
description = "Run fermipy analysis chain."
parser = argparse.ArgumentParser(usage=usage,description=description)

parser.add_argument('--config', default = 'sample_config.yaml')
parser.add_argument('--source', default = None)

args = parser.parse_args()

config = yaml.load(open(args.config,'r'))

# Remove point sources from model
#config['model']['catalogs'] = []

gta = GTAnalysis(args.config,logging={'verbosity' : 3})

gta.setup()

for s in gta.roi.point_sources:

    if s.name == 'testsource':
        continue
    
    gta.delete_source(s.name)

ext_fit_data = []

gta.write_roi('base_model',save_model_map=False,make_plots=False)

for i in range(100):

    gta.load_roi('base_model')
    
    gta.simulate_roi()

    gta.free_source('testsource')
    gta.free_source('galdiff',pars='norm')
    gta.free_source('isodiff',pars='norm')
    
    gta.fit()
    gta.free_sources(free=False)

    gta.extension('testsource',width=np.logspace(-2.5,-0.5,9))

    ext_fit_data += [copy.deepcopy(gta.roi['testsource'])]

    gta.write_roi('fit%04i'%i,save_model_map=False,make_plots=False,
                  format='npy')    

np.save(os.path.join(gta._savedir,'ext_fit_data.npy'),ext_fit_data)
    
