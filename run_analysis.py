import os
import sys
import copy
from fermipy.gtanalysis import GTAnalysis
import numpy as np
import itertools
import argparse

usage = "usage: %(prog)s [config file]"
description = "Run fermipy analysis chain."
parser = argparse.ArgumentParser(usage=usage,description=description)

parser.add_argument('--config', default = 'sample_config.yaml')
parser.add_argument('--source', default = None)

args = parser.parse_args()

gta = GTAnalysis(args.config,logging={'verbosity' : 3})

gta.setup()

gta.update_source_map(args.source)

# Get a reasonable starting point for the spectral model
gta.free_source(args.source)
gta.fit()
gta.free_source(args.source,False)

# -----------------------------------
# Pass 0 - Source at Nominal Position
# -----------------------------------

gta.optimize()

gta.free_source(args.source)
gta.free_source('galdiff',pars='norm')
gta.free_source('isodiff',pars='norm')
gta.fit()
gta.free_sources(free=False)

gta.extension(args.source)
gta.sed(args.source)

model = { 'SpatialModel' : 'PointSource', 'Index' : 2.0 }

# Baseline Fit
gta.write_roi('fit0')
gta.tsmap('fit0',model=model)
gta.tsmap('fit0_nosource',model=model,exclude=[args.source])
gta.residmap('fit0',model=model)

# -------------------------------------
# Pass 1 - Source at Localized Position
# -------------------------------------

if gta.roi[args.source]['ts'] > 25.:
    gta.localize(args.source,update=True)

#gta.optimize()

gta.free_source(args.source)
gta.free_source('galdiff',pars='norm')
gta.free_source('isodiff',pars='norm')
gta.fit()
gta.free_sources(free=False)

gta.extension(args.source)
gta.sed(args.source)

# Post-Relocalization Fit
gta.write_roi('fit1')
gta.tsmap('fit1',model=model)
gta.tsmap('fit1_nosource',model=model,exclude=[args.source])
gta.residmap('fit1',model=model)

skydir = gta.roi.get_source_by_name(args.source)[0].skydir

halo_source_dict = {
    'ra' : skydir.ra.deg,
    'dec' : skydir.dec.deg,
    'SpectrumType' : 'PowerLaw', 
    'Index' : 2.0, 
    'Scale' : 1000,
    'Prefactor' : { 'value' : 1E-5, 'scale' : 1e-13 },
    'SpatialModel' : 'GaussianSource',
    'SpatialWidth' : 1.0
    }


halo_width = np.logspace(-1,0,9)
halo_index = np.array([1.5,2.0,2.5,3.0])

halo_data = []

#for i,w in enumerate(halo_width):
for i, (w,idx) in enumerate(itertools.product(halo_width,halo_index)):
    halo_source_dict['SpatialWidth'] = w
    halo_source_dict['Index'] = idx
    
    halo_source_name = 'halo_gauss'
    
    gta.add_source(halo_source_name,halo_source_dict)
    gta.free_norm(halo_source_name)
    gta.free_norm(args.source)
    gta.free_norm('galdiff')
    gta.free_norm('isodiff')
    gta.fit()

    gta.update_source(halo_source_name,reoptimize=True,
                      npts=20)
    
    gta.write_roi('halo_gauss_%02i'%i,make_plots=False,
                  save_model_map=False,format='npy')
    halo_data += [copy.deepcopy(gta.roi['halo_gauss'])]
    
    gta.delete_source(halo_source_name,save_template=False)    
    gta.load_roi('fit1')

    
np.save(os.path.join(gta._savedir,'halo_data.npy'),halo_data)


for i, w in enumerate(halo_width):
    halo_source_dict['SpatialWidth'] = w
    halo_source_dict['Index'] = 2.0
    
    halo_source_name = 'halo_gauss'
    
    gta.add_source(halo_source_name,halo_source_dict)
    gta.free_norm(halo_source_name)
    gta.free_norm(args.source)
    gta.free_norm('galdiff')
    gta.free_norm('isodiff')
    gta.fit()

#    gta.update_source(halo_source_name,reoptimize=True)
    
    gta.sed(halo_source_name)    
    gta.write_roi('halo_gauss_sed_%02i'%i,make_plots=False,
                  save_model_map=False,format='npy')
    gta.delete_source(halo_source_name,save_template=False)    
    gta.load_roi('fit1')
