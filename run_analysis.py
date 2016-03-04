import os
import sys
import copy
from fermipy.gtanalysis import GTAnalysis
import numpy as np
import itertools
import argparse

os.environ['USE_ADAPTIVE_PSF_ESTIMATOR']='1'

usage = "usage: %(prog)s [config file]"
description = "Run fermipy analysis chain."
parser = argparse.ArgumentParser(usage=usage,description=description)

parser.add_argument('--config', default = 'sample_config.yaml')
parser.add_argument('--source', default = None)

args = parser.parse_args()
gta = GTAnalysis(args.config,logging={'verbosity' : 3})

gta.setup()

model0 = { 'SpatialModel' : 'PointSource', 'Index' : 1.5 }
model1 = { 'SpatialModel' : 'PointSource', 'Index' : 2.0 }
model2 = { 'SpatialModel' : 'PointSource', 'Index' : 2.5 }
src_name = gta.roi.sources[0].name

# reset the source map
gta.reload_source(src_name)

# -----------------------------------
# Get a Baseline Model
# -----------------------------------

# Get a reasonable starting point for the spectral model
gta.free_source(src_name)
gta.fit()
gta.free_source(src_name,False)

gta.optimize()

gta.tsmap('base',model=model1)

# Look for new point sources outside the inner 0.5 deg
gta.find_sources('base',model=model0,
                 search_skydir=gta.roi.skydir,
                 max_iter=3,min_separation=1.0,
                 sqrt_ts_threshold=5,
                 search_minmax_radius=[0.5,None])

gta.find_sources('base',model=model1,
                 search_skydir=gta.roi.skydir,
                 max_iter=3,min_separation=1.0,
                 sqrt_ts_threshold=5,
                 search_minmax_radius=[0.5,None])

gta.find_sources('base',model=model2,
                 search_skydir=gta.roi.skydir,
                 max_iter=3,min_separation=1.0,
                 sqrt_ts_threshold=5,
                 search_minmax_radius=[0.5,None])

gta.optimize()

gta.write_roi('base')

# -----------------------------------
# Pass 0 - Source at Nominal Position
# -----------------------------------

gta.free_sources(distance=0.5,exclude_diffuse=True)
gta.free_sources(distance=0.5,pars='norm')
gta.fit()
#gta.free_sources(free=False)

gta.extension(src_name)
gta.sed(src_name)

# Baseline Fit
gta.write_roi('fit0')
gta.tsmap('fit0',model=model1)
gta.tsmap('fit0_nosource',model=model1,exclude=[src_name])
gta.residmap('fit0',model=model1)

# Likelihood for Model
like0 = -gta.like()

gta.logger.info('Fit0 Model Likelihood: %f'%like0)

# -------------------------------------
# Pass 1 - Source at Localized Position
# -------------------------------------

if gta.roi[src_name]['ts'] > 9.:
    gta.localize(src_name,nstep=5,dtheta_max=0.5,update=True)

gta.optimize()

gta.free_sources(distance=0.5,exclude_diffuse=True)
gta.free_sources(distance=0.5,pars='norm')
gta.fit()
#gta.free_sources(free=False)

gta.extension(src_name)
gta.sed(src_name)

# Post-Relocalization Fit
gta.write_roi('fit1')
gta.tsmap('fit1',model=model0)
gta.tsmap('fit1',model=model1)
gta.tsmap('fit1',model=model2)
gta.tsmap('fit1_nosource',model=model1,exclude=[src_name])
gta.residmap('fit1',model=model1)

# Likelihood for Model
like1 = -gta.like()

gta.logger.info('Fit1 Model Likelihood: %f'%like1)

# -------------------------------------
# Pass 2 - Multiple Point Sources
# -------------------------------------

# Delete the central source
src1 = gta.delete_source(src_name)

src_dict = copy.deepcopy(src1.data)

if 'extension' in src_dict:
    del src_dict['extension']

if 'sed' in src_dict:
    del src_dict['sed']

# Search for up to 2 sources in a circle of 0.5 deg radius
srcs_pass0 = gta.find_sources('fit2_pass0',model=src_dict,
                         search_skydir=gta.roi.skydir,
                         max_iter=1,
                         sources_per_iter=1,
                         sqrt_ts_threshold=3,
                         min_separation=0.5,
                         search_minmax_radius=[None,0.5])

srcs_pass1 = gta.find_sources('fit2_pass1',
                         search_skydir=gta.roi.skydir,
                         max_iter=1,
                         sources_per_iter=1,
                         sqrt_ts_threshold=3,
                         min_separation=0.5,
                         search_minmax_radius=[None,0.5])

srcs = srcs_pass0['sources'] + srcs_pass1['sources']

for s in srcs:
    gta.localize(s.name,nstep=5,dtheta_max=0.4,
                 update=True)

# If no sources were found then re-add the central source
if len(srcs) == 0:
    gta.add_source(src1.name,src1)

gta.optimize()

gta.free_sources(distance=0.5,exclude_diffuse=True)
gta.free_sources(distance=0.5,pars='norm')
gta.fit()

src_name = gta.roi.sources[0].name

for s in srcs:
    gta.extension(s.name)
    gta.sed(s.name)

new_source_data = []
for s in srcs:
    src_data = gta.roi[s.name].data
    new_source_data.append(copy.deepcopy(src_data))

np.save(os.path.join(gta._savedir,'new_source_data.npy'),
        new_source_data)
        
    
gta.write_roi('fit2')
gta.tsmap('fit2',model=model0)
gta.tsmap('fit2',model=model1)
gta.tsmap('fit2',model=model2)
gta.tsmap('fit2_nosource',model=model1,exclude=[src_name])
gta.residmap('fit2',model=model1)

like2 = -gta.like()

gta.logger.info('Fit2 Model Likelihood: %f'%like2)

halo_source_dict = {
    'SpectrumType' : 'PowerLaw', 
    'Index' : 2.0, 
    'Scale' : 1000,
    'Prefactor' : { 'value' : 1E-5, 'scale' : 1e-13 },
    'SpatialModel' : 'GaussianSource',
    'SpatialWidth' : 1.0
    }


halo_width = np.logspace(-1,0,9)
halo_index = np.array([1.5,2.0,2.5,3.0])

fit1_halo_data = []
fit2_halo_data = []

for i, (w,idx) in enumerate(itertools.product(halo_width,halo_index)):
    halo_source_dict['SpatialWidth'] = w
    halo_source_dict['Index'] = idx    
    halo_source_name = 'halo_gauss'

    # Halo Search with Single Source Model
    gta.load_roi('fit1')

    halo_source_dict['ra'] = gta.roi.sources[0]['ra']
    halo_source_dict['dec'] = gta.roi.sources[0]['dec']
    
    gta.add_source(halo_source_name,halo_source_dict)
    gta.free_sources(distance=1.0,pars='norm')
    gta.fit()
    
    gta.update_source(halo_source_name,reoptimize=True,
                      npts=10)

    gta.logger.info('Fit1 Halo Width: %6.3f Index: %6.2f TS: %6.2f'%(w,idx,gta.roi[halo_source_name]['ts']))
    
    gta.write_roi('fit1_halo_gauss_%02i'%i,make_plots=False,
                  save_model_map=False,format='npy')
    fit1_halo_data += [copy.deepcopy(gta.roi['halo_gauss'].data)]    
    gta.delete_source(halo_source_name,save_template=False)    

    # Halo Search with Multiple Source Model
    gta.load_roi('fit2')

    halo_source_dict['ra'] = gta.roi.sources[0]['ra']
    halo_source_dict['dec'] = gta.roi.sources[0]['dec']
    
    gta.add_source(halo_source_name,halo_source_dict)
    gta.free_sources(distance=1.0,pars='norm')
    gta.fit()

    gta.update_source(halo_source_name,reoptimize=True,
                      npts=10)

    gta.logger.info('Fit2 Halo Width: %6.3f Index: %6.2f TS: %6.2f'%(w,idx,gta.roi[halo_source_name]['ts']))
    
    gta.write_roi('fit2_halo_gauss_%02i'%i,make_plots=False,
                  save_model_map=False,format='npy')
    fit2_halo_data += [copy.deepcopy(gta.roi['halo_gauss'].data)]                        
    gta.delete_source(halo_source_name,save_template=False)
    
    
np.save(os.path.join(gta._savedir,'fit1_halo_data.npy'),fit1_halo_data)
np.save(os.path.join(gta._savedir,'fit2_halo_data.npy'),fit2_halo_data)

for i, w in enumerate(halo_width):

    gta.load_roi('fit1')
    
    halo_source_dict['SpatialWidth'] = w
    halo_source_dict['Index'] = 2.0    
    halo_source_name = 'halo_gauss'

    halo_source_dict['ra'] = gta.roi.sources[0]['ra']
    halo_source_dict['dec'] = gta.roi.sources[0]['dec']
    
    gta.add_source(halo_source_name,halo_source_dict)
    gta.free_sources(distance=1.0,pars='norm')
    gta.fit()

#    gta.update_source(halo_source_name,reoptimize=True)
    
    gta.sed(halo_source_name)    
    gta.write_roi('fit1_halo_gauss_sed_%02i'%i,make_plots=False,
                  save_model_map=False,format='npy')
    gta.delete_source(halo_source_name,save_template=False) 
