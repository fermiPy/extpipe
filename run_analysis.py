import os
import sys
import copy
from fermipy.gtanalysis import GTAnalysis
import numpy as np
import itertools
import argparse


def fit_region(gta,modelname,erange=None):

    if erange is not None:
        gta.setEnergyRange(erange[0],erange[1])
    
    model0 = { 'SpatialModel' : 'PointSource', 'Index' : 1.5 }
    model1 = { 'SpatialModel' : 'PointSource', 'Index' : 2.0 }
    model2 = { 'SpatialModel' : 'PointSource', 'Index' : 2.5 }
    
    gta.optimize()

    gta.free_sources(False)
    gta.free_sources(distance=1.0,exclude_diffuse=True)
    gta.free_sources(distance=1.0,pars='norm')
    gta.fit()

    src_name = gta.roi.sources[0].name

    for s in gta.roi.sources:

        if s['offset'] > 1.0:
            continue

        if not s['SpatialModel'] == 'PointSource':
            continue
        
        gta.extension(s.name)
        gta.sed(s.name)
    
    gta.write_roi(modelname)
    gta.tsmap(modelname,model=model0,erange=erange)
    gta.tsmap(modelname,model=model1,erange=erange)
    gta.tsmap(modelname,model=model2,erange=erange)
    gta.tsmap('%s_nosource'%modelname,
              model=model1,exclude=[src_name],erange=erange)
    gta.residmap(modelname,model=model1,erange=erange)

    lnl = -gta.like()

    gta.logger.info('%s Model Likelihood: %f'%(modelname,lnl))

def fit_halo(gta,modelname,src_name,erange=None):

    halo_source_dict = {
        'SpectrumType' : 'PowerLaw', 
        'Index' : 2.0, 
        'Scale' : 1000,
        'Prefactor' : { 'value' : 1E-5, 'scale' : 1e-13 },
        'SpatialModel' : 'GaussianSource',
        'SpatialWidth' : 1.0
        }
    
    halo_width = np.logspace(-1,0,9)
    halo_index = np.array([1.5,1.75,2.0,2.25,2.5,2.75,3.0])

    halo_data = []

    gta.load_roi(modelname)
    if erange is not None:
        gta.setEnergyRange(erange[0],erange[1])

    gta.free_sources(False)
    gta.free_sources(distance=1.0,pars='norm')
    gta.write_xml(modelname + '_base')
    
    for i, (w,idx) in enumerate(itertools.product(halo_width,halo_index)):
        halo_source_dict['SpatialWidth'] = w
        halo_source_dict['Index'] = idx    
        halo_source_name = 'halo_gauss'

        gta.load_xml(modelname + '_base')

        halo_source_dict['ra'] = gta.roi[src_name]['ra']
        halo_source_dict['dec'] = gta.roi[src_name]['dec']
    
        gta.add_source(halo_source_name,halo_source_dict,free=True)
#        gta.free_sources(False)
#        gta.free_sources(distance=1.0,pars='norm')
        gta.fit(update=False)
    
        gta.update_source(halo_source_name,reoptimize=True,
                          npts=10)

        gta.logger.info('%s Halo Width: %6.3f Index: %6.2f TS: %6.2f'%(modelname,w,idx,
                                                                       gta.roi[halo_source_name]['ts']))
    
        gta.write_roi('%s_halo_gauss_%02i'%(modelname,i),make_plots=False,
                      save_model_map=False,format='npy')
        halo_data += [copy.deepcopy(gta.roi['halo_gauss'].data)]    
        gta.delete_source(halo_source_name,save_template=False) 


    np.save(os.path.join(gta.workdir,'%s_halo_data.npy'%modelname),halo_data)

    #gta.load_roi(modelname)
    if erange is not None:
        gta.setEnergyRange(erange[0],erange[1])

    for i, w in enumerate(halo_width):        
        halo_source_dict['SpatialWidth'] = w
        halo_source_dict['Index'] = 2.0    
        halo_source_name = 'halo_gauss'

        gta.load_xml(modelname + '_base')
        
        halo_source_dict['ra'] = gta.roi[src_name]['ra']
        halo_source_dict['dec'] = gta.roi[src_name]['dec']
        
        gta.add_source(halo_source_name,halo_source_dict,free=True)
#        gta.free_sources(distance=1.0,pars='norm')
        gta.fit()
    
        gta.sed(halo_source_name)    
        gta.write_roi('%s_halo_gauss_sed_%02i'%(modelname,i),make_plots=False,
                      save_model_map=False,format='npy')
        gta.delete_source(halo_source_name,save_template=False) 


usage = "usage: %(prog)s [config file]"
description = "Run fermipy analysis chain."
parser = argparse.ArgumentParser(usage=usage,description=description)

parser.add_argument('--config', default = 'sample_config.yaml')
parser.add_argument('--source', default = None)

args = parser.parse_args()
gta = GTAnalysis(args.config,logging={'verbosity' : 3})

gta.setup(overwrite=True)

sqrt_ts_threshold=3

model0 = { 'SpatialModel' : 'PointSource', 'Index' : 1.5 }
model1 = { 'SpatialModel' : 'PointSource', 'Index' : 2.0 }
model2 = { 'SpatialModel' : 'PointSource', 'Index' : 2.5 }
src_name = gta.roi.sources[0].name

# -----------------------------------
# Get a Baseline Model
# -----------------------------------

# Get a reasonable starting point for the spectral model
gta.free_source(src_name)
gta.fit()
gta.free_source(src_name,False)

gta.optimize()

# Localize 3FGL sources
for s in gta.roi.sources:

    if not s['SpatialModel'] == 'PointSource':
        continue

    if s['offset'] < 0.5 or s['ts'] < 25.:
        continue

    if np.abs(s['offset_glon']) > 2.8 or np.abs(s['offset_glat']) > 2.8:
        continue

    gta.localize(s.name,nstep=5,dtheta_max=0.5,update=True,
                 prefix='base')

    gta.free_source(s.name,False)
        
gta.tsmap('base',model=model1)
gta.tsmap('base_emin40',model=model1,erange=[4.0,5.5])

# Look for new point sources outside the inner 1.0 deg

gta.find_sources('base',model=model1,
                 search_skydir=gta.roi.skydir,
                 max_iter=4,min_separation=0.5,
                 sqrt_ts_threshold=sqrt_ts_threshold,
                 search_minmax_radius=[1.0,None])
gta.optimize()

gta.write_roi('base')

# -----------------------------------
# Pass 0 - Source at Nominal Position
# -----------------------------------

fit_region(gta,'fit0')
fit_region(gta,'fit0_emin40',erange=[4.0,5.5])

gta.load_roi('fit0')

# -------------------------------------
# Pass 1 - Source at Localized Position
# -------------------------------------

gta.localize(src_name,nstep=5,dtheta_max=0.5,update=True,
             prefix='fit1')

fit_region(gta,'fit1')
fit_halo(gta,'fit1',src_name)

fit_region(gta,'fit1_emin40',erange=[4.0,5.5])
fit_halo(gta,'fit1_emin40',src_name,erange=[4.0,5.5])

gta.load_roi('fit1')

# -------------------------------------
# Pass 2 - 2+ Point Sources
# -------------------------------------

srcs = []

for i in range(2,6):

    srcs_fit = gta.find_sources('fit%i'%i,
                                search_skydir=gta.roi.skydir,
                                max_iter=1,
                                sources_per_iter=1,
                                sqrt_ts_threshold=3,
                                min_separation=0.5,
                                search_minmax_radius=[None,1.0])

    if len(srcs_fit['sources']) == 0:
        break

    srcs += srcs_fit['sources']

    gta.localize(src_name,nstep=5,dtheta_max=0.4,
                 update=True,prefix='fit%i'%i)
    
    # Relocalize new sources
    for s in sorted(srcs, key=lambda t: t['ts'],reverse=True):        
        gta.localize(s.name,nstep=5,dtheta_max=0.4,
                     update=True,prefix='fit%i'%i)

    fit_region(gta,'fit%i'%i)
    fit_halo(gta,'fit%i'%i,src_name)

    fit_region(gta,'fit%i_emin40'%i,erange=[4.0,5.5])
    fit_halo(gta,'fit%i_emin40'%i,src_name,erange=[4.0,5.5])
    
    gta.load_roi('fit%i'%i)


new_source_data = []
for s in srcs:
    src_data = gta.roi[s.name].data
    new_source_data.append(copy.deepcopy(src_data))

np.save(os.path.join(gta.workdir,'new_source_data.npy'),
        new_source_data)
