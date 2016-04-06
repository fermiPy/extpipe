import os
import sys
import copy
from fermipy.gtanalysis import GTAnalysis
import numpy as np
import itertools
import argparse
from haloanalysis.fit_funcs import fit_region, fit_halo
    



if __name__ == '__main__':
        
    usage = "usage: %(prog)s [config file]"
    description = "Run fermipy analysis chain."
    parser = argparse.ArgumentParser(usage=usage,description=description)

    parser.add_argument('--config', default = 'sample_config.yaml')
    parser.add_argument('--source', default = None)

    args = parser.parse_args()
    gta = GTAnalysis(args.config,logging={'verbosity' : 3})

    gta.setup(overwrite=True)

    sqrt_ts_threshold=3

    halo_width = np.logspace(-1,0,9)
    halo_index = np.array([1.5,1.75,2.0,2.25,2.5,2.75,3.0])
    
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

    fit_region(gta,'fit0',src_name)
    fit_region(gta,'fit0_emin40',src_name,erange=[4.0,5.5])

    gta.load_roi('fit0')

    # -------------------------------------
    # Pass 1 - Source at Localized Position
    # -------------------------------------

    gta.localize(src_name,nstep=5,dtheta_max=0.5,update=True,
                 prefix='fit1')

    fit_region(gta,'fit1',src_name)
    fit_halo(gta,'fit1',src_name,halo_width,halo_index)

    fit_region(gta,'fit1_emin40',src_name,erange=[4.0,5.5])
    fit_halo(gta,'fit1_emin40',src_name,halo_width,halo_index,erange=[4.0,5.5])

    gta.load_roi('fit1')

    # -------------------------------------
    # Pass 2 - 2+ Point Sources
    # -------------------------------------

    srcs = []
    best_fit_idx = 1

    # Fit up to 4 sources
    for i in range(2,5):

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
        best_fit_idx = i
        
        gta.localize(src_name,nstep=5,dtheta_max=0.4,
                     update=True,prefix='fit%i'%i)

        # Relocalize new sources
        for s in sorted(srcs, key=lambda t: t['ts'],reverse=True):        
            gta.localize(s.name,nstep=5,dtheta_max=0.4,
                         update=True,prefix='fit%i'%i)

        fit_region(gta,'fit%i'%i,src_name)
        fit_halo(gta,'fit%i'%i,src_name,halo_width,halo_index)

        fit_region(gta,'fit%i_emin40'%i,src_name,erange=[4.0,5.5])
        fit_halo(gta,'fit%i_emin40'%i,src_name,halo_width,halo_index,erange=[4.0,5.5])

        gta.load_roi('fit%i'%i)

    # Only Run Halo Fit for Best-fit Model
    #if best_fit_idx > 1:
    #    fit_halo(gta,'fit%i'%best_fit_idx,src_name)
    #    fit_halo(gta,'fit%i_emin40'%best_fit_idx,src_name,erange=[4.0,5.5])
        
    new_source_data = []
    for s in srcs:
        src_data = gta.roi[s.name].data
        new_source_data.append(copy.deepcopy(src_data))

    np.save(os.path.join(gta.workdir,'new_source_data.npy'),
            new_source_data)
