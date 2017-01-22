import os
import sys
import copy
import numpy as np
import itertools
import argparse

from fermipy.utils import init_matplotlib_backend

init_matplotlib_backend()

from fermipy.gtanalysis import GTAnalysis
from haloanalysis.fit_funcs import fit_region, fit_halo

def main():
        
    usage = "usage: %(prog)s [config file]"
    description = "Run fermipy analysis chain."
    parser = argparse.ArgumentParser(usage=usage,description=description)

    parser.add_argument('--config', default = 'sample_config.yaml')
    parser.add_argument('--source', default = None)

    args = parser.parse_args()
    gta = GTAnalysis(args.config,logging={'verbosity' : 3},
                     fileio={'workdir_regex' : '\.xml$|\.npy$'})

    if args.source is None:
        src_name = gta.config['selection']['target']
    else:
        src_name = args.source
    
    # Delete unassociated sources
    for s in gta.roi.sources:
        if s.name == src_name:
            continue        
        if s['class'] == '' and not s.diffuse:
            gta.delete_source(s.name)
    
    gta.setup()

    names = [s.name for s in gta.roi.sources if not s.diffuse]
    gta.reload_sources(names)

    gta.print_roi()
    
    sqrt_ts_threshold=3
    
    model0 = { 'SpatialModel' : 'PointSource', 'Index' : 1.5 }
    model1 = { 'SpatialModel' : 'PointSource', 'Index' : 2.0 }
    model2 = { 'SpatialModel' : 'PointSource', 'Index' : 2.5 }

    newsrc_model = { 'SpectrumType' : 'LogParabola',
                     'alpha' : 2.0,
                     'beta' : {'value' : 0.0, 'min' : 0.0, 'max' : 5.0} }
    
    # -----------------------------------
    # Fit the Baseline Model
    # -----------------------------------

    gta.optimize()

    gta.print_roi()
    
    # Localize all point sources
    for s in sorted(gta.roi.sources, key=lambda t: t['ts'],reverse=True):
#    for s in gta.roi.sources:

        if not s['SpatialModel'] == 'PointSource':
            continue

#        if s['offset'] < 0.5 or s['ts'] < 25.:
#            continue

        if s['offset_roi_edge'] > -0.1:
            continue
        
        gta.localize(s.name,nstep=5,dtheta_max=0.5,update=True,
                     prefix='base', make_plots=True)

        gta.free_source(s.name,False)

    gta.tsmap('base_nosrcs',model=model1, make_plots=True)

    # Look for new point sources outside the inner 1.0 deg

    gta.find_sources('base_pass0',model=newsrc_model,
                     search_skydir=gta.roi.skydir,
                     max_iter=5,min_separation=0.5,
                     sqrt_ts_threshold=sqrt_ts_threshold,
                     search_minmax_radius=[1.0,None],
                     free_params=['alpha','norm'])
    gta.optimize()

    gta.find_sources('base_pass1',model=newsrc_model,
                     search_skydir=gta.roi.skydir,
                     max_iter=5,min_separation=0.5,
                     sqrt_ts_threshold=sqrt_ts_threshold,
                     search_minmax_radius=[1.0,None],
                     free_params=['alpha','norm'])

    
    gta.print_roi()

    gta.tsmap('base',model=model1, make_plots=True)

    gta.write_roi('base')

    # -------------------------------------
    # Pass 1 - Source at Localized Position
    # -------------------------------------

    fit_region(gta,'fit0',src_name)
    #fit_halo(gta,'fit0',src_name)
    gta.load_roi('fit0')

    # -------------------------------------
    # Pass 2 - 2+ Point Sources
    # -------------------------------------

    srcs = []

    # Fit up to 4 sources
    for i in range(1,5):

        srcs_fit = gta.find_sources('fit%i'%i,
                                    model=newsrc_model,
                                    search_skydir=gta.roi.skydir,
                                    max_iter=1,
                                    sources_per_iter=1,
                                    sqrt_ts_threshold=3,
                                    min_separation=0.5,
                                    search_minmax_radius=[None,1.0],
                                    free_params=['alpha','norm'])

        if len(srcs_fit['sources']) == 0:
            break

        srcs += srcs_fit['sources']
        
        gta.localize(src_name,nstep=5,dtheta_max=0.4,
                     update=True,prefix='fit%i'%i,
                     free_radius=0.5, make_plots=True)

        # Relocalize new sources
        for s in sorted(srcs, key=lambda t: t['ts'],reverse=True):        
            gta.localize(s.name,nstep=5,dtheta_max=0.4,
                         update=True,prefix='fit%i'%i,
                         free_radius=0.5)

        fit_region(gta,'fit%i'%i,src_name)
        #fit_halo(gta,'fit%i'%i,src_name)

        gta.load_roi('fit%i'%i)
        
    new_source_data = []
    for s in srcs:
        src_data = gta.roi[s.name].data
        new_source_data.append(copy.deepcopy(src_data))

    np.save(os.path.join(gta.workdir,'new_source_data.npy'),
            new_source_data)


if __name__ == '__main__':

    main()
