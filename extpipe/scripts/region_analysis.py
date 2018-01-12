import os
import sys
import copy
import numpy as np
import itertools
import argparse
import logging

from fermipy.utils import init_matplotlib_backend

init_matplotlib_backend()

from fermipy.gtanalysis import GTAnalysis
from fermipy.catalog import Catalog3FGL
from haloanalysis.fit_funcs import fit_region, fit_halo

def localize(gta, spatial_models, prefix, skip_loc, src_name):

    for src in sorted(gta.roi.sources, key=lambda t: t['ts'],reverse=True):

        if not src['SpatialModel'] in spatial_models:
            continue

        if src['offset_roi_edge'] > -0.1:
            continue

        if src.name in skip_loc:
            continue
        
        o = gta.localize(src.name,nstep=5,
                         dtheta_max=max(0.5,src['SpatialWidth']),
                         fix_shape=False if src['ts'] >= 25.0 else True,
                         update=True,
                         prefix=prefix, make_plots=True)

        if src.name == src_name and ((not o['fit_success']) or (not o['fit_inbounds'])):
            gta.localize(src.name,nstep=5,
                         dtheta_max=max(1.0,src['SpatialWidth']),
                         fix_shape=False if src['ts'] >= 25.0 else True,
                         update=True,
                         prefix=prefix, make_plots=True)

def optimize_source(gta, name):
    """Optimize and individual source while keeping all other
    sources in the model fixed."""
        
    # preserve free parameters
    free = gta.get_free_param_vector()

    # Fix all parameters
    gta.free_sources(free=False, loglevel=logging.DEBUG)
    gta.free_source(name)
    gta.fit()
    gta.set_free_param_vector(free)

def update_to_lp(gta, ts_thresh, names=None, exclude=None):

    for src in sorted(gta.roi.sources, key=lambda t: t['ts'],reverse=True):

        if names is not None and src.name not in names:
            continue

        if exclude is not None and src.name in exclude:
            continue
        
        if src.diffuse or src['SpectrumType'] != 'PowerLaw':
            continue
                
        if src['ts'] < ts_thresh:
            continue

        gta.set_source_spectrum(src.name, spectrum_type='LogParabola',
                                spectrum_pars={'beta' : {'value' : 0.0, 'scale' : 1.0,
                                                         'min' : 0.0, 'max' : 1.0}})
        optimize_source(gta, src.name)

def main():
        
    usage = "usage: %(prog)s [config file]"
    description = "Run fermipy analysis chain."
    parser = argparse.ArgumentParser(usage=usage,description=description)

    parser.add_argument('--config', default = 'sample_config.yaml')
    parser.add_argument('--source', default = None)
    parser.add_argument('--radius', default = 1.0, type=float)

    args = parser.parse_args()
    gta = GTAnalysis(args.config,logging={'verbosity' : 3},
                     fileio={'workdir_regex' : '\.xml$|\.npy$'})

    if args.source is None:
        src_name = gta.config['selection']['target']
    else:
        src_name = args.source


    max_iter = 5
    if src_name in ['3FGL J0425.8+5600','3FGL J2125.8+5832','FHES J1741.4-3857e']:
        args.radius = 2.0        
    elif 'FHES' in src_name or src_name in ['3FGL J0007.0+7302','3FGL J1725.0-0513']:
        args.radius = 1.5
        
    cat = Catalog3FGL('/u/gl/mdwood/fermi/catalogs/gll_psc_v16_ext.fit')

    flag_mask = 0
    flags = [5,6,8]
    for t in flags:
        flag_mask += 2**t
    
    # Delete unassociated and low TS 3FGL sources
    for s in gta.roi.sources:
        if s.name == src_name:
            continue
        if s.diffuse or s.extended:
            continue
        if s['class']:
            continue
        if not '3FGL' in s.name:
            continue
        if s['offset'] < 0.01:
            continue
        
        m = cat.table['Source_Name'] == s.name
        
        if np.sum(m) and cat.table[m]['TS'] > 100 and (cat.table[m]['Flags']&flag_mask)==0:
            continue

        if np.sum(m):
            gta.logger.info('Deleting %s TS = %.2f Flags = %i',
                            s.name,cat.table[m]['TS'],cat.table[m]['Flags'])
        
        gta.delete_source(s.name)

    if '3FGL J0534.5+2201s' in gta.roi:
        gta.delete_source('3FGL J0534.5+2201s')
        
    gta.setup(overwrite=True)
    
    if '3FGL J0534.5+2201' in gta.roi:
        gta.set_parameter('3FGL J0534.5+2201','Prefactor',0.6057114175,scale=1E-9,true_value=False)
        gta.set_parameter('3FGL J0534.5+2201','Index1',2.237024994,scale=-1.0,true_value=False)
        gta.set_parameter('3FGL J0534.5+2201','Cutoff',1.540488107,scale=10000.,true_value=False)
        gta.set_parameter('3FGL J0534.5+2201','Scale',635.5911255,scale=1.0,true_value=False)
        gta.lock_source('3FGL J0534.5+2201')

    #if '3FGL J0534.5+2201s' in gta.roi:
    #    gta.lock_source('3FGL J0534.5+2201s')
        
    gta.print_roi()
    
    sqrt_ts_threshold=3
    
    model0 = { 'SpatialModel' : 'PointSource', 'Index' : 1.5 }
    model1 = { 'SpatialModel' : 'PointSource', 'Index' : 2.0 }
    model2 = { 'SpatialModel' : 'PointSource', 'Index' : 2.7 }
    model3 = { 'SpatialModel' : 'RadialDisk', 'Index' : 2.0, 'SpatialWidth' : 0.1 }

    newsrc_model = { 'SpectrumType' : 'PowerLaw', 'Index' : 2.0 }
#                     'beta' : {'value' : 0.0, 'min' : 0.0, 'max' : 1.0} }

    skip_loc = ['3FGL J0534.5+2201s']
    
    for s in gta.roi.sources:
        if s['SpectrumType'] == 'LogParabola':
            gta.set_parameter_bounds(s.name, 'beta', [0.0,1.0])
        
    # -----------------------------------
    # Fit the Baseline Model
    # -----------------------------------
    gta.free_sources(skydir=gta.roi[src_name].skydir,distance=1.5, pars='norm')
    gta.fit()
    gta.free_sources(False, pars='norm')
    
    gta.optimize()
    gta.print_roi()
    update_to_lp(gta,100.,exclude=[src_name])

    # Update source of interest to LP if curvature is significant or TS > 1000
    curv = gta.curvature(src_name)
    if curv.ts_curv > 9.0 or gta.roi[src_name]['ts'] > 1000.:
        update_to_lp(gta,100.,names=[src_name])
    
    gta.optimize()

    gta.free_sources(False)
    gta.free_sources(skydir=gta.roi[src_name].skydir,distance=1.5, pars='norm')
    gta.free_source(src_name)
    gta.fit()
    gta.free_sources(False)
    
    gta.print_roi()
    
    # Localize all point sources
    localize(gta, ['PointSource'], 'base0', skip_loc, src_name)
    
    gta.write_roi('base0_roi')            
    gta.tsmap('base0',model=model1, make_plots=True)
    gta.tsmap('base0_nosource',model=model1, exclude=[src_name], make_plots=True)
    gta.tsmap('base0_nosource',model=model2, exclude=[src_name], make_plots=True)

    # Look for new point sources outside the inner 1.0 deg
    o = gta.find_sources('base1_pass0',model=newsrc_model,
                     search_skydir=gta.roi.skydir,
                     max_iter=5,min_separation=0.5,
                     sqrt_ts_threshold=sqrt_ts_threshold,
                     search_minmax_radius=[args.radius,None])
    update_to_lp(gta,100., names=[src.name for src in o['sources']])
    gta.optimize()

    o = gta.find_sources('base1_pass1',model=newsrc_model,
                     search_skydir=gta.roi.skydir,
                     max_iter=5,min_separation=0.5,
                     sqrt_ts_threshold=sqrt_ts_threshold,
                     search_minmax_radius=[args.radius,None])
    update_to_lp(gta,100., names=[src.name for src in o['sources']])

    localize(gta, ['RadialGaussian','RadialDisk'], 'base1', skip_loc, src_name)
    
    gta.print_roi()

    gta.tsmap('base1',model=model1, make_plots=True)
    gta.tsmap('base1',model=model2, make_plots=True)

    gta.write_roi('base1_roi')

    # -------------------------------------
    # Pass 1 - Source at Localized Position
    # -------------------------------------
    if (not src_name in skip_loc and
        gta.roi[src_name]['SpatialModel'] in ['PointSource','RadialGaussian','RadialDisk']):
        gta.localize(src_name,nstep=7,
                     dtheta_max=max(0.4,gta.roi[src_name]['SpatialWidth']),
                     update=True,prefix='fit0',
                     free_radius=max(0.5,gta.roi[src_name]['SpatialWidth']),
                     make_plots=True)
    
    fit_region(gta,'fit0',src_name)
    #gta.load_roi('fit0_roi')

    # -------------------------------------
    # Pass 2 - 2+ Point Sources
    # -------------------------------------

    srcs = []

    # Fit up to 4 additional sources
    for i in range(1,max_iter):

        srcs_fit = gta.find_sources('fit%i'%i,
                                    model=newsrc_model,
                                    search_skydir=gta.roi.skydir,
                                    max_iter=1,
                                    sources_per_iter=1,
                                    sqrt_ts_threshold=3,
                                    min_separation=0.5,
                                    search_minmax_radius=[None,args.radius])
        update_to_lp(gta,100., names=[src.name for src in srcs_fit['sources']])
        
        if len(srcs_fit['sources']) == 0:
            break

        srcs += srcs_fit['sources']
        if (not src_name in skip_loc and
            gta.roi[src_name]['SpatialModel'] in ['PointSource','RadialGaussian','RadialDisk']):
            gta.localize(src_name,nstep=7,
                         dtheta_max=max(0.4,gta.roi[src_name]['SpatialWidth']),
                         update=True,prefix='fit%i'%i,
                         free_radius=max(0.5,gta.roi[src_name]['SpatialWidth']),
                         make_plots=True)

        # Relocalize new sources
        for src in sorted(srcs, key=lambda t: t['ts'],reverse=True):

            gta.localize(src.name,nstep=7,dtheta_max=0.4,
                         update=True,
                         fix_shape=False if src['ts'] >= 25.0 else True,
                         prefix='fit%i'%i,
                         free_radius=0.5, make_plots=True)

        fit_region(gta,'fit%i'%i,src_name)
        tab = gta.roi.create_table([s.name for s in srcs])
        tab.write(os.path.join(gta.workdir,'fit%i_new_source_data.fits'%i),
                  overwrite=True)

        #gta.load_roi('fit%i_roi'%i)
        
    new_source_data = []
    for s in srcs:
        src_data = gta.roi[s.name].data
        new_source_data.append(copy.deepcopy(src_data))

    np.save(os.path.join(gta.workdir,'new_source_data.npy'),
            new_source_data)
    tab = gta.roi.create_table([s.name for s in srcs])
    tab.write(os.path.join(gta.workdir,'new_source_data.fits'),overwrite=True)
    
    

if __name__ == '__main__':

    main()
