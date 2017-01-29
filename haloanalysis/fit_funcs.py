import os
import sys
import copy
import itertools
import logging

import numpy as np
    
def fit_region(gta,modelname,src_name,loge_bounds=None):

    gta.logger.info('Starting Region Fit %s'%(modelname))
    lnl0 = -gta.like()    
    gta.logger.info('%s Model Likelihood: %f'%(modelname,lnl0))
    gta.print_params()
    
    if loge_bounds is not None:
        gta.set_energy_range(loge_bounds[0],loge_bounds[1])
    
    model0 = { 'SpatialModel' : 'PointSource', 'Index' : 1.5 }
    model1 = { 'SpatialModel' : 'PointSource', 'Index' : 2.0 }
    model2 = { 'SpatialModel' : 'PointSource', 'Index' : 2.7 }

    model3 = { 'SpatialModel' : 'Gaussian', 'Index' : 2.0,
               'SpatialWidth' : 0.15 }
    
    gta.optimize()

    diff_sources = [s.name for s in gta.roi.sources if s.diffuse]
    skydir = gta.roi[src_name].skydir
    
    gta.free_sources(False)
    #gta.free_sources(skydir=skydir,distance=1.0, exclude=diff_sources)
    gta.free_sources(skydir=skydir,distance=1.0, pars='norm')
    gta.fit()

    gta.free_sources(skydir=skydir,distance=1.0, pars='norm',
                     exclude=diff_sources)

    gta.extension(src_name, outfile=modelname + '_ext_gauss_fixed',
                  spatial_model='RadialGaussian',
                  prefix=modelname + '_gauss_fixed',
                  optimizer={'optimizer' : 'NEWTON'},
                  fit_position=False, free_radius=1.0,
                  make_plots=True)
    
    gta.extension(src_name, outfile=modelname + '_ext_gauss',
                  spatial_model='RadialGaussian',
                  prefix=modelname + '_gauss',
                  optimizer={'optimizer' : 'NEWTON'},
                  fit_position=True, free_radius=1.0,
                  make_plots=True)

    gta.extension(src_name, outfile=modelname + '_ext_disk_fixed',
                  spatial_model='RadialDisk',
                  prefix=modelname + '_disk_fixed',
                  optimizer={'optimizer' : 'NEWTON'},
                  fit_position=False, free_radius=1.0,
                  make_plots=True)
    
    gta.extension(src_name, outfile=modelname + '_ext_disk',
                  spatial_model='RadialDisk',
                  prefix=modelname + '_disk',
                  optimizer={'optimizer' : 'NEWTON'},
                  fit_position=True, free_radius=1.0,
                  make_plots=True)

    gta.sed(src_name, outfile=modelname + '_sed_fixed',
            prefix=modelname + '_fixed',
            optimizer={'optimizer' : 'MINUIT'},
            make_plots=True)
    
    gta.sed(src_name, outfile=modelname + '_sed',
            prefix=modelname,
            optimizer={'optimizer' : 'MINUIT'},
            free_radius=1.0, make_plots=True)

    gta.write_roi(modelname, make_plots=True)
    gta.tsmap(modelname, model=model0,
              loge_bounds=loge_bounds, make_plots=True)
    maps_model1 = gta.tsmap(modelname, model=model1,
                            loge_bounds=loge_bounds, make_plots=True)
    gta.tsmap(modelname, model=model2,
              loge_bounds=loge_bounds, make_plots=True)
    maps_model1_nosource = gta.tsmap('%s_nosource'%modelname,
                                     model=model1, exclude=[src_name],
                                     loge_bounds=loge_bounds, make_plots=True)
    gta.residmap(modelname, model=model3,
                 loge_bounds=loge_bounds, make_plots=True)

    # Make zoom plots
    gta.plotter.make_tsmap_plots(maps_model1, gta.roi,
                                  zoom=2,suffix='tsmap_zoom')
    gta.plotter.make_tsmap_plots(maps_model1_nosource, gta.roi,
                                  zoom=2,suffix='tsmap_zoom')    

    lnl1 = -gta.like()

    gta.print_roi()
    gta.print_params()
    
    gta.logger.info('%s Model Likelihood: %f'%(modelname,lnl1))
    gta.logger.info('%s Model Likelihood Delta: %f'%(modelname,lnl1-lnl0))
    gta.logger.info('Finished Region Fit %s'%(modelname))


def fit_halo_sed(gta,modelname,src_name,halo_width,
                 halo_index,spatial_model='RadialGaussian',
                 loge_bounds=None):

    gta.logger.info('Starting Halo SED Fit %s'%(modelname))
    
    halo_source_name = 'halo_' + spatial_model
    halo_source_dict = {
        'SpectrumType' : 'PowerLaw', 
        'Index' : { 'value' : 2.0, 'scale' : -1.0, 'min' : 1.0, 'max' : 4.5 }, 
        'Scale' : 1000,
        'Prefactor' : { 'value' : 1E-5, 'scale' : 1e-13, 'min' : 1E-5, 'max' : 1E4 },
        'SpatialModel' : spatial_model,
        'SpatialWidth' : 1.0
        }

    halo_source_dict['ra'] = gta.roi[src_name]['ra']
    halo_source_dict['dec'] = gta.roi[src_name]['dec']

    gta.load_roi(modelname)
    if loge_bounds is not None:
        gta.set_energy_range(loge_bounds[0],loge_bounds[1])


    diff_sources = [s.name for s in gta.roi.sources if s.diffuse]
        
    gta.free_sources(False)
    gta.free_sources(distance=1.0,pars='norm', exclude=diff_sources)
    gta.write_xml(modelname + '_base')
    
    for i, w in enumerate(halo_width):

        halo_source_dict['SpatialWidth'] = w        
        gta.load_xml(modelname + '_base')
        
        gta.add_source(halo_source_name,halo_source_dict,free=True)

        # Do one fit with index free
        gta.set_parameter(halo_source_name,'Index',-2.0,
                          update_source=False)

        gta.fit()
        
        # SED w/ Index = 2.0
        gta.sed(halo_source_name,prefix='%s_%02i'%(modelname,i),
                fix_background=False, cov_scale=5.0)
        gta.write_roi('%s_halo_gauss_sed_%02i'%(modelname,i),
                      make_plots=False)

    gta.logger.info('Finished Halo SED Fit %s'%(modelname))

def fit_halo_scan(gta, modelname, src_name, halo_width,
                  halo_index, spatial_model='RadialGaussian',
                  loge_bounds=None, optimizer='NEWTON'):

    gta.logger.info('Starting Halo Scan %s'%(modelname))

    halo_source_name = 'halo_' + spatial_model
    halo_source_dict = {
        'SpectrumType' : 'PowerLaw', 
        'Index' : { 'value' : 2.0, 'scale' : -1.0, 'min' : 1.0, 'max' : 4.5 }, 
        'Scale' : 1000,
        'Prefactor' : { 'value' : 1E-5, 'scale' : 1e-13,
                        'min' : 1E-5, 'max' : 1E4 },
        'SpatialModel' : spatial_model,
        'SpatialWidth' : 1.0
        }

    outprefix = '%s_%s'%(modelname,halo_source_name)
    
    halo_source_dict['ra'] = gta.roi[src_name]['ra']
    halo_source_dict['dec'] = gta.roi[src_name]['dec']

    #gta.load_roi(modelname)
    #if loge_bounds is not None:
    #    gta.set_energy_range(loge_bounds[0],loge_bounds[1])

    skydir = gta.roi[src_name].skydir
    diff_sources = [s.name for s in gta.roi.sources if s.diffuse]
    
    gta.free_sources(False)
    gta.free_sources(skydir=skydir,distance=1.0,pars='norm',
                     exclude=diff_sources)
    gta.write_xml(modelname + '_base')

    halo_tab = gta.roi.create_table([])
    halo_data = []
    halo_data_idx_free = []
        
    for i, w in enumerate(halo_width):

        gta.logger.info('Fitting Halo Width %.3f',w)
        
        halo_source_dict['SpatialWidth'] = w
        gta.load_xml(modelname + '_base')
        
        gta.add_source(halo_source_name, halo_source_dict, free=True)

        # Free Index
        gta.free_norm(halo_source_name)
        gta.fit(optimizer=optimizer)
        gta.sed(halo_source_name, prefix='%s_cov05_%02i'%(modelname,i),
                free_radius=1.0, cov_scale=5.0,
                optimizer={'optimizer' : 'MINUIT'},
                make_plots=True)
        
        gta.free_parameter(halo_source_name,'Index')
        gta.fit(optimizer=optimizer)
        gta.free_parameter(halo_source_name,'Index',False)
        gta.update_source(halo_source_name,reoptimize=True,
                          optimizer={'optimizer' : optimizer})

        halo_data_idx_free += [copy.deepcopy(gta.roi[halo_source_name].data)]
        gta.write_roi('%s_%02i'%(outprefix,i),make_plots=False)

        gta.print_params(loglevel=logging.DEBUG)
        
        # Scan over fixed index
        for j, idx in enumerate(halo_index):

            gta.logger.info('Fitting Halo Index %.3f',idx)
            
            model_idx = i*len(halo_index) + j
            gta.set_norm(halo_source_name, 1.0, update_source=False)            
            gta.set_parameter(halo_source_name, 'Index', -1.0*idx,
                              update_source=False)
            
            gta.fit(update=False, optimizer=optimizer)

            gta.print_params(loglevel=logging.DEBUG)
            
            gta.update_source(halo_source_name,reoptimize=True,
                              optimizer={'optimizer' : optimizer})

            gta.logger.info('%s Halo Width: %6.3f Index: %6.2f TS: %6.2f',
                            modelname,w,idx,gta.roi[halo_source_name]['ts'])
    
            #gta.write_roi('%s_%02i_%02i'%(outprefix,i,j),make_plots=False)
            halo_data += [copy.deepcopy(gta.roi[halo_source_name].data)]
            gta.roi[halo_source_name].add_to_table(halo_tab)
            
        gta.delete_source(halo_source_name,save_template=False) 

    np.save(os.path.join(gta.workdir,'%s_data.npy'%outprefix),halo_data)
    np.save(os.path.join(gta.workdir,'%s_data_idx_free.npy'%outprefix),
            halo_data_idx_free)

    tab_halo_width, tab_halo_index = np.meshgrid(halo_width,halo_index,indexing='ij')
    halo_tab['halo_width'] = np.ravel(tab_halo_width)
    halo_tab['halo_index'] = np.ravel(tab_halo_index)
    
    halo_tab.write(os.path.join(gta.workdir,'%s_data.fits'%outprefix),overwrite=True)
    gta.logger.info('Finished Halo Scan %s'%(modelname))
    
    
def fit_halo(gta, modelname, src_name,
             spatial_model='RadialGaussian',
             loge_bounds=None, optimizer='NEWTON'):

    gta.logger.info('Starting Halo Fit %s'%(modelname))
    
    halo_source_name = 'halo_' + spatial_model
    halo_source_dict = {
        'SpectrumType' : 'PowerLaw', 
        'Index' : { 'value' : 2.0, 'scale' : -1.0, 'min' : 1.0, 'max' : 4.5 }, 
        'Scale' : 1000,
        'Prefactor' : { 'value' : 1E-5, 'scale' : 1e-13,
                        'min' : 1E-5, 'max' : 1E4 },
        'SpatialModel' : spatial_model,
        'SpatialWidth' : 1.0
        }

    outprefix = '%s_%s'%(modelname,halo_source_name)
    
    halo_source_dict['ra'] = gta.roi[src_name]['ra']
    halo_source_dict['dec'] = gta.roi[src_name]['dec']

#    gta.load_roi(modelname)
#    if loge_bounds is not None:
#        gta.set_energy_range(loge_bounds[0],loge_bounds[1])

    diff_sources = [s.name for s in gta.roi.sources if s.diffuse]
    
    gta.free_sources(False)
    gta.free_sources(distance=1.0,pars='norm',
                     exclude=diff_sources)

    # Find best-fit halo model
    halo_source_dict['SpatialWidth'] = 0.1
    gta.add_source(halo_source_name,halo_source_dict)
    gta.free_norm(halo_source_name)
    gta.extension(halo_source_name,update=True,
                  optimizer={'optimizer' : optimizer},
                  free_radius=1.0)

    # Fit spectrum
    gta.free_parameter(halo_source_name,'Index')
    gta.fit()

    # Re-fit extension
    gta.extension(halo_source_name,update=True,
                  optimizer={'optimizer' : optimizer},
                  free_radius=1.0)    

    # Re-fit Spectrum
    gta.fit()

    gta.update_source(halo_source_name,reoptimize=True,
                      optimizer={'optimizer' : optimizer})

    gta.print_params()
    
    gta.write_roi(outprefix,make_plots=False)    
    np.save(os.path.join(gta.workdir,'%s_data.npy'%outprefix),
            copy.deepcopy(gta.roi[halo_source_name].data))
    gta.delete_source(halo_source_name,save_template=False)
    
    gta.logger.info('Finished Halo Fit %s'%(modelname))
