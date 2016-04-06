import os
import sys
import copy
import itertools

import numpy as np
    
def fit_region(gta,modelname,src_name,erange=None):

    if erange is not None:
        gta.setEnergyRange(erange[0],erange[1])
    
    model0 = { 'SpatialModel' : 'PointSource', 'Index' : 1.5 }
    model1 = { 'SpatialModel' : 'PointSource', 'Index' : 2.0 }
    model2 = { 'SpatialModel' : 'PointSource', 'Index' : 2.5 }

    model3 = { 'SpatialModel' : 'Gaussian', 'Index' : 2.0,
               'SpatialWidth' : 0.2 }
    
    gta.optimize()

    gta.free_sources(False)
    gta.free_sources(distance=1.0,exclude_diffuse=True)
    gta.free_sources(distance=1.0,pars='norm')
    gta.fit()

    gta.free_sources(False)
    gta.free_sources(distance=1.0,pars='norm',exclude_diffuse=True)
    
    gta.extension(src_name)
    gta.sed(src_name)
    
    gta.write_roi(modelname)
    gta.tsmap(modelname,model=model0,erange=erange)
    maps_model1 = gta.tsmap(modelname,model=model1,erange=erange)
    gta.tsmap(modelname,model=model2,erange=erange)
    maps_model1_nosource = gta.tsmap('%s_nosource'%modelname,
              model=model1,exclude=[src_name],erange=erange)
    gta.residmap(modelname,model=model3,erange=erange)

    # Make zoom plots
    gta._plotter.make_tsmap_plots(gta,maps_model1,
                                  zoom=2,suffix='tsmap_zoom')
    gta._plotter.make_tsmap_plots(gta,maps_model1_nosource,
                                  zoom=2,suffix='tsmap_zoom')    

    lnl = -gta.like()

    gta.logger.info('%s Model Likelihood: %f'%(modelname,lnl))

def fit_halo(gta,modelname,src_name,halo_width,halo_index,erange=None):

    halo_source_name = 'halo_gauss'
    halo_source_dict = {
        'SpectrumType' : 'PowerLaw', 
        'Index' : 2.0, 
        'Scale' : 1000,
        'Prefactor' : { 'value' : 1E-5, 'scale' : 1e-13 },
        'SpatialModel' : 'GaussianSource',
        'SpatialWidth' : 1.0
        }

    halo_source_dict['ra'] = gta.roi[src_name]['ra']
    halo_source_dict['dec'] = gta.roi[src_name]['dec']
    
    #halo_width = np.logspace(-1,0,9)
    #halo_index = np.array([1.5,1.75,2.0,2.25,2.5,2.75,3.0])
    halo_data = []
    halo_data_idx_free = []

    gta.load_roi(modelname)
    if erange is not None:
        gta.setEnergyRange(erange[0],erange[1])

    gta.free_sources(False)
    gta.free_sources(distance=1.0,pars='norm', exclude_diffuse=True)
    gta.write_xml(modelname + '_base')

    # Find best-fit halo model
    halo_source_dict['SpatialWidth'] = 0.1
    gta.add_source(halo_source_name,halo_source_dict,free=True)
    gta.free_parameter(halo_source_name,'Index')
    gta.extension(halo_source_name,update=True)
    gta.free_parameter(halo_source_name,'Index')
    gta.fit()
    gta.free_parameter(halo_source_name,'Index',False)
    gta.update_source(halo_source_name,reoptimize=True,npts=9)

    gta.write_roi('%s_halo_gauss'%(modelname),make_plots=False,
                  save_model_map=False,format='npy')    
    np.save(os.path.join(gta.workdir,'%s_halo_data.npy'%modelname),
            copy.deepcopy(gta.roi['halo_gauss'].data))
    gta.delete_source(halo_source_name,save_template=False) 
    
    #for i, (w,idx) in enumerate(itertools.product(halo_width,halo_index)):
    for i, w in enumerate(halo_width):
        halo_source_dict['SpatialWidth'] = w
        gta.load_xml(modelname + '_base')
        
        gta.add_source(halo_source_name,halo_source_dict,free=True)

        # Do one fit with index free
        gta.set_parameter(halo_source_name,'Index',-2.0,
                          update_source=False)

        gta.fit()
        
        # SED
        gta.sed(halo_source_name)    
        gta.write_roi('%s_halo_gauss_sed_%02i'%(modelname,i),make_plots=False,
                      save_model_map=False,format='npy')

        # Free Index
        gta.free_parameter(halo_source_name,'Index')
        gta.fit()
        gta.free_parameter(halo_source_name,'Index',False)
        gta.update_source(halo_source_name,reoptimize=True,npts=9)

        halo_data_idx_free += [copy.deepcopy(gta.roi['halo_gauss'].data)]
        gta.write_roi('%s_halo_gauss_%02i'%(modelname,i),make_plots=False,
                      save_model_map=False,format='npy')
        
        # Scan over fixed index
        for j, idx in enumerate(halo_index):

            model_idx = i*len(halo_index) + j

            gta.set_parameter(halo_source_name,'Index',-1.0*idx,
                              update_source=False)
            gta.fit(update=False)
    
            gta.update_source(halo_source_name,reoptimize=True, npts=9)

            gta.logger.info('%s Halo Width: %6.3f Index: %6.2f TS: %6.2f'%(modelname,w,idx,
                                                                           gta.roi[halo_source_name]['ts']))
    
            gta.write_roi('%s_halo_gauss_%02i_%02i'%(modelname,i,j),make_plots=False,
                          save_model_map=False,format='npy')
            halo_data += [copy.deepcopy(gta.roi['halo_gauss'].data)]
            
        gta.delete_source(halo_source_name,save_template=False) 

    np.save(os.path.join(gta.workdir,'%s_halo_data.npy'%modelname),halo_data)
    np.save(os.path.join(gta.workdir,'%s_halo_data_idx_free.npy'%modelname),
            halo_data_idx_free)
