import os
import sys
import copy
from fermipy.gtanalysis import GTAnalysis
import numpy as np
import itertools
import argparse
from haloanalysis.fit_funcs import fit_region, fit_halo
from haloanalysis.batch import check_log

if __name__ == '__main__':
        
    usage = "usage: %(prog)s [config file]"
    description = "Run fermipy analysis chain."
    parser = argparse.ArgumentParser(usage=usage,description=description)

    parser.add_argument('--config', default = 'sample_config.yaml')
    parser.add_argument('--source', default = None)

    args = parser.parse_args()
    gta = GTAnalysis(args.config,logging={'verbosity' : 3})

    logfile = os.path.join(gta.outdir,'run_analysis.log')

    if not check_log(logfile)=='Successful':
        sys.exit(1)
    
    gta.setup()

    sqrt_ts_threshold=3

    halo_width = np.logspace(-1,0,9)
    halo_index = np.array([1.5,1.75,2.0,2.25,2.5,2.75,3.0])
    
    model0 = { 'SpatialModel' : 'PointSource', 'Index' : 1.5 }
    model1 = { 'SpatialModel' : 'PointSource', 'Index' : 2.0 }
    model2 = { 'SpatialModel' : 'PointSource', 'Index' : 2.5 }
    src_name = gta.roi.sources[0].name
    
    gta.load_roi('base',reload_sources=True)
    #gta.tsmap('base',model=model1)
    gta.tsmap('base_emin40',model=model1,erange=[4.0,5.5])

    gta.print_roi()

    # -----------------------------------
    # Pass 0 - Source at Nominal Position
    # -----------------------------------

    gta.load_roi('fit0',reload_sources=True)
    
    #fit_region(gta,'fit0',src_name)
    fit_region(gta,'fit0_emin40',src_name,erange=[4.0,5.5])

    # -------------------------------------
    # Pass 1 - Source at Localized Position
    # -------------------------------------

    gta.load_roi('fit1',reload_sources=True)
    
    #fit_region(gta,'fit1',src_name)
    #fit_halo(gta,'fit1',src_name,halo_width,halo_index)

    fit_region(gta,'fit1_emin40',src_name,erange=[4.0,5.5])
    fit_halo(gta,'fit1_emin40',src_name,halo_width,halo_index,erange=[4.0,5.5])

    # -------------------------------------
    # Pass 2 - 2+ Point Sources
    # -------------------------------------

    best_fit_idx = 1

    # Fit up to 4 sources
    for i in range(2,5):

        roi_file = 'fit%i.npy'%i

        if not os.path.isfile(os.path.join(gta.workdir,roi_file)):
            continue
        
        best_fit_idx = i
        
#        fit_region(gta,'fit%i'%i,src_name)
#        fit_halo(gta,'fit%i'%i,src_name,halo_width,halo_index,
#                 do_scan=False)
        gta.load_roi('fit%i'%i,reload_sources=True)
        fit_region(gta,'fit%i_emin40'%i,src_name,erange=[4.0,5.5])
        fit_halo(gta,'fit%i_emin40'%i,src_name,halo_width,
                 halo_index,erange=[4.0,5.5],
                 do_scan=False)



    # Only Run Halo Fit for Best-fit Model
    if best_fit_idx > 1:
#        fit_halo(gta,'fit%i'%best_fit_idx,src_name,
#                 halo_width,halo_index)
        fit_halo(gta,'fit%i_emin40'%best_fit_idx,src_name,
                 halo_width,halo_index,
                 erange=[4.0,5.5])
        
