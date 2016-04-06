import os
import sys
import copy
from fermipy.gtanalysis import GTAnalysis
from haloanalysis.fit_funcs import fit_halo

import numpy as np
import itertools
import argparse
import yaml

if __name__ == '__main__':

    usage = "usage: %(prog)s [config file]"
    description = "Run fermipy analysis chain."
    parser = argparse.ArgumentParser(usage=usage,description=description)

    parser.add_argument('--config', default = 'sample_config.yaml')
    parser.add_argument('--source', default = None)

    args = parser.parse_args()

    config = yaml.load(open(args.config,'r'))

    halo_width = np.logspace(-1.5,0,4)
    halo_index = np.array([1.5,2.0,2.5])

    
#    halo_source_dict = {
#        'SpectrumType' : 'PowerLaw', 
#        'Index' : 2.0, 
#        'Scale' : 1000,
#        'Prefactor' : { 'value' : 1E-5, 'scale' : 1e-13 },
#        'SpatialModel' : 'GaussianSource',
#        'SpatialWidth' : 1.0
#        }

    gta = GTAnalysis(args.config,logging={'verbosity' : 3})

    gta.setup()
    gta.simulate_roi(restore=True)

    ext_fit_data = []
    halo_fit_data = []

    gta.write_roi('base_model',save_model_map=False,make_plots=False)

    for i in range(10):

        gta.load_xml('base_model')

        gta.simulate_roi()

        gta.free_source('testsource')
        gta.free_source('galdiff',pars='norm')
        gta.free_source('isodiff',pars='norm')

        gta.fit()
        gta.update_source('testsource',reoptimize=True,npts=9)
        gta.free_sources(free=False)

        gta.extension('testsource',width=np.logspace(-2.5,-0.5,9))

        ext_fit_data += [copy.deepcopy(gta.roi['testsource'])]

        gta.write_roi('fit%04i'%i,save_model_map=False,make_plots=False,
                      format='npy')

        fit_halo(gta,'fit%04i'%i,'testsource',halo_width,halo_index)
        
    np.save(os.path.join(gta._savedir,'ext_fit_data.npy'),ext_fit_data)

