import os
import sys
import copy
import numpy as np
import itertools
import argparse

from fermipy.utils import init_matplotlib_backend

init_matplotlib_backend()

from fermipy.batch import check_log
from fermipy.gtanalysis import GTAnalysis
from haloanalysis.fit_funcs import fit_halo_scan
    
def main():
        
    usage = "usage: %(prog)s [config file]"
    description = "Run fermipy analysis chain."
    parser = argparse.ArgumentParser(usage=usage,description=description)

    parser.add_argument('--config', default = 'sample_config.yaml')
    parser.add_argument('--source', default = None)

    args = parser.parse_args()
    gta = GTAnalysis(args.config,logging={'verbosity' : 3},
                     fileio={'workdir_regex' : ['\.xml$|\.npy$']})

    logfile = os.path.join(gta.outdir,'run-region-analysis.log')
    if not check_log(logfile)=='Successful':
        print('Region analysis incomplete.  Exiting.')
        sys.exit(1)
        
    gta.setup()

    halo_width = np.logspace(-1.25,0.25,13)
    halo_index = np.array([1.5,1.75,2.0,2.25,2.5,2.75,3.0,3.25,3.5])
    src_name = gta.roi.sources[0].name

    for i in range(1,6):

        npy_file = os.path.join(gta.workdir,'fit%i.npy'%i)

        if not os.path.isfile(npy_file):
            continue

        gta.load_roi('fit%i'%i,reload_sources=True)
        fit_halo_scan(gta,'fit%i'%i,src_name,
                      halo_width,halo_index,optimizer='NEWTON')

if __name__ == '__main__':

    main()
