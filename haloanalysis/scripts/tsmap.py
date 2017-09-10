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
from astropy.table import Table
    
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
    tab = Table.read('/u/gl/mdwood/fermi/ext_analysis/v20/std_psf0123_joint2a_stdmodel/table_std_psf0123_joint2a_stdmodel_cat.fits')
    src_name = gta.config['selection']['target']
    #codename = src_name.lower().replace(' ','_').replace('_off','')
    row = tab[tab['name_roi'] == src_name][0]

    if row['fit_ext_ts_ext'] < 9.0:
        print('Source not extended.  Exiting.')
        sys.exit(0)

    ext_model = row['fit_ext_model']
    idx = row['fit_idx_ext']
    model = { 'SpatialModel' : 'PointSource', 'Index' : min(np.abs(row['fit_ext_index']),4.0) }
    
    gta.load_roi('fit%i_ext_%s_roi'%(idx,ext_model),reload_sources=True)
    
    gta.tsmap('fit_ext_nosource', outfile='fit_ext_nosource_tsmap',
              model=model, exclude=[src_name],
              make_plots=True)

    gta.tsmap('fit_ext', outfile='fit_ext_tsmap',
              model=model, 
              make_plots=True)

if __name__ == '__main__':

    main()
