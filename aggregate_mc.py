
import glob
import sys
import os
import argparse

import numpy as np
from astropy.table import Table, Column
from fermipy.roi_model import ROIModel
from haloanalysis.utils import collect_dirs

usage = "usage: %(prog)s"
description = "Aggregate analysis output."
parser = argparse.ArgumentParser(usage=usage,description=description)


parser.add_argument('--output', default = 'test.fits')
parser.add_argument('dirs', nargs='+', default = None,
                    help='Run analyses in all subdirectories of this '
                    'directory.')

args = parser.parse_args()
dirs = sorted(collect_dirs(args.dirs))

cols = [Column([], name='name', dtype='S20', format='%s'),        
        Column([], name='ra', dtype='f8', format='%.3f',unit='deg'),
        Column([], name='dec', dtype='f8', format='%.3f',unit='deg'),
        Column([], name='glon', dtype='f8', format='%.3f',unit='deg'),
        Column([], name='glat', dtype='f8', format='%.3f',unit='deg'),
        Column([], name='index', dtype='f8', format='%.3f',unit='deg'),
        Column([], name='prefactor', dtype='f8', format='%.4g',unit='deg'),
        Column([], name='eflux', dtype='f8', format='%.4g',unit='deg'),
        Column([], name='ts', dtype='f8', format='%.2f'),
        Column([], name='npred', dtype='f8', format='%.2f'),
        Column([], name='ext_ts', dtype='f8', format='%.2f'),
        Column([], name='ext_mle', dtype='f8', format='%.3f',unit='deg'),
        Column([], name='ext_err', dtype='f8', format='%.3f',unit='deg'),
        Column([], name='ext_ul95', dtype='f8', format='%.3f',unit='deg'),
        ]

tab = Table(cols)

for d in dirs:

    print d
    
    file0 = os.path.join(d,'base_model.npy')
    file1 = os.path.join(d,'ext_fit_data.npy')
    
    if not os.path.isfile(file0) or not os.path.isfile(file1):
        continue

    
    data0 = np.load(file0).flat[0]
    data1 = np.load(file1)
    src0 = data0['sources']['testsource']
    
    for s in data1:
        
        name = '%03.1f_%03.1f'%(-src0['params']['Index'][0],
                                np.log10(src0['params']['Prefactor'][0]/1E-13))

        ext = s['extension']
        row = [name,s['ra'],s['dec'],s['glon'],s['glat'],
               -src0['params']['Index'][0],src0['params']['Prefactor'][0],
               src0['eflux'][0],s['ts'],s['Npred'],
               max(ext['ts_ext'],0),ext['ext'],ext['ext_err'],ext['ext_ul95']]
        tab.add_row(row)
        
tab.write(args.output,format='fits',overwrite=True)
