import numpy as np
import glob
import sys
import os
import argparse

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

cols = [Column([], name='name', dtype='S20', format='%s',description='Source Name'),
        Column([], name='codename', dtype='S20', format='%s'),
        Column([], name='linkname', dtype='S20', format='%s'),
        Column([], name='assoc', dtype='S20', format='%s'),
        Column([], name='class', dtype='S20', format='%s'),
        Column([], name='ra', dtype='f8', format='%.3f',unit='deg'),
        Column([], name='dec', dtype='f8', format='%.3f',unit='deg'),
        Column([], name='glon', dtype='f8', format='%.3f',unit='deg'),
        Column([], name='glat', dtype='f8', format='%.3f',unit='deg'),
        Column([], name='ts', dtype='f8', format='%.2f'),
        Column([], name='npred', dtype='f8', format='%.2f'),
        Column([], name='ext0_ts', dtype='f8', format='%.2f'),
        Column([], name='ext0_mle', dtype='f8', format='%.3f',unit='deg'),
        Column([], name='ext0_err', dtype='f8', format='%.3f',unit='deg'),
        Column([], name='ext0_ul95', dtype='f8', format='%.3f',unit='deg',description='Extension 95% CL UL (pre-localization)'),
        Column([], name='ext1_ts', dtype='f8', format='%.2f'),
        Column([], name='ext1_mle', dtype='f8', format='%.3f',unit='deg'),
        Column([], name='ext1_err', dtype='f8', format='%.3f',unit='deg'),
        Column([], name='ext1_ul95', dtype='f8', format='%.3f',unit='deg',description='Extension 95% CL UL (post-localization)'),
        ]


tab = Table(cols)

for d in dirs:

    print d
    
    file0 = os.path.join(d,'fit0.npy')
    file1 = os.path.join(d,'fit1.npy')
    file2 = os.path.join(d,'halo_data.npy')

    if not os.path.isfile(file2):
        continue
    
    if not os.path.isfile(file0) or not os.path.isfile(file1):
        continue
    
    data0 = np.load(file0).flat[0]
    data1 = np.load(file1).flat[0]
    data2 = np.load(file2)
    
    for k,v in data0['sources'].items():
        if v['extension'] is not None:
            name = v['name']
            break
        
    src0 = data0['sources'][name]
    src1 = data1['sources'][name]

    ext0 = data0['sources'][name]['extension']
    ext1 = data1['sources'][name]['extension']

    row = [src0['assoc']['ASSOC1'],src0['class'],
           src0['ra'],src0['dec'],src0['glon'],src0['glat'],
           src0['ts'],src0['Npred'],
           max(ext0['ts_ext'],0),ext0['ext'],ext0['ext_err'],ext0['ext_ul95'],
           max(ext1['ts_ext'],0),ext1['ext'],ext1['ext_err'],ext1['ext_ul95'],]


    codename = os.path.basename(d)
    linkname = '{%s}'%codename.replace('+','p').replace('.','_')
    
    for hd in data2:
        colname = 'halo_%.2f_%.3f_ts'%(np.abs(hd.data['params']['Index'][0]),
                                       hd.data['SpatialWidth'])
        if not colname in tab.colnames:
            tab.add_column(Column([], name=colname, dtype='f8', format='%.3f'))        
        row += [hd.data['ts']]

        colname = 'halo_%.2f_%.3f_eflux_ul95'%(np.abs(hd.data['params']['Index'][0]),
                                               hd.data['SpatialWidth'])
        if not colname in tab.colnames:
            tab.add_column(Column([], name=colname, dtype='f8', format='%.3f'))
        
        row += [hd.data['eflux_ul95']]
            
    tab.add_row([name,codename,linkname] + row)

    
m = tab['class']==''
tab['class'][m] = 'unkn'
    
tab.write(args.output,format='fits',overwrite=True)
