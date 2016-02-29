import numpy as np
import glob
import sys
import os
import traceback
import argparse

from astropy.table import Table, Column
from fermipy.roi_model import ROIModel
from haloanalysis.utils import collect_dirs

def find_source(data):

    srcs = sorted(data['sources'].values(), key=lambda t: t['offset'])

    for s in srcs:
        if s['SourceType'] == 'DiffuseSource':
            continue

        return s
    
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
        Column([], name='fit0_ts', dtype='f8', format='%.2f'),
        Column([], name='fit1_ts', dtype='f8', format='%.2f'),
        Column([], name='fit2_ts', dtype='f8', format='%.2f'),
        Column([], name='fit0_offset', dtype='f8', format='%.2f'),
        Column([], name='fit1_offset', dtype='f8', format='%.2f'),
        Column([], name='fit2_offset', dtype='f8', format='%.2f'),
        #block added by RC
        Column([], name='dfde1000', dtype='f8', format='%.2f'),
        Column([], name='dfde1000_err', dtype='f8', format='%.2f'),
        Column([], name='dfde1000_index', dtype='f8', format='%.2f'),
        Column([], name='dfde1000_index_err', dtype='f8', format='%.2f'),
        Column([], name='flux1000', dtype='f8', format='%.2f'),
        Column([], name='flux1000_err', dtype='f8', format='%.2f'),
        Column([], name='eflux1000', dtype='f8', format='%.2f'),
        Column([], name='eflux1000_err', dtype='f8', format='%.2f'),
        Column([], name='flux100', dtype='f8', format='%.2f'),
        Column([], name='flux100_err', dtype='f8', format='%.2f'),
        Column([], name='eflux100', dtype='f8', format='%.2f'),
        Column([], name='eflux100_err', dtype='f8', format='%.2f'),
        Column([], name='spectrum_type', dtype='S20', format='%s'),
        #Column([], name='lnlprofile', dtype='f8', format='%.2f'),
        Column([], name='ext0_ts', dtype='f8', format='%.2f'),
        Column([], name='ext0_mle', dtype='f8', format='%.3f',unit='deg'),
        Column([], name='ext0_err', dtype='f8', format='%.3f',unit='deg'),
        Column([], name='ext0_ul95', dtype='f8', format='%.3f',unit='deg',description='Extension 95% CL UL (pre-localization)'),
        #block added by RC
        #Column([], name='ext0_dlogLike', dtype='f8', format='%.2f'),
        #Column([], name='ext0_width', dtype='f8', format='%.3f',unit='deg'),
        Column([], name='ext1_ts', dtype='f8', format='%.2f'),
        Column([], name='ext1_mle', dtype='f8', format='%.3f',unit='deg'),
        Column([], name='ext1_err', dtype='f8', format='%.3f',unit='deg'),
        Column([], name='ext1_ul95', dtype='f8', format='%.3f',unit='deg',description='Extension 95% CL UL (post-localization)'),
        #block added by RC
        #Column([], name='ext1_dlogLike', dtype='f8', format='%.2f'),
        #Column([], name='ext1_width', dtype='f8', format='%.3f',unit='deg'),
        Column([], name='ext2_ts', dtype='f8', format='%.2f'),
        Column([], name='ext2_mle', dtype='f8', format='%.3f',unit='deg'),
        Column([], name='ext2_err', dtype='f8', format='%.3f',unit='deg'),
        Column([], name='ext2_ul95', dtype='f8', format='%.3f',unit='deg',description='Extension 95% CL UL (post-localization)'),
        Column([], name='fit1_dlike', dtype='f8', format='%.2f'),
        Column([], name='fit2_dlike', dtype='f8', format='%.2f'),
        Column([], name='fit2_nsrc', dtype='i8'),
        ]


tab = Table(cols)

for d in dirs:

    print d
    
    file0 = os.path.join(d,'fit0.npy')
    file1 = os.path.join(d,'fit1.npy')
    file2 = os.path.join(d,'fit2.npy')
    file3 = os.path.join(d,'fit1_halo_data.npy')
    file4 = os.path.join(d,'fit2_halo_data.npy')
    file5 = os.path.join(d,'new_source_data.npy')

    if not os.path.isfile(file4):
        continue
    
    if not os.path.isfile(file0) or not os.path.isfile(file1):
        continue

    try:
        data0 = np.load(file0).flat[0]
        data1 = np.load(file1).flat[0]
        data2 = np.load(file2).flat[0]
        halo_data1 = np.load(file3)
        halo_data2 = np.load(file4)
        new_srcs = np.load(file5)
    except Exception as e:
        traceback.print_exc()
        continue
        
    src0 = find_source(data0)
    src1 = find_source(data1)
    src2 = find_source(data2)
        
    ext0 = src0['extension']
    ext1 = src1['extension']
    ext2 = src2['extension']

    ext0_results = [np.nan,np.nan,np.nan,np.nan]
    ext1_results = [np.nan,np.nan,np.nan,np.nan]
    ext2_results = [np.nan,np.nan,np.nan,np.nan]

    if ext0 is not None:
        ext0_results = [max(ext0['ts_ext'],0),ext0['ext'],ext0['ext_err'],ext0['ext_ul95']]

    if ext1 is not None:
        ext1_results = [max(ext1['ts_ext'],0),ext1['ext'],ext1['ext_err'],ext1['ext_ul95']]

    if ext2 is not None:
        ext2_results = [max(ext2['ts_ext'],0),ext2['ext'],ext2['ext_err'],ext2['ext_ul95']]
    
    fit1_dlike = data1['roi']['logLike'] - data0['roi']['logLike']
    fit2_dlike = data2['roi']['logLike'] - data0['roi']['logLike']
    nsrc = len(new_srcs)
    
    row = [src0['assoc']['ASSOC1'],src0['class'],
           src1['ra'],src1['dec'],src1['glon'],src1['glat'],
           src1['ts'],src1['Npred'],
           src0['ts'],src1['ts'],src2['ts'],
           src0['offset'],src1['offset'],src2['offset'],
           src1['dfde1000'][0],src1['dfde1000'][1],
           src1['dfde1000_index'][0],src1['dfde1000_index'][1],
           src1['flux1000'][0],src0['flux1000'][1],
           src1['eflux1000'][0],src0['eflux1000'][1],
           src1['flux100'][0],src0['flux100'][1],
           src1['eflux100'][0],src0['eflux100'][1],
           src1['SpectrumType']]

    row += ext0_results
    row += ext1_results
    row += ext2_results    
    row += [fit1_dlike,fit2_dlike,nsrc]

    codename = os.path.basename(d)
    linkname = '{%s}'%codename.replace('+','p').replace('.','_')
    
    for hd in halo_data1:
        colname = 'fit1_halo_%.2f_%.3f_ts'%(np.abs(hd['params']['Index'][0]),
                                            hd['SpatialWidth'])
        if not colname in tab.colnames:
            tab.add_column(Column([], name=colname, dtype='f8', format='%.3f'))        
        row += [hd['ts']]

        colname = 'fit1_halo_%.2f_%.3f_eflux_ul95'%(np.abs(hd['params']['Index'][0]),
                                                    hd['SpatialWidth'])
        if not colname in tab.colnames:
            tab.add_column(Column([], name=colname, dtype='f8', format='%.3f'))
        
        row += [hd['eflux_ul95']]

    for hd in halo_data2:
        colname = 'fit2_halo_%.2f_%.3f_ts'%(np.abs(hd['params']['Index'][0]),
                                            hd['SpatialWidth'])
        if not colname in tab.colnames:
            tab.add_column(Column([], name=colname, dtype='f8', format='%.3f'))        
        row += [hd['ts']]

        colname = 'fit2_halo_%.2f_%.3f_eflux_ul95'%(np.abs(hd['params']['Index'][0]),
                                                    hd['SpatialWidth'])
        if not colname in tab.colnames:
            tab.add_column(Column([], name=colname, dtype='f8', format='%.3f'))
        
        row += [hd['eflux_ul95']]
            
    tab.add_row([src0['name'],codename,linkname] + row)

    
m = tab['class']==''
tab['class'][m] = 'unkn'
    
tab.write(args.output,format='fits',overwrite=True)
