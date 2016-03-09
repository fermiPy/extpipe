import numpy as np
import glob
import sys
import os
import traceback
import argparse

from astropy.table import Table, Column
from fermipy.roi_model import ROIModel
from haloanalysis.utils import collect_dirs
from haloanalysis.batch import *

def find_source(data):

    srcs = sorted(data['sources'].values(), key=lambda t: t['offset'])

    for s in srcs:
        if s['SourceType'] == 'DiffuseSource':
            continue

        return s

def find_new_sources(data):

    srcs = sorted(data['sources'].values(), key=lambda t: t['offset'])

    new_srcs = []
    
    for s in srcs:
        if s['SourceType'] == 'DiffuseSource':
            continue

        if s['offset'] > 1.0:
            continue

        if '3FGL' in s['name']:
            continue
        
        new_srcs += [s]

    return new_srcs
    
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

        Column([], name='ext_ts', dtype='f8', format='%.2f'),
        Column([], name='ext_mle', dtype='f8', format='%.3f',unit='deg'),
        Column([], name='ext_err', dtype='f8', format='%.3f',unit='deg'),
        Column([], name='ext_ul95', dtype='f8', format='%.3f',unit='deg',description='Extension 95% CL UL (post-localization)'),
        
        Column([], name='fit1_dlike', dtype='f8', format='%.2f'),
        Column([], name='fit2_dlike', dtype='f8', format='%.2f'),
        Column([], name='fit3_dlike', dtype='f8', format='%.2f'),
        Column([], name='fit4_dlike', dtype='f8', format='%.2f'),
        Column([], name='fit_dlike', dtype='f8', format='%.2f'),
        Column([], name='fit2_nsrc', dtype='i8'),
        Column([], name='fit_nsrc', dtype='i8'),
        ]


tab = Table(cols)

for d in dirs:

    logfile = os.path.join(d,'run_analysis.log')

    if not os.path.isfile(logfile):
        continue

    if not check_log(logfile)=='Successful':
        continue
    
    print d

    fit_data = []
    halo_data = []

    for i in range(5):

        file1 = os.path.join(d,'fit%i.npy'%i)
        file2 = os.path.join(d,'fit%i_halo_data.npy'%i)
        if os.path.isfile(file1):
            fit_data += [np.load(file1).flat[0]]
        
        if os.path.isfile(file2):
            halo_data += [np.load(file2)]

#    try:
#        data0 = np.load(file0).flat[0]
#        data1 = np.load(file1).flat[0]
#        data2 = np.load(file2).flat[0]
#        halo_data1 = np.load(file3)
#        halo_data2 = np.load(file4)
#    except Exception as e:
#        traceback.print_exc()
#        continue

    src_data = []
    ext_data = []
    ext_results = []
    new_srcs = []
    
    for fd in fit_data:

        src = find_source(fd)
        srcs = find_new_sources(fd)
        ext = src['extension']
        src_data += [src]
        ext_data += [src['extension']]
        new_srcs += [srcs]
        
        if src['extension'] is None:
            ext_results += [[np.nan,np.nan,np.nan,np.nan]]
        else:
            ext_results += [[max(ext['ts_ext'],0),ext['ext'],ext['ext_err'],ext['ext_ul95']]]
                        
    fit1_dlike = fit_data[1]['roi']['logLike'] - fit_data[0]['roi']['logLike']
    fit2_dlike = fit_data[2]['roi']['logLike'] - fit_data[0]['roi']['logLike']
    fit_dlike = fit_data[-1]['roi']['logLike'] - fit_data[0]['roi']['logLike']
    
    row = [src_data[0]['assoc']['ASSOC1'],src_data[0]['class'],
           src_data[1]['ra'],src_data[1]['dec'],src_data[1]['glon'],src_data[1]['glat'],
           src_data[1]['ts'],src_data[1]['Npred'],
           src_data[0]['ts'],src_data[1]['ts'],src_data[2]['ts'],
           src_data[0]['offset'],src_data[1]['offset'],src_data[2]['offset'],
           src_data[1]['dfde1000'][0],src_data[1]['dfde1000'][1],
           src_data[1]['dfde1000_index'][0],src_data[1]['dfde1000_index'][1],
           src_data[1]['flux1000'][0],src_data[0]['flux1000'][1],
           src_data[1]['eflux1000'][0],src_data[0]['eflux1000'][1],
           src_data[1]['flux100'][0],src_data[0]['flux100'][1],
           src_data[1]['eflux100'][0],src_data[0]['eflux100'][1],
           src_data[1]['SpectrumType']]

    print len(new_srcs[0]), len(new_srcs[1]), len(new_srcs[2])
    fit2_nsrc = len(new_srcs[2])
    fit_nsrc = len(new_srcs[-1])
    
    for i in range(3):
        row += ext_results[i]

    row += ext_results[-1]
    row += [fit1_dlike,fit2_dlike,np.nan,np.nan,fit_dlike,fit2_nsrc,fit_nsrc]

    
    
    codename = os.path.basename(d)
    linkname = '{%s}'%codename.replace('+','p').replace('.','_')
    
    for hd in halo_data[0]:
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

    for hd in halo_data[1]:
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
            
    tab.add_row([src_data[0]['name'],codename,linkname] + row)

    
m = tab['class']==''
tab['class'][m] = 'unkn'
    
tab.write(args.output,format='fits',overwrite=True)
