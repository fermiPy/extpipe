
import glob
import sys
import os
import argparse
from collections import OrderedDict

import numpy as np
from astropy.table import Table, Column
from haloanalysis.utils import collect_dirs

def extract_halo_data(halo_data,shape=(4,3)):

    o = dict(index = [],
             width = [],
             ts = [],
             eflux = [],
             dlnl = [],
             dlnl_eflux = [])
    
    for hd in halo_data:

        o['index'] += [np.abs(hd['params']['Index'][0])]
        o['width'] += [hd['SpatialWidth']]
        o['ts'] += [hd['ts']]
        o['eflux'] += [hd['eflux_ul95']]
        o['dlnl'] += list(hd['lnlprofile']['dlogLike'])
        o['dlnl_eflux'] += list(hd['lnlprofile']['eflux'])

    for k, v in o.items():
        o[k] = np.array(o[k])

    o['width'] = o['width'].reshape(shape)
    o['index'] = o['index'].reshape(shape)
        
    o['ts'] = o['ts'].reshape(shape)
    o['eflux'] = o['eflux'].reshape(shape)
    o['dlnl'] = o['dlnl'].reshape(shape + (9,))
    o['dlnl_eflux'] = o['dlnl_eflux'].reshape(shape + (9,))
        
    return o

if __name__ == '__main__':

    usage = "usage: %(prog)s"
    description = "Aggregate analysis output."
    parser = argparse.ArgumentParser(usage=usage,description=description)


    parser.add_argument('--output', default = 'test.fits')
    parser.add_argument('dirs', nargs='+', default = None,
                        help='Run analyses in all subdirectories of this '
                        'directory.')

    args = parser.parse_args()
    dirs = sorted(collect_dirs(args.dirs))

    cols_dict = OrderedDict()

    cols_dict['name'] = dict(dtype='S30', format='%s',description='Source Name')
    cols_dict['cfgname'] = dict(dtype='S30', format='%s',description='Source Name')
    cols_dict['ext_ul95'] = dict(dtype='f8', format='%.3f')
    cols_dict['halo_eflux_ul95'] = dict(dtype='f8', format='%.4g',shape=(4,3))
    cols_dict['index'] = dict(dtype='f8', format='%.3f')
    cols_dict['amp'] = dict(dtype='f8', format='%.3f')
    cols_dict['ts'] = dict(dtype='f8', format='%.3f')

    cols = [Column(name=k, **v) for k,v in cols_dict.items()]

    tab = Table(cols)



    for d in dirs:

        print d, os.path.basename(d)
        cfgname = d.split('/')[0]

        file0 = os.path.join(d,'base_model.npy')
        file1 = os.path.join(d,'ext_fit_data.npy')

        if not os.path.isfile(file0) or not os.path.isfile(file1):
            continue


        data0 = np.load(file0).flat[0]
        data1 = np.load(file1)

        halo_files = glob.glob(os.path.join(d,'fit*halo_data.npy'))

        if len(halo_files) == 0:
            continue
        
        halo_data = []
        for f in sorted(halo_files):

            print f

            halo_data += [np.load(f)]

        src0 = data0['sources']['testsource']

        row_dict = {}

        for i, s in enumerate(data1):

            halo_pars = extract_halo_data(halo_data[i])
            
            name = '%03.1f_%03.1f'%(-src0['params']['Index'][0],
                                    np.log10(src0['params']['Prefactor'][0]/1E-14))

            ext = s['extension']

            row_dict['name'] = name
            row_dict['cfgname'] = cfgname
            row_dict['index'] = -src0['params']['Index'][0]
            row_dict['amp'] = np.log10(src0['params']['Prefactor'][0]/1E-14)
            row_dict['eflux'] = src0['eflux']
            row_dict['ts'] = s['ts']
            row_dict['ext_ul95'] = ext['ext_ul95']
            row_dict['halo_eflux_ul95'] = halo_pars['eflux']
            #row_dict['halo_eflux_ul95'] = np.zeros((4,1))

            tab.add_row([row_dict[k] for k in cols_dict.keys()])

            #row = [name,s['ra'],s['dec'],s['glon'],s['glat'],
            #       -src0['params']['Index'][0],src0['params']['Prefactor'][0],
            #       src0['eflux'][0],s['ts'],s['Npred'],
            #       max(ext['ts_ext'],0),ext['ext'],ext['ext_err'],ext['ext_ul95']]
            #tab.add_row(row)

    tab.write(args.output,format='fits',overwrite=True)
