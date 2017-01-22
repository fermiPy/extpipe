import numpy as np
import glob
import sys
import os
import traceback
import argparse
from collections import OrderedDict

from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.table import Table, Column

from fermipy.utils import collect_dirs
from haloanalysis.batch import *



def extract_sources(tab_srcs):

    # Add all sources
    srcs_delete = []

    skydir = SkyCoord(tab_srcs['ra'],tab_srcs['dec'],unit='deg',
                      frame='icrs')

    srcs = []

    for k, v in fit_data[-1]['sources'].items():

        if v['name'] == 'galdiff':
            continue
        if v['name'] == 'isodiff':
            continue

        skydir0 = SkyCoord(v['ra'],v['dec'],unit='deg',
                           frame='icrs')

        sep = skydir0.separation(skydir).deg

        row_dict_src = {}
        # Look for source of the same name
        m = (tab_srcs['name'] == v['name'])

        if np.sum(m):

            mrow = tab_srcs[m][0]

            print 'match name ', v['name'], mrow['name'], np.min(sep)

            if v['offset'] < mrow['offset']:
                srcs_delete += [np.argmax(m)]
            else:
                continue
        elif len(sep) and np.min(sep) < 0.1:

            mrow = tab_srcs[np.argmin(sep)]

            print 'match sep ', v['name'], mrow['name'], np.min(sep)

            if v['offset'] < mrow['offset']:
                srcs_delete += [np.argmin(sep)]
            else:
                continue

        row_dict_src['name'] = v['name']

        if 'assoc' in v and v['assoc']:
            row_dict_src['assoc'] = v['assoc']['ASSOC1']
        else:
            row_dict_src['assoc'] = ''

        if v['class'] is None or v['class'] == '':
            row_dict_src['class'] = 'unkn'
        else:
            row_dict_src['class'] = v['class']

        row_dict_src['ra'] = v['ra']
        row_dict_src['dec'] = v['dec']
        row_dict_src['glon'] = v['glon']
        row_dict_src['glat'] = v['glat']
        row_dict_src['ts'] = v['ts']
        row_dict_src['npred'] = v['npred']
        row_dict_src['offset'] = v['offset']
        row_dict_src['offset_roi'] = max(np.abs(v['offset_glon']),np.abs(v['offset_glat']))

        if '3FGL' in v['name']:
            row_dict_src['in3fgl'] = 1
        else:
            row_dict_src['in3fgl'] = 0

        srcs += [row_dict_src]

    # remove source pairs
    print 'Deleting rows ', srcs_delete
    tab_srcs.remove_rows(srcs_delete)

    for row in srcs:
        tab_srcs.add_row([row[k] for k in cols_dict_srcs.keys()])

def extract_halo_sed(halo_data,halo_scan_shape):

    o = dict(dlnl = [], loglike = [],
             dlnl_eflux = [])

    if len(halo_data) == 0:
        return {}
    
    for hd in halo_data:

        eflux_scan = hd['sed']['norm_scan']*hd['sed']['ref_eflux'][:,np.newaxis]    
        o['dlnl'] += [hd['sed']['dloglike_scan']]
        o['loglike'] += [hd['sed']['loglike_scan']]
        o['dlnl_eflux'] += [eflux_scan]

    o['dlnl'] = np.stack(o['dlnl'])
    o['loglike'] = np.stack(o['loglike'])
    #o['dlnl'] = o['loglike'] - np.max(np.max(o['loglike'],axis=2),axis=0)[None,:,None]
    #o['dlnl'] = o['loglike'] - o['loglike']
    o['dlnl'] = o['loglike'] - o['loglike'][:,:,0][:,:,None]
    o['dlnl_eflux'] = np.stack(o['dlnl_eflux'])        
    return o

def extract_halo_data(halo_data, halo_scan_shape):
    
    o = dict(index = [],
             width = [],
             ts = [],
             eflux = [],
             eflux_ul95 = [],
             dlnl = [],
             loglike = [],
             dlnl_eflux = [])

    if halo_data.ndim == 0 or len(halo_data) == 0:
        return {}
    
    for hd in halo_data:

        o['index'] += [np.abs(hd['params']['Index'][0])]
        o['width'] += [hd['SpatialWidth']]

        ts = 2.0*(np.max(hd['lnlprofile']['loglike']) - hd['lnlprofile']['loglike'][0])        
        o['ts'] += [ts]
#        print '%10.3f %10.3f %10.3f %10.3f'%(np.abs(hd['params']['Index'][0]),
#                                             hd['SpatialWidth'], ts, hd['ts'])

        o['eflux'] += [hd['eflux'][0]]
        o['eflux_ul95'] += [hd['eflux_ul95']]
        o['dlnl'] += list(hd['lnlprofile']['dloglike'])
        o['loglike'] += list(hd['lnlprofile']['loglike'])
        o['dlnl_eflux'] += list(hd['lnlprofile']['eflux'])

    for k, v in o.items():
        o[k] = np.array(o[k])

    o['width'] = o['width'].reshape(halo_scan_shape)
    o['index'] = o['index'].reshape(halo_scan_shape)
        
    o['ts'] = o['ts'].reshape(halo_scan_shape)
    o['eflux'] = o['eflux'].reshape(halo_scan_shape)
    o['eflux_ul95'] = o['eflux_ul95'].reshape(halo_scan_shape)
    o['dlnl'] = o['dlnl'].reshape(halo_scan_shape + (9,))
    o['loglike'] = o['loglike'].reshape(halo_scan_shape + (9,))
    o['dlnl_eflux'] = o['dlnl_eflux'].reshape(halo_scan_shape + (9,))
    o['dlnl'] =  o['loglike'] - o['loglike'][:,:,0][:,:,None]
    
    return o

def centroid(x,y):

    return np.sum(x)/len(x), np.sum(y)/len(y)

def mean_separation(x,y,xc,yc):

    return np.sum(((x-xc)**2+(y-yc)**2)/len(x))**0.5

def find_source(data):

    srcs = [s for s in data['sources'].values() if np.isfinite(s['offset'])]        
    srcs = sorted(srcs, key=lambda t: float(t['offset']))
    return srcs[0]


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


def create_tables(nebins, next_bins, halo_scan_shape, eflux_scan_pts):

    nfit = 6
    nscan_pts = 9
    
    cols_dict_sed = OrderedDict()
    cols_dict_sed['NAME'] = dict(dtype='S20', format='%s',description='Source Name')
    cols_dict_sed['NORM'] = dict(dtype='f8', format='%.3f',shape=(nebins))
    cols_dict_sed['NORM_ERR'] = dict(dtype='f8', format='%.3f',shape=(nebins))
    cols_dict_sed['NORM_ERRP'] = dict(dtype='f8', format='%.3f',shape=(nebins))
    cols_dict_sed['NORM_ERRN'] = dict(dtype='f8', format='%.3f',shape=(nebins))
    cols_dict_sed['NORM_UL'] = dict(dtype='f8', format='%.3f',shape=(nebins))
    cols_dict_sed['REF_FLUX'] = dict(dtype='f8', format='%.3f',shape=(nebins))
    cols_dict_sed['REF_EFLUX'] = dict(dtype='f8', format='%.3f',shape=(nebins))
    cols_dict_sed['REF_NPRED'] = dict(dtype='f8', format='%.3f',shape=(nebins))
    cols_dict_sed['REF_DFDE'] = dict(dtype='f8', format='%.3f',shape=(nebins), unit='1 / (MeV cm2 s)')
    cols_dict_sed['NORM_SCAN'] = dict(dtype='f8', format='%.3f',shape=(nebins,nscan_pts))
    cols_dict_sed['DLOGLIKE_SCAN'] = dict(dtype='f8', format='%.3f',shape=(nebins,nscan_pts))

    
    cols_dict_lnl = OrderedDict()
    cols_dict_lnl['name'] = dict(dtype='S20', format='%s',description='Source Name')

    cols_dict_lnl['fit_ext_scan_dlnl'] = dict(dtype='f8', format='%.3f',shape=(next_bins))
    cols_dict_lnl['fit1_ext_scan_dlnl'] = dict(dtype='f8', format='%.3f',shape=(next_bins))

    cols_dict_lnl['fit_halo_scan_dlnl'] = dict(dtype='f8', format='%.3f',
                                               shape=halo_scan_shape + (len(eflux_scan_pts),))
    cols_dict_lnl['fit1_halo_scan_dlnl'] = dict(dtype='f8', format='%.3f',
                                                shape=halo_scan_shape + (len(eflux_scan_pts),))
    
    cols_dict_lnl['fit_halo_sed_scan_dlnl'] = dict(dtype='f8', format='%.3f',
                                                   shape=(halo_scan_shape[0],nebins,len(eflux_scan_pts)))

    
    cols_dict_lnl['fit_src_sed_scan_dlnl'] = dict(dtype='f8', format='%.3f',
                                                   shape=(nebins,nscan_pts))
    cols_dict_lnl['fit_src_sed_scan_eflux'] = dict(dtype='f8', format='%.4g',
                                                   shape=(nebins,nscan_pts))
    
    cols_dict = OrderedDict()
    cols_dict['name'] = dict(dtype='S20', format='%s',description='Source Name')
    cols_dict['codename'] = dict(dtype='S20', format='%s')
    cols_dict['linkname'] = dict(dtype='S20', format='%s')
    cols_dict['assoc'] = dict(dtype='S20', format='%s')
    cols_dict['class'] = dict(dtype='S20', format='%s')
    cols_dict['ra'] = dict(dtype='f8', format='%.3f',unit='deg')
    cols_dict['dec'] = dict(dtype='f8', format='%.3f',unit='deg')
    cols_dict['glon'] = dict(dtype='f8', format='%.3f',unit='deg')
    cols_dict['glat'] = dict(dtype='f8', format='%.3f',unit='deg')
    cols_dict['ts'] = dict(dtype='f8', format='%.2f')
    cols_dict['npred'] = dict(dtype='f8', format='%.2f')
    cols_dict['fitn_ts'] = dict(dtype='f8', format='%.2f',shape=(nfit,))
    cols_dict['fit_ts'] = dict(dtype='f8', format='%.2f')
    cols_dict['fitn_offset'] = dict(dtype='f8', format='%.2f',shape=(nfit,))
    cols_dict['fit_offset'] = dict(dtype='f8', format='%.2f')
    cols_dict['fitn_ra'] = dict(dtype='f8', format='%.2f',shape=(nfit,))
    cols_dict['fit_ra'] = dict(dtype='f8', format='%.2f')
    cols_dict['fitn_dec'] = dict(dtype='f8', format='%.2f',shape=(nfit,))
    cols_dict['fit_dec'] = dict(dtype='f8', format='%.2f')

    cols_dict['dfde10000'] = dict(dtype='f8', format='%.4g')
    cols_dict['dfde10000_err'] = dict(dtype='f8', format='%.4g')
    cols_dict['flux10000'] = dict(dtype='f8', format='%.4g')
    cols_dict['flux10000_err'] = dict(dtype='f8', format='%.4g')
    cols_dict['eflux10000'] = dict(dtype='f8', format='%.4g')
    cols_dict['eflux10000_err'] = dict(dtype='f8', format='%.4g')    
    cols_dict['dfde1000'] = dict(dtype='f8', format='%.4g')
    cols_dict['dfde1000_err'] = dict(dtype='f8', format='%.4g')
    cols_dict['dfde1000_index'] = dict(dtype='f8', format='%.2f')
    cols_dict['dfde1000_index_err'] = dict(dtype='f8', format='%.2f')
    cols_dict['flux1000'] = dict(dtype='f8', format='%.4g')
    cols_dict['flux1000_err'] = dict(dtype='f8', format='%.4g')
    cols_dict['eflux1000'] = dict(dtype='f8', format='%.4g')
    cols_dict['eflux1000_err'] = dict(dtype='f8', format='%.4g')
    cols_dict['dfde100'] = dict(dtype='f8', format='%.4g')
    cols_dict['dfde100_err'] = dict(dtype='f8', format='%.4g')
    cols_dict['flux100'] = dict(dtype='f8', format='%.4g')
    cols_dict['flux100_err'] = dict(dtype='f8', format='%.4g')
    cols_dict['eflux100'] = dict(dtype='f8', format='%.4g')
    cols_dict['eflux100_err'] = dict(dtype='f8', format='%.4g')

    cols_dict['spectrum_type'] = dict(dtype='S20', format='%s')

    # Halo Fits
    cols_dict['fitn_halo_ts'] = dict(dtype='f8', format='%.2f',shape=(nfit,))
    cols_dict['fitn_halo_eflux'] = dict(dtype='f8', format='%.2f',shape=(nfit,))
    cols_dict['fitn_halo_eflux_ul95'] = dict(dtype='f8', format='%.2f',shape=(nfit,))
    cols_dict['fitn_halo_width'] = dict(dtype='f8', format='%.2f',shape=(nfit,))
    cols_dict['fitn_halo_index'] = dict(dtype='f8', format='%.2f',shape=(nfit,))

    cols_dict['fit_halo_ts'] = dict(dtype='f8', format='%.2f')
    cols_dict['fit_halo_eflux'] = dict(dtype='f8', format='%.2f')
    cols_dict['fit_halo_eflux_ul95'] = dict(dtype='f8', format='%.2f')
    cols_dict['fit_halo_width'] = dict(dtype='f8', format='%.2f')
    cols_dict['fit_halo_index'] = dict(dtype='f8', format='%.2f')

    cols_dict['fitm_halo_ts'] = dict(dtype='f8', format='%.2f')
    cols_dict['fitm_halo_eflux'] = dict(dtype='f8', format='%.2f')
    cols_dict['fitm_halo_eflux_ul95'] = dict(dtype='f8', format='%.2f')
    cols_dict['fitm_halo_width'] = dict(dtype='f8', format='%.2f')
    cols_dict['fitm_halo_index'] = dict(dtype='f8', format='%.2f')

    # Ext Fits
    cols_dict['fitn_ext_ts'] = dict(dtype='f8', format='%.2f',shape=(nfit,))
    cols_dict['fitn_ext_mle'] = dict(dtype='f8', format='%.2f',shape=(nfit,))
    cols_dict['fitn_ext_err'] = dict(dtype='f8', format='%.2f',shape=(nfit,))
    cols_dict['fitn_ext_ul95'] = dict(dtype='f8', format='%.2f',shape=(nfit,))

    cols_dict['fit_ext_ts'] = dict(dtype='f8', format='%.2f')
    cols_dict['fit_ext_mle'] = dict(dtype='f8', format='%.2f')
    cols_dict['fit_ext_err'] = dict(dtype='f8', format='%.2f')
    cols_dict['fit_ext_ul95'] = dict(dtype='f8', format='%.2f')

    cols_dict['fitm_ext_ts'] = dict(dtype='f8', format='%.2f')
    cols_dict['fitm_ext_mle'] = dict(dtype='f8', format='%.2f')
    cols_dict['fitm_ext_err'] = dict(dtype='f8', format='%.2f')
    cols_dict['fitm_ext_ul95'] = dict(dtype='f8', format='%.2f')

    # Delta LogLikes
    cols_dict['fitn_dlike'] = dict(dtype='f8', format='%.2f',shape=(nfit,))
    cols_dict['fitn_dlike1'] = dict(dtype='f8', format='%.2f',shape=(nfit,))
    cols_dict['fitn_dlike_ps_ext'] = dict(dtype='f8', format='%.2f',shape=(nfit,))
    cols_dict['fitn_dlike_ps_halo'] = dict(dtype='f8', format='%.2f',shape=(nfit,))
    cols_dict['fitn_dlike1_ext'] = dict(dtype='f8', format='%.2f',shape=(nfit,))
    cols_dict['fitn_dlike1_halo'] = dict(dtype='f8', format='%.2f',shape=(nfit,))    

    cols_dict['fitm_dlike'] = dict(dtype='f8', format='%.2f')
    cols_dict['fitm_dlike1'] = dict(dtype='f8', format='%.2f')
    cols_dict['fitm_dlike_ps_ext'] = dict(dtype='f8', format='%.2f')
    cols_dict['fitm_dlike_ps_halo'] = dict(dtype='f8', format='%.2f')
    cols_dict['fitm_dlike1_ext'] = dict(dtype='f8', format='%.2f')
    cols_dict['fitm_dlike1_halo'] = dict(dtype='f8', format='%.2f')
    
    cols_dict['fit_dlike_ext'] = dict(dtype='f8', format='%.2f')
    cols_dict['fit_dlike_halo'] = dict(dtype='f8', format='%.2f')
    cols_dict['fit_dlike_ps_ext'] = dict(dtype='f8', format='%.2f')
    cols_dict['fit_dlike_ps_halo'] = dict(dtype='f8', format='%.2f')
    cols_dict['fit_daic_ps_ext'] = dict(dtype='f8', format='%.2f')
    cols_dict['fit_daic_ps_halo'] = dict(dtype='f8', format='%.2f')
    cols_dict['fit_dlike1_ext'] = dict(dtype='f8', format='%.2f')
    cols_dict['fit_dlike1_halo'] = dict(dtype='f8', format='%.2f')
    
    cols_dict['fit1_halo_scan_ts'] = dict(dtype='f8', format='%.2f',shape=halo_scan_shape)
    cols_dict['fit1_halo_scan_eflux_ul95'] = dict(dtype='f8', format='%.4g',shape=halo_scan_shape)

    cols_dict['fit_halo_scan_ts'] = dict(dtype='f8', format='%.2f',shape=halo_scan_shape)
    cols_dict['fit_halo_scan_eflux_ul95'] = dict(dtype='f8', format='%.4g',shape=halo_scan_shape)

    cols_dict['fit_nsrc'] = dict(dtype='i8')
    cols_dict['fit_nsrc_ext'] = dict(dtype='i8')
    cols_dict['fit_nsrc_halo'] = dict(dtype='i8')
    cols_dict['fit_mean_sep'] = dict(dtype='f8', format='%.3f')

    tab = Table([Column(name=k, **v) for k,v in cols_dict.items()])
    tab_lnl = Table([Column(name=k, **v) for k,v in cols_dict_lnl.items()])
    tab_sed = Table([Column(name=k, **v) for k,v in cols_dict_sed.items()])

    return tab, tab_lnl, tab_sed

    
def aggregate(dirs,output,suffix=''):

    dirs = [d for argdir in dirs for d in collect_dirs(argdir)]

    nfit = 6
    nscan_pts = 9
    nebins = 20
    next_bins = 30
    halo_scan_shape = (13,9)
    eflux_scan_pts = np.logspace(-9,-5,41)
    
    row_dict = {}
    row_dict_lnl = {}
    row_dict_sed = {}
    tab = None
    tab_lnl = None
    tab_sed = None
    
    halo_name = 'halo_RadialGaussian'
    
    for d in dirs:

        logfile0 = os.path.join(d,'run-region-analysis.log')
        logfile1 = os.path.join(d,'run-halo-analysis.log')

        if not os.path.isfile(logfile0):
            continue

        if not check_log(logfile0)=='Successful':
            continue

        if os.path.isfile(logfile1) and check_log(logfile1)=='Exited':
            continue
        
        print d

        fit_data = []
        halo_data = []
        halo_data_sed = []
        halo_fit = []
        new_src_data = []

        if not os.path.isfile(os.path.join(d,'new_source_data.npy')):
            print 'skipping'
            continue

        new_src_data = np.load(os.path.join(d,'new_source_data.npy'))

        for i in range(nfit):

            file1 = os.path.join(d,'fit%i%s.npy'%(i,suffix))
            file2 = os.path.join(d,'fit%i%s_%s_data.npy'%(i,suffix,halo_name))
            file3 = os.path.join(d,'fit%i%s_%s.npy'%(i,suffix,halo_name))
            file4 = os.path.join(d,'fit%i%s_%s_data_idx_free.npy'%(i,suffix,halo_name))
            
            if os.path.isfile(file1):
                fit_data += [np.load(file1).flat[0]]

            if os.path.isfile(file2):
                halo_data += [np.load(file2)]

            if os.path.isfile(file4):
                halo_data_sed += [np.load(file4)]

            if os.path.isfile(file3):
                halo_fit += [np.load(file3).flat[0]]

        if len(fit_data) == 0:
            print 'skipping'
            continue

        src_name = fit_data[0]['config']['selection']['target']
        src = fit_data[0]['sources'][src_name]

        

        src_data = []
        ext_data = []
        new_srcs = []

        src_ts = [np.nan]*nfit
        src_offset = [np.nan]*nfit
        fit_dlike = [np.nan]*nfit
        fit_dlike_ext = [np.nan]*nfit
        fit_dlike_halo = [np.nan]*nfit
        fit_dlike_ps_ext = [np.nan]*nfit
        fit_dlike_ps_halo = [np.nan]*nfit
        fit_daic_ps_ext = [np.nan]*nfit
        fit_daic_ps_halo = [np.nan]*nfit
        
        fit_dlike1 = [np.nan]*nfit
        fit_dlike1_ext = [np.nan]*nfit
        fit_dlike1_halo = [np.nan]*nfit

        extn_ts = [np.nan]*nfit
        extn_mle = [np.nan]*nfit
        extn_err = [np.nan]*nfit
        extn_ul95 = [np.nan]*nfit

        halon_ts = [np.nan]*nfit
        halon_width = [np.nan]*nfit
        halon_index = [np.nan]*nfit
        halon_eflux = [np.nan]*nfit
        halon_eflux_ul95 = [np.nan]*nfit
        halon_loglike = [np.nan]*nfit
        
        halo_pars = []
        for hd in halo_data:
            halo_pars += [extract_halo_data(hd,halo_scan_shape)]

        halo_seds = []
        for hd in halo_data_sed:
            halo_seds += [extract_halo_sed(hd,halo_scan_shape)]
            
        ext_width = []
            
        for i, fd in enumerate(fit_data):

            src = fd['sources'][src_name]
            ext = src['extension']
            src_data += [src]
            ext_data += [src['extension']]

            if src['extension'] is not None:
                extn_ts[i] = max(ext['ts_ext'],0)
                extn_mle[i] = ext['ext']
                extn_err[i] = ext['ext_err']
                extn_ul95[i] = ext['ext_ul95']
                ext_width += [ext['width']]

        for i, fd in enumerate(halo_fit):
            continue
            halon_ts[i+1] = fd['sources'][halo_name]['ts']
            halon_eflux_ul95[i+1] = fd['sources'][halo_name]['eflux_ul95']
            halon_width[i+1] = fd['sources'][halo_name]['SpatialWidth']
            halon_index[i+1] = np.abs(fd['sources'][halo_name]['params']['Index'][0])            

        for i in range(len(halo_pars)):

            if not halo_pars[0]:
                continue
            
            halo_scan_ts = np.array(halo_pars[i]['ts']).reshape(halo_scan_shape)
            halo_scan_loglike = np.max(np.array(halo_pars[i]['loglike']),axis=2)
            halo_scan_eflux = np.array(halo_pars[i]['eflux']).reshape(halo_scan_shape)
            halo_scan_eflux_ul95 = np.array(halo_pars[i]['eflux_ul95']).reshape(halo_scan_shape)
            halo_scan_idx = np.argmax(halo_scan_ts)
            halo_scan_max_ts = halo_scan_ts.flat[halo_scan_idx]
            halo_scan_max_eflux = halo_scan_eflux.flat[halo_scan_idx]
            halo_scan_max_eflux_ul95 = halo_scan_eflux_ul95.flat[halo_scan_idx]
            halo_scan_max_width = halo_pars[0]['width'].flat[halo_scan_idx]
            halo_scan_max_index = halo_pars[0]['index'].flat[halo_scan_idx] 

            halon_ts[i+1] = halo_scan_max_ts
            halon_eflux[i+1] = halo_scan_max_eflux
            halon_eflux_ul95[i+1] = halo_scan_max_eflux_ul95
            halon_width[i+1] = halo_scan_max_width
            halon_index[i+1] = halo_scan_max_index
            halon_loglike[i+1] = halo_scan_loglike.flat[halo_scan_idx]
            
#            print i, halo_scan_max_ts, halo_scan_max_width, halo_scan_max_index
#            print i, halon_ts[i+1], halon_width[i+1], halon_index[i+1]
            
        for i in range(len(src_data)):
            src_ts[i] = src_data[i]['ts']
            src_offset[i] = src_data[i]['offset']
            
        fit_nsrc = len(fit_data)-1
        for i in range(len(fit_data)):

            loglike = fit_data[i]['roi']['loglike']
#            loglike0 = fit_data[0]['roi']['loglike']
            loglike1 = fit_data[1]['roi']['loglike']

            if i >= 1:
                fit_dlike[i] = loglike - fit_data[i-1]['roi']['loglike']
                fit_dlike1[i] = loglike - loglike1
                fit_dlike1_ext[i] = (loglike+extn_ts[i]/2.) - loglike1
                fit_dlike1_halo[i] = (loglike+halon_ts[i]/2.) - loglike1
                fit_dlike_ext[i] = fit_dlike1_ext[i] - fit_dlike1_ext[i-1]
                fit_dlike_halo[i] = fit_dlike1_halo[i] - fit_dlike1_halo[i-1]

            # L_(N+1) - L_(ext)
            if i >= 1 and i < len(fit_data)-1:
                fit_dlike_ps_ext[i] = fit_data[i+1]['roi']['loglike'] - (fit_data[i]['roi']['loglike']+extn_ts[i]/2.)
                fit_dlike_ps_halo[i] = fit_data[i+1]['roi']['loglike'] - (fit_data[i]['roi']['loglike']+halon_ts[i]/2.)
                fit_daic_ps_ext[i] = -2.0*(fit_dlike_ps_ext[i] - 3.0)
                fit_daic_ps_halo[i] = -2.0*(fit_dlike_ps_halo[i] - 1.0)
                
        print '-'*80
        fit_nsrc_ext = fit_nsrc
        fit_nsrc_halo = fit_nsrc
        for i in range(1,len(fit_data)-1):

            delta_ext_ts = extn_ts[i+1] - extn_ts[i]
            print '%4i %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f'%(i, fit_dlike1[i], fit_dlike1_ext[i], 0.5*extn_ts[i], fit_dlike[i],
                                                                   fit_dlike_ext[i], fit_dlike_ps_ext[i],delta_ext_ts)
            if fit_dlike_ps_ext[i] < 3.0 and delta_ext_ts < -1.0 \
                    and fit_nsrc_ext == fit_nsrc:
                fit_nsrc_ext = i

        for i in range(1,len(fit_data)-1):
            
            delta_halo_ts = halon_ts[i+1] - halon_ts[i]
            print '%4i %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f'%(i, fit_dlike1[i], fit_dlike1_halo[i], 0.5*halon_ts[i], fit_dlike[i],
                                                                   fit_dlike_halo[i], fit_dlike_ps_halo[i],delta_halo_ts)
            if fit_dlike_ps_halo[i] < 1.0 and delta_halo_ts < -1.0 \
                    and fit_nsrc_halo == fit_nsrc:
                fit_nsrc_halo = i
            
        x = [src_data[-1]['offset_glon']]
        y = [src_data[-1]['offset_glat']]

        for j in range(len(new_src_data)):
            x += [new_src_data[j]['offset_glon']]
            y += [new_src_data[j]['offset_glat']]

        xc, yc = centroid(x,y)
        fit_mean_sep = mean_separation(x,y,xc,yc)

        if fit_nsrc == 1:
            fit_mean_sep = np.nan

        codename = os.path.basename(d)
        linkname = '{%s}'%codename.replace('+','p').replace('.','_')

        #src_sed_tab = Table.read('fit%i_%s_sed.fits'%(fit_nsrc,codename))
        
        row_dict['fit1_halo_scan_ts'] = np.nan*np.ones(halo_scan_shape)
        row_dict['fit1_halo_scan_eflux_ul95'] = np.nan*np.ones(halo_scan_shape)
        row_dict['fit_halo_scan_ts'] = np.nan*np.ones(halo_scan_shape)
        row_dict['fit_halo_scan_eflux_ul95'] = np.nan*np.ones(halo_scan_shape)

        if len(halo_pars) and halo_pars[0]:
            row_dict['fit1_halo_scan_ts'] = np.array(halo_pars[0]['ts']).reshape(halo_scan_shape)
            row_dict['fit1_halo_scan_eflux_ul95'] = np.array(halo_pars[0]['eflux']).reshape(halo_scan_shape)
            row_dict['fit_halo_scan_ts'] = np.array(halo_pars[fit_nsrc_halo-1]['ts']).reshape(halo_scan_shape)
            row_dict['fit_halo_scan_eflux_ul95'] = np.array(halo_pars[fit_nsrc_halo-1]['eflux']).reshape(halo_scan_shape)

        fitn_ra = np.array([np.nan]*nfit)
        fitn_dec = np.array([np.nan]*nfit)
        fitn_ra[:len(src_data)] = np.array([sd['ra'] for sd in src_data])
        fitn_dec[:len(src_data)] = np.array([sd['dec'] for sd in src_data])
            
        row_dict['name'] = src_data[0]['name']
        row_dict['codename'] = codename
        row_dict['linkname'] = linkname
        row_dict['assoc'] = src_data[0]['assoc']['ASSOC1']
        row_dict['class'] = src_data[0]['class']
        row_dict['ra'] = src_data[1]['ra']
        row_dict['dec'] = src_data[1]['dec']
        row_dict['glon'] = src_data[1]['glon']
        row_dict['glat'] = src_data[1]['glat']
        row_dict['ts'] = src_data[1]['ts']
        row_dict['npred'] = src_data[1]['npred']
        row_dict['fitn_ts'] = src_ts
        row_dict['fit_ts'] = src_data[-1]['ts']
        row_dict['fitn_offset'] = src_offset
        row_dict['fit_offset'] = src_data[-1]['offset']
        row_dict['fitn_ra'] = fitn_ra
        row_dict['fit_ra'] = src_data[-1]['ra']
        row_dict['fitn_dec'] = fitn_dec
        row_dict['fit_dec'] = src_data[-1]['dec']

        
        for t in ['dfde','flux','eflux']:
            row_dict['%s10000'%t] = src_data[1]['%s10000'%t][0]
            row_dict['%s10000_err'%t] = src_data[1]['%s10000'%t][1]
            row_dict['%s1000'%t] = src_data[1]['%s1000'%t][0]
            row_dict['%s1000_err'%t] = src_data[1]['%s1000'%t][1]
            row_dict['%s100'%t] = src_data[1]['%s100'%t][0]
            row_dict['%s100_err'%t] = src_data[1]['%s100'%t][1]

        row_dict['dfde1000_index'] = src_data[1]['dfde1000_index'][0]
        row_dict['dfde1000_index_err'] = src_data[1]['dfde1000_index'][1]
        row_dict['spectrum_type'] = src_data[1]['SpectrumType']

        row_dict['fitn_dlike'] = fit_dlike
        row_dict['fitn_dlike1'] = fit_dlike1
        row_dict['fitn_dlike_ps_ext'] = fit_dlike_ps_ext
        row_dict['fitn_dlike_ps_halo'] = fit_dlike_ps_halo
        row_dict['fitn_dlike1_ext'] = fit_dlike1_ext
        row_dict['fitn_dlike1_halo'] = fit_dlike1_halo        
        
        row_dict['fitm_dlike'] = fit_dlike[fit_nsrc]
        row_dict['fitm_dlike1'] = fit_dlike1[fit_nsrc]
        row_dict['fitm_dlike_ps_ext'] = fit_dlike_ps_ext[fit_nsrc]
        row_dict['fitm_dlike_ps_halo'] = fit_dlike_ps_halo[fit_nsrc]
        row_dict['fitm_dlike1_ext'] = fit_dlike1_ext[fit_nsrc]
        row_dict['fitm_dlike1_halo'] = fit_dlike1_halo[fit_nsrc]
        
        row_dict['fit_dlike_ext'] = fit_dlike[fit_nsrc_ext]
        row_dict['fit_dlike_halo'] = fit_dlike[fit_nsrc_halo]        
        row_dict['fit_dlike_ps_ext'] = fit_dlike_ps_ext[fit_nsrc_ext]
        row_dict['fit_dlike_ps_halo'] = fit_dlike_ps_halo[fit_nsrc_halo]
        row_dict['fit_daic_ps_ext'] = fit_daic_ps_ext[fit_nsrc_ext]
        row_dict['fit_daic_ps_halo'] = fit_daic_ps_halo[fit_nsrc_halo]
        row_dict['fit_dlike1_ext'] = fit_dlike1_ext[fit_nsrc_ext]
        row_dict['fit_dlike1_halo'] = fit_dlike1_halo[fit_nsrc_halo]
        
        
        row_dict['fit_mean_sep'] = fit_mean_sep
        row_dict['fit_nsrc'] = fit_nsrc
        row_dict['fit_nsrc_halo'] = fit_nsrc_halo
        row_dict['fit_nsrc_ext'] = fit_nsrc_ext

        row_dict['fitn_halo_ts'] = halon_ts
        row_dict['fitn_halo_width'] = halon_width
        row_dict['fitn_halo_index'] = halon_index
        row_dict['fitn_halo_eflux'] = halon_eflux
        row_dict['fitn_halo_eflux_ul95'] = halon_eflux_ul95
        row_dict['fitm_halo_ts'] = halon_ts[fit_nsrc]
        row_dict['fitm_halo_width'] = halon_width[fit_nsrc]
        row_dict['fitm_halo_index'] = halon_index[fit_nsrc]
        row_dict['fitm_halo_eflux'] = halon_eflux[fit_nsrc]
        row_dict['fitm_halo_eflux_ul95'] = halon_eflux_ul95[fit_nsrc]
        row_dict['fit_halo_ts'] = halon_ts[fit_nsrc_halo]
        row_dict['fit_halo_width'] = halon_width[fit_nsrc_halo]
        row_dict['fit_halo_index'] = halon_index[fit_nsrc_halo]
        row_dict['fit_halo_eflux'] = halon_eflux[fit_nsrc_halo]
        row_dict['fit_halo_eflux_ul95'] = halon_eflux_ul95[fit_nsrc_halo]
        
        row_dict['fitn_ext_ts'] = extn_ts
        row_dict['fitn_ext_mle'] = extn_mle
        row_dict['fitn_ext_err'] = extn_err
        row_dict['fitn_ext_ul95'] = extn_ul95
        row_dict['fitm_ext_ts'] = extn_ts[fit_nsrc]
        row_dict['fitm_ext_mle'] = extn_mle[fit_nsrc]
        row_dict['fitm_ext_err'] = extn_err[fit_nsrc]
        row_dict['fitm_ext_ul95'] = extn_ul95[fit_nsrc]    
        row_dict['fit_ext_ts'] = extn_ts[fit_nsrc_ext]
        row_dict['fit_ext_mle'] = extn_mle[fit_nsrc_ext]
        row_dict['fit_ext_err'] = extn_err[fit_nsrc_ext]
        row_dict['fit_ext_ul95'] = extn_ul95[fit_nsrc_ext]        
                
        from scipy.interpolate import UnivariateSpline

        fit1_dlnl_interp = np.zeros(halo_scan_shape + (len(eflux_scan_pts),))
        fit_dlnl_interp = np.zeros(halo_scan_shape + (len(eflux_scan_pts),))

        fit_sed_dlnl_interp = np.zeros((halo_scan_shape[0],nebins,len(eflux_scan_pts)))
                
        for i in range(halo_scan_shape[0]):

            if not len(halo_seds):
                continue
            
            for j in range(nebins):
                sp = UnivariateSpline(halo_seds[fit_nsrc_halo-1]['dlnl_eflux'][i,j],
                                      halo_seds[fit_nsrc_halo-1]['dlnl'][i,j],k=2,s=0.001)
                fit_sed_dlnl_interp[i,j] = sp(eflux_scan_pts)
                
            
            for j in range(halo_scan_shape[1]):
                sp = UnivariateSpline(halo_pars[0]['dlnl_eflux'][i,j],
                                      halo_pars[0]['dlnl'][i,j],k=2,s=0.001)
                fit1_dlnl_interp[i,j] = sp(eflux_scan_pts)

                sp = UnivariateSpline(halo_pars[fit_nsrc_halo-1]['dlnl_eflux'][i,j],
                                      halo_pars[fit_nsrc_halo-1]['dlnl'][i,j],k=2,s=0.001)
                fit_dlnl_interp[i,j] = sp(eflux_scan_pts)

        row_dict_lnl['fit_src_sed_scan_dlnl'] = np.zeros((nebins,nscan_pts))
        row_dict_lnl['fit_src_sed_scan_eflux'] = np.zeros((nebins,nscan_pts))
                
        for i in range(nebins):
            row_dict_lnl['fit_src_sed_scan_dlnl'][i,:] = src_data[fit_nsrc_halo]['sed']['lnlprofile'][i]['dloglike']
            row_dict_lnl['fit_src_sed_scan_eflux'][i,:] = src_data[fit_nsrc_halo]['sed']['lnlprofile'][i]['eflux']

        row_dict_sed['NAME'] = src_data[-1]['name']
        row_dict_sed['NORM'] = src_data[-1]['sed']['norm']
        row_dict_sed['NORM_ERR'] = src_data[-1]['sed']['norm_err']
        row_dict_sed['NORM_ERRP'] = src_data[-1]['sed']['norm_err_hi']
        row_dict_sed['NORM_ERRN'] = src_data[-1]['sed']['norm_err_lo']
        row_dict_sed['NORM_UL'] = src_data[-1]['sed']['norm_ul']
        row_dict_sed['NORM_SCAN'] = src_data[-1]['sed']['norm_scan']
        row_dict_sed['DLOGLIKE_SCAN'] = src_data[-1]['sed']['dloglike_scan']
        row_dict_sed['REF_FLUX'] = src_data[-1]['sed']['ref_flux']
        row_dict_sed['REF_EFLUX'] = src_data[-1]['sed']['ref_eflux']
        row_dict_sed['REF_NPRED'] = src_data[-1]['sed']['ref_npred']
        row_dict_sed['REF_DFDE'] = src_data[-1]['sed']['ref_dfde']
            
        row_dict_lnl['name'] = src_data[0]['name']
        
        row_dict_lnl['fit1_ext_scan_dlnl'] = ext_data[1]['dloglike']
        row_dict_lnl['fit_ext_scan_dlnl'] = ext_data[fit_nsrc_ext]['dloglike']

        row_dict_lnl['fit1_halo_scan_dlnl'] = fit1_dlnl_interp
        row_dict_lnl['fit_halo_scan_dlnl'] = fit_dlnl_interp

        row_dict_lnl['fit_halo_sed_scan_dlnl'] = fit_sed_dlnl_interp

        if tab is None:
            tab, tab_lnl, tab_sed = create_tables(len(src_data[0]['sed']['emin']),
                                                  len(ext_width[0]),
                                                  halo_scan_shape,
                                                  eflux_scan_pts)
        
        tab.add_row([row_dict[k] for k in tab.columns])
        tab_lnl.add_row([row_dict_lnl[k] for k in tab_lnl.columns])
        tab_sed.add_row([row_dict_sed[k] for k in tab_sed.columns])
                

    m = tab['class']==''
    tab['class'][m] = 'unkn'

    output_lnl = os.path.splitext(output)[0] + '_lnl.fits'
    output_sed = os.path.splitext(output)[0] + '_sed.fits'

    cols_dict = OrderedDict()
    if len(halo_pars) and halo_pars[0]:
        cols_dict['halo_scan_width'] = dict(dtype='f8', format='%.3f',
                                            data=halo_pars[0]['width'][:,0][np.newaxis,:])
        cols_dict['halo_scan_index'] = dict(dtype='f8', format='%.3f',
                                            data=halo_pars[0]['index'][0,:][np.newaxis,:])
        cols_dict['halo_scan_eflux'] = dict(dtype='f8', format='%.3f',
                                            data=eflux_scan_pts[np.newaxis,:])
    cols_dict['ext_width'] = dict(dtype='f8', format='%.3f',
                                  data=ext_width[0][np.newaxis,:])
    
    tab_grid = Table([Column(name=k, **v) for k, v in cols_dict.items()])

    cols_dict = OrderedDict()
    cols_dict['E_MIN'] = dict(dtype='f8', format='%.3f',
                              data=src_data[0]['sed']['emin'],unit='MeV')
    cols_dict['E_MAX'] = dict(dtype='f8', format='%.3f',
                              data=src_data[0]['sed']['emax'],unit='MeV')
    cols_dict['E_REF'] = dict(dtype='f8', format='%.3f',
                              data=src_data[0]['sed']['ectr'],unit='MeV')
    cols_dict['REF_FLUX'] = dict(dtype='f8', format='%.3f',
                                 data=src_data[0]['sed']['ref_flux'],
                                 unit='ph / (cm2 s)')
    cols_dict['REF_EFLUX'] = dict(dtype='f8', format='%.3f',
                                  data=src_data[0]['sed']['ref_eflux'],
                                  unit='MeV / (cm2 s)')
    cols_dict['REF_DFDE'] = dict(dtype='f8', format='%.3f',
                                 data=src_data[0]['sed']['ref_dfde'],
                                 unit='1 / (MeV cm2 s)')
        
    tab_ebounds = Table([Column(name=k, **v) for k, v in cols_dict.items()])
    
    hdulist = fits.HDUList()
    hdulist.append(fits.table_to_hdu(tab))
    hdulist.append(fits.table_to_hdu(tab_grid))
    hdulist.append(fits.table_to_hdu(tab_ebounds))
    hdulist[1].name = 'CATALOG'
    hdulist[2].name = 'SCAN_PARS'
    hdulist[3].name = 'EBOUNDS'    
    hdulist.writeto(output,clobber=True)

    hdulist = fits.HDUList()
    hdulist.append(fits.table_to_hdu(tab_lnl))
    hdulist.append(fits.table_to_hdu(tab_grid))
    hdulist.append(fits.table_to_hdu(tab_ebounds))
    hdulist[1].name = 'CATALOG'
    hdulist[2].name = 'SCAN_PARS'
    hdulist[3].name = 'EBOUNDS'    
    hdulist.writeto(output_lnl,clobber=True)

    hdulist = fits.HDUList()
    hdulist.append(fits.table_to_hdu(tab_sed))
    hdulist.append(fits.table_to_hdu(tab_grid))
    hdulist.append(fits.table_to_hdu(tab_ebounds))
    hdulist[1].name = 'CATALOG'
    hdulist[2].name = 'SCAN_PARS'
    hdulist[3].name = 'EBOUNDS'    
    hdulist.writeto(output_sed,clobber=True)
    

def main():

    usage = "usage: %(prog)s"
    description = "Aggregate analysis output."
    parser = argparse.ArgumentParser(usage=usage,description=description)


    parser.add_argument('--output', default = 'test.fits')
    parser.add_argument('--suffix', default = '')


    parser.add_argument('dirs', nargs='+', default = None,
                        help='Run analyses in all subdirectories of this '
                        'directory.')    
    args = parser.parse_args()
    print(args.dirs)
    aggregate(args.dirs,args.output,args.suffix)
    

if __name__ == '__main__':

    main()
