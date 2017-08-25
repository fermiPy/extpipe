import numpy as np
import glob
import sys
import os
import traceback
import argparse
import yaml
from collections import OrderedDict

from scipy.ndimage.interpolation import map_coordinates
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.table import Table, Column

from fermipy.catalog import *
from fermipy.utils import collect_dirs
from fermipy.utils import fit_parabola, get_parameter_limits, create_source_name
from haloanalysis.batch import *


src_colnames = ['ts','npred','offset','ra','dec','glon','glat','index','index_err',
                'pos_err','pos_r68','pos_r95','pos_r99']

ext_colnames = ['ts_ext','sys_ts_ext',
                'r68', 'r68_err', 'r68_sys_err', 'r68_ul95', 'r68_sys_ul95',
                'flux','flux_err','eflux','eflux_err',
                'flux100','flux100_err','eflux100','eflux100_err',
                'flux1000','flux1000_err','eflux1000','eflux1000_err',
                'flux10000','flux10000_err','eflux10000','eflux10000_err',
                'ts','npred','offset',
                'loglike','index','index_err',
                'ra','dec','glon','glat',
                'ra_err','dec_err','glon_err','glat_err', 'pos_err', 'pos_r68', 'pos_r95', 'pos_r99']

halo_colnames = ['ts','r68','index','r68_err','index_err', 'loglike',
                 'eflux','eflux_err', 'eflux_ul95', 'flux','flux_err','flux_ul95']

def match_catalog(row_dict, cat, sigma95_fhes, sigma95_cat, skydir_fhes, skydir_cat, cat_name):

    sigma95 = np.sqrt(sigma95_fhes**2 + sigma95_cat**2)

    sep = skydir_cat.separation(skydir_fhes).deg
    sep_sigma95 = sep/sigma95
    m = (sep < 1.5*sigma95) & np.isfinite(sigma95) 
    
    assoc = np.array(cat['Source_Name'][m])
    assoc = assoc[np.argsort(sep[m])]
    assoc_sep_sigma95 = sep_sigma95[m][np.argsort(sep[m])]
    assoc_sep = np.sort(sep[m])            
    nassoc = min(5,len(assoc))
    row_dict['assoc_%s_num'%cat_name] = nassoc
    row_dict['assoc_%s'%cat_name][:nassoc] = assoc[:nassoc]
    row_dict['assoc_%s_sep'%cat_name][:nassoc] = assoc_sep[:nassoc]
    row_dict['assoc_%s_sep_sigma95'%cat_name][:nassoc] = assoc_sep_sigma95[:nassoc]
    if len(assoc):
        row_dict['name_%s'%cat_name] = assoc[0]
        row_dict['sep_%s'%cat_name] = assoc_sep[0]
        row_dict['sep_%s_sigma95'%cat_name] = assoc_sep_sigma95[0]

    print('nassoc',nassoc)
    print('sigma95_fhes',sigma95_fhes)
    print('sigma95_cat',sigma95_cat[m])
    print('assoc_%s'%cat_name,row_dict['assoc_%s'%cat_name])
    print('assoc_%s_sep'%cat_name,row_dict['assoc_%s_sep'%cat_name])
    print('assoc_%s_sep_sigma95'%cat_name,row_dict['assoc_%s_sep_sigma95'%cat_name])
    print('name_%s'%cat_name,row_dict['name_%s'%cat_name])
    print('sep_%s'%cat_name,row_dict['sep_%s'%cat_name])
    print('sep_%s_sigma95'%cat_name,row_dict['sep_%s_sigma95'%cat_name])


def extract_index(src, egy=1E3):

    if src['SpectrumType'] == 'PowerLaw':
        return src['param_values'][1], src['param_errors'][1]
    elif src['SpectrumType'] == 'LogParabola':

        alpha, alpha_err = src['param_values'][1], src['param_errors'][1]
        beta, beta_err = src['param_values'][2], src['param_errors'][2]
        eb = src['param_values'][3]
        index = -alpha - 2*beta*(np.log(egy) - np.log(eb))
        index_err = alpha_err + 2*beta_err*np.abs(np.log(egy) - np.log(eb))
        return index, index_err
    elif src['SpectrumType'] == 'PLSuperExpCutoff':

        index1, index1_err = src['param_values'][1], src['param_errors'][1]
        cutoff, cutoff_err = src['param_values'][3], src['param_errors'][3]
        index2, index2_err = src['param_values'][4], src['param_errors'][4]

        index = index1 - egy/cutoff
        index_err = index1_err + cutoff_err*(index**2/cutoff*(egy/cutoff)**index2)
        return index, index_err
    else:
        raise Exception

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

        eflux_scan = hd['norm_scan']*hd['ref_eflux'][:,np.newaxis]    
        o['dlnl'] += [hd['dloglike_scan']]
        o['loglike'] += [hd['loglike_scan']]
        o['dlnl_eflux'] += [eflux_scan]

    o['dlnl'] = np.stack(o['dlnl'])
    o['loglike'] = np.stack(o['loglike'])
    #o['dlnl'] = o['loglike'] - np.max(np.max(o['loglike'],axis=2),axis=0)[None,:,None]
    #o['dlnl'] = o['loglike'] - o['loglike']
    o['dlnl'] = o['loglike'] - o['loglike'][:,:,0][:,:,None]
    o['dlnl_eflux'] = np.stack(o['dlnl_eflux'])        
    return o


def extract_ext_data(row_dict, prefix, src_name, ext_data, tab_roi, ext_data_psfhi, ext_data_psflo):

    for i, (ext,ext_roi) in enumerate(zip(ext_data,tab_roi)):

        row_dict['fitn_%s_ts_ext'%prefix][i] = max(ext['ts_ext'],0)
        row_dict['fitn_%s_r68'%prefix][i] = ext['ext']
        row_dict['fitn_%s_r68_err'%prefix][i] = ext['ext_err']
        row_dict['fitn_%s_r68_ul95'%prefix][i] = ext['ext_ul95']            
        row_dict['fitn_%s_ra'%prefix][i] = ext['ra']
        row_dict['fitn_%s_dec'%prefix][i] = ext['dec']
        row_dict['fitn_%s_glon'%prefix][i] = ext['glon']
        row_dict['fitn_%s_glat'%prefix][i] = ext['glat']
        row_dict['fitn_%s_ra_err'%prefix][i] = ext['ra_err']
        row_dict['fitn_%s_dec_err'%prefix][i] = ext['dec_err']
        row_dict['fitn_%s_glon_err'%prefix][i] = ext['glon_err']
        row_dict['fitn_%s_glat_err'%prefix][i] = ext['glat_err']
        row_dict['fitn_%s_pos_err'%prefix][i] = ext['pos_err']
        row_dict['fitn_%s_pos_r68'%prefix][i] = ext['pos_r68']
        row_dict['fitn_%s_pos_r95'%prefix][i] = ext['pos_r95']
        row_dict['fitn_%s_pos_r99'%prefix][i] = ext['pos_r99']
                
        #row_dict['fitn_%s_loglike'%prefix][i] = ext_roi['roi']['loglike']
        row_dict['fitn_%s_loglike'%prefix][i] = ext['loglike_ext']
        if ext_data_psfhi:
            row_dict['fitn_%s_sys_ts_ext'%prefix][i] = max(0,min(ext_data_psfhi[i]['ts_ext'],
                                                                 ext['ts_ext'],
                                                                 ext_data_psflo[i]['ts_ext']))
            row_dict['fitn_%s_r68_sys_err'%prefix][i] = 0.5*np.abs(ext_data_psfhi[i]['ext'] - ext_data_psflo[i]['ext'])
            row_dict['fitn_%s_r68_sys_ul95'%prefix][i] = max(ext_data_psfhi[i]['ext_ul95'],
                                                            ext['ext_ul95'],
                                                            ext_data_psflo[i]['ext_ul95'])
        
        #src = ext_roi['sources'][src_name]
        src = ext_roi[ext_roi['name'] == src_name][0]

        for x in ['','100','1000','10000']:
            row_dict['fitn_%s_flux%s'%(prefix,x)][i] = src['flux%s'%x]
            row_dict['fitn_%s_flux%s_err'%(prefix,x)][i] = src['flux%s_err'%x]
            row_dict['fitn_%s_eflux%s'%(prefix,x)][i] = src['eflux%s'%x]
            row_dict['fitn_%s_eflux%s_err'%(prefix,x)][i] = src['eflux%s_err'%x]
        row_dict['fitn_%s_ts'%prefix][i] = src['ts']
        row_dict['fitn_%s_npred'%prefix][i] = src['npred']

        index, index_err = extract_index(src)
        row_dict['fitn_%s_index'%prefix][i] = np.abs(index)
        row_dict['fitn_%s_index_err'%prefix][i] = np.abs(index_err)

def extract_halo_data(halo_data, halo_scan_shape):
    
    o = dict(index = [],
             r68 = [],
             ts = [],
             eflux = [],
             eflux_ul95 = [],
             dlnl = [],
             loglike = [],
             dlnl_eflux = [])

    if halo_data.ndim == 0 or len(halo_data) == 0:
        return {}
    
    for hd in halo_data:

        o['index'] += [np.abs(hd['param_values'][hd['param_names']=='Index'][0])]
        o['r68'] += [hd['SpatialWidth']]

        ts = 2.0*(np.max(hd['dloglike_scan']) - hd['dloglike_scan'][0])        
        o['ts'] += [ts]
#        print '%10.3f %10.3f %10.3f %10.3f'%(np.abs(hd['params']['Index'][0]),
#                                             hd['SpatialWidth'], ts, hd['ts'])

        o['eflux'] += [hd['eflux']]
        o['eflux_ul95'] += [hd['eflux_ul95']]
        o['dlnl'] += list(hd['dloglike_scan'])
        o['loglike'] += list(hd['loglike_scan'])
        o['dlnl_eflux'] += list(hd['eflux_scan'])

    for k, v in o.items():
        o[k] = np.array(o[k])

    o['r68'] = o['r68'].reshape(halo_scan_shape)
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

def min_separation(x,y):

    sep = []
    
    for i in range(len(x)):

        for j in range(i+1,len(x)):

            sep += [np.sqrt((x[i]-x[j])**2 + (y[i]-y[j])**2)]

    if not len(sep):
        return np.nan
            
    return np.min(np.array(sep))

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

    nfit = 5
    nscan_pts = 9

    # SED Table
    
    cols_dict_sed = OrderedDict()
    cols_dict_sed['name'] = dict(dtype='S20', format='%s',description='Source Name')
    cols_dict_sed['norm'] = dict(dtype='f8', format='%.3f',shape=(nebins))
    cols_dict_sed['norm_err'] = dict(dtype='f8', format='%.3f',shape=(nebins))
    cols_dict_sed['norm_errp'] = dict(dtype='f8', format='%.3f',shape=(nebins))
    cols_dict_sed['norm_errn'] = dict(dtype='f8', format='%.3f',shape=(nebins))
    cols_dict_sed['norm_ul'] = dict(dtype='f8', format='%.3f',shape=(nebins))
    cols_dict_sed['ref_flux'] = dict(dtype='f8', format='%.3f',shape=(nebins))
    cols_dict_sed['ref_eflux'] = dict(dtype='f8', format='%.3f',shape=(nebins))
    cols_dict_sed['ref_npred'] = dict(dtype='f8', format='%.3f',shape=(nebins))
    cols_dict_sed['ref_dnde'] = dict(dtype='f8', format='%.3f',shape=(nebins), unit='1 / (MeV cm2 s)')
    cols_dict_sed['norm_scan'] = dict(dtype='f8', format='%.3f',shape=(nebins,nscan_pts))
    cols_dict_sed['dloglike_scan'] = dict(dtype='f8', format='%.3f',shape=(nebins,nscan_pts))

    # LNL Table
    
    cols_dict_lnl = OrderedDict()
    cols_dict_lnl['name'] = dict(dtype='S20', format='%s',description='Source Name')

    cols_dict_lnl['fit_ext_gauss_scan_dlnl'] = dict(dtype='f8', format='%.3f',shape=(next_bins))
    cols_dict_lnl['fit_ext_gauss_scan_psfhi_dlnl'] = dict(dtype='f8', format='%.3f',shape=(next_bins))
    cols_dict_lnl['fit_ext_gauss_scan_psflo_dlnl'] = dict(dtype='f8', format='%.3f',shape=(next_bins))
    cols_dict_lnl['fit_ext_disk_scan_dlnl'] = dict(dtype='f8', format='%.3f',shape=(next_bins))
    cols_dict_lnl['fit_ext_disk_scan_psfhi_dlnl'] = dict(dtype='f8', format='%.3f',shape=(next_bins))
    cols_dict_lnl['fit_ext_disk_scan_psflo_dlnl'] = dict(dtype='f8', format='%.3f',shape=(next_bins))
    cols_dict_lnl['fit_halo_scan_dlnl'] = dict(dtype='f8', format='%.3f',
                                               shape=halo_scan_shape + (len(eflux_scan_pts),))
    cols_dict_lnl['fit_halo_sed_scan_dlnl'] = dict(dtype='f8', format='%.3f',
                                                   shape=(halo_scan_shape[0],nebins,len(eflux_scan_pts)))    
    cols_dict_lnl['fit_src_sed_scan_dlnl'] = dict(dtype='f8', format='%.3f',
                                                   shape=(nebins,nscan_pts))
    cols_dict_lnl['fit_src_sed_scan_eflux'] = dict(dtype='f8', format='%.4g',
                                                   shape=(nebins,nscan_pts))

    # CATALOG Table
    
    cols_dict = OrderedDict()
    cols_dict['name'] = dict(dtype='S32', format='%s',description='Source Name')
    cols_dict['codename'] = dict(dtype='S32', format='%s')
    cols_dict['linkname'] = dict(dtype='S32', format='%s')
    cols_dict['assoc'] = dict(dtype='S32', format='%s')
    cols_dict['redshift'] = dict(dtype='f8')
    cols_dict['var_index'] = dict(dtype='f8')
    cols_dict['3lac_radio_flux'] = dict(dtype='f8')
    cols_dict['3lac_radio_nu'] = dict(dtype='f8')
    cols_dict['3lac_radio_nufnu'] = dict(dtype='f8')
    cols_dict['3lac_xray_flux'] = dict(dtype='f8')
    cols_dict['3lac_fx_fr'] = dict(dtype='f8')
    cols_dict['name_3fgl'] = dict(dtype='S32', format='%s')
    cols_dict['assoc_3fgl'] = dict(dtype='S20', format='%s',description='Source Name',shape=(5,))
    cols_dict['assoc_3fgl_sep'] = dict(dtype='f8', format='%s',shape=(5,))
    cols_dict['assoc_3fgl_sep_sigma95'] = dict(dtype='f8', format='%s',shape=(5,))
    cols_dict['assoc_3fgl_num'] = dict(dtype='i8')
    
    cols_dict['name_3fhl'] = dict(dtype='S32', format='%s')
    cols_dict['assoc_3fhl'] = dict(dtype='S20', format='%s',description='Source Name',shape=(5,))
    cols_dict['assoc_3fhl_sep'] = dict(dtype='f8', format='%s',shape=(5,))
    cols_dict['assoc_3fhl_sep_sigma95'] = dict(dtype='f8', format='%s',shape=(5,))
    cols_dict['assoc_3fhl_num'] = dict(dtype='i8')

    cols_dict['class'] = dict(dtype='S32', format='%s')
    cols_dict['class_optical'] = dict(dtype='S32', format='%s')
    cols_dict['class_sed'] = dict(dtype='S32', format='%s')
    cols_dict['nupeak'] = dict(dtype='f8')
    cols_dict['sep_3fgl'] = dict(dtype='f8', format='%.2f')
    cols_dict['sep_3fgl_sigma95'] = dict(dtype='f8', format='%.2f')
    cols_dict['sep_3fhl'] = dict(dtype='f8', format='%.2f')
    cols_dict['sep_3fhl_sigma95'] = dict(dtype='f8', format='%.2f')

    for k in src_colnames:
        cols_dict[k] = dict(dtype='f8', format='%.3f')
        cols_dict['fitn_%s'%k] = dict(dtype='f8', format='%.2f',shape=(nfit,))
        cols_dict['fit_%s'%k] = dict(dtype='f8', format='%.2f')
    
    cols_dict['spatial_model'] = dict(dtype='S32', format='%s')
    cols_dict['fit_ext_model'] = dict(dtype='S32', format='%s')
    
    for t in ['dnde','flux','eflux']:        
        for x in ['','100','1000','10000']:        
            cols_dict['%s%s'%(t,x)] = dict(dtype='f8', format='%.4g')
            cols_dict['%s%s_err'%(t,x)] = dict(dtype='f8', format='%.4g')

    cols_dict['spectrum_type'] = dict(dtype='S32', format='%s')

    # Halo Fits
    for k in halo_colnames:    
        cols_dict['fitn_halo_%s'%k] = dict(dtype='f8', format='%.2f',shape=(nfit,))
        cols_dict['fit_halo_%s'%k] = dict(dtype='f8', format='%.2f')

    # Ext Fits
    for k in ext_colnames:

        if 'flux' in k:
            fmt = '%.5g'
        else:
            fmt = '%.2f'

        cols_dict['fit_ext_%s'%k] = dict(dtype='f8', format=fmt)
        cols_dict['fitn_ext_gauss_%s'%k] = dict(dtype='f8', format=fmt,shape=(nfit,))    
        cols_dict['fit_ext_gauss_%s'%k] = dict(dtype='f8', format=fmt)
        cols_dict['fitn_ext_disk_%s'%k] = dict(dtype='f8', format=fmt,shape=(nfit,))    
        cols_dict['fit_ext_disk_%s'%k] = dict(dtype='f8', format=fmt)
        
        
    # Delta LogLikes
    cols_dict['fitn_dlike'] = dict(dtype='f8', format='%.2f',shape=(nfit,))
    cols_dict['fitn_dlike1'] = dict(dtype='f8', format='%.2f',shape=(nfit,))

    for k in ['ext_gauss','ext_disk','halo']:    
        cols_dict['fitn_dlike_ps_%s'%k] = dict(dtype='f8', format='%.2f',shape=(nfit,))
        cols_dict['fitn_daic_ps_%s'%k] = dict(dtype='f8', format='%.2f',shape=(nfit,))
        cols_dict['fitn_aic_%s'%k] = dict(dtype='f8', format='%.2f',shape=(nfit,))
        cols_dict['fitn_dlike1_%s'%k] = dict(dtype='f8', format='%.2f',shape=(nfit,))
        cols_dict['fitn_dlike_%s'%k] = dict(dtype='f8', format='%.2f',shape=(nfit,))
        
    cols_dict['fitn_aic_ps'] = dict(dtype='f8', format='%.2f',shape=(nfit,))

    # Best-fit model
    for k in ['ext_gauss','ext_disk','ext','halo']:    
        cols_dict['fit_dlike_%s'%k] = dict(dtype='f8', format='%.2f')
        cols_dict['fit_dlike1_%s'%k] = dict(dtype='f8', format='%.2f')
        cols_dict['fit_dlike_ps_%s'%k] = dict(dtype='f8', format='%.2f')
        cols_dict['fit_daic_ps_%s'%k] = dict(dtype='f8', format='%.2f')
        cols_dict['fit_aic_%s'%k] = dict(dtype='f8', format='%.2f')
    
    cols_dict['fit_halo_scan_ts'] = dict(dtype='f8', format='%.2f',shape=halo_scan_shape)
    cols_dict['fit_halo_scan_eflux_ul95'] = dict(dtype='f8', format='%.4g',shape=halo_scan_shape)

    cols_dict['fitn_nparam_ps'] = dict(dtype='i8',shape=(nfit,))
    cols_dict['fitn_nparam_ext_gauss'] = dict(dtype='i8',shape=(nfit,))
    cols_dict['fitn_nparam_ext_disk'] = dict(dtype='i8',shape=(nfit,))
    cols_dict['fitn_nparam_halo'] = dict(dtype='i8',shape=(nfit,))
    cols_dict['fit_idx'] = dict(dtype='i8')
    cols_dict['fit_idx_ext_gauss'] = dict(dtype='i8')
    cols_dict['fit_idx_ext_disk'] = dict(dtype='i8')
    cols_dict['fit_idx_ext'] = dict(dtype='i8')
    cols_dict['fit_idx_halo'] = dict(dtype='i8')
    cols_dict['fit_mean_sep'] = dict(dtype='f8', format='%.3f')
    cols_dict['fit_min_sep'] = dict(dtype='f8', format='%.3f')

    tab = Table([Column(name=k, **v) for k,v in cols_dict.items()])
    tab_lnl = Table([Column(name=k, **v) for k,v in cols_dict_lnl.items()])
    tab_sed = Table([Column(name=k, **v) for k,v in cols_dict_sed.items()])

    return tab, tab_lnl, tab_sed


def load_table(filename):

    if not os.path.isfile(filename):
        return None
    else:
        return Table.read(filename)


def load_npy_files(dirname, suffix, wildcard=False):

    o = []
    for i in range(5):
        infile = os.path.join(dirname,'fit%i%s.npy'%(i,suffix))
        if wildcard:
            files = sorted(glob.glob(infile))
            if len(files) == 0:
                break
            
            o += [ [np.load(f).flat[0] for f in files] ]
        else:
        
            if not os.path.isfile(infile):
                break

            data = np.load(infile)
            if data.ndim == 0:
                o += [data.flat[0]]
            else:
                o += [data]
    return o

def table_to_dict(tab):

    o = {}
    for c in tab.columns:
        o[c] = tab[c][0]

    return o
        
def load_tables(dirname, suffix, todict=False, hdu=None):

    o = []
    for i in range(5):
        infile = os.path.join(dirname,'fit%i%s.fits'%(i,suffix))
        if not os.path.isfile(infile):
            continue

        tab = Table.read(infile,hdu)
        if todict:
            tab = table_to_dict(tab)
        o += [tab]
        
    return o
    

def count_free_params(tab):

    if len(tab) == 0:
        return 0
    
    nfree = 0    
    for row in tab:

        # for spatial parameters
        nfree += 2
        
        # for spectral parameters
        if row['SpectrumType'] == 'PowerLaw':
            nfree += 2
        else:
            nfree += 3

    return nfree


def find_best_model(row, model_name, par_name, daic_threshold=0.0):

    idx = 0
    
    for i in range(5):

        if not np.isfinite(row['fitn_dlike1'][i]):
            break
        if not np.isfinite(row['fitn_daic_ps_%s'%model_name][i]):
            idx = i
            break
            
        delta_ts = row['fitn_%s_%s'%(model_name,par_name)][i+1] - row['fitn_%s_%s'%(model_name,par_name)][i]
        #print(i,row['fitn_daic_ps_%s'%model_name][i],delta_ts)        
        
        if row['fitn_daic_ps_%s'%model_name][i] < daic_threshold and delta_ts < -1.0:
            idx = i
            break

    return idx

def aggregate(dirs,output,suffix=''):

    dirs = [d for argdir in dirs for d in collect_dirs(argdir)]

    nfit = 5
    nscan_pts = 9
    nebins = 20
    next_bins = 30
    halo_scan_shape = (15,9)
    eflux_scan_pts = np.logspace(-10,-4,61)
    
    row_dict = {}
    row_dict_lnl = {}
    row_dict_sed = {}
    tab = None
    tab_lnl = None
    tab_sed = None
    halo_name = 'halo_RadialGaussian'
    
    #roi = ROIModel({'catalogs' : ['/u/gl/mdwood/fermi/catalogs/gll_psc_v16_ext.fit','lmc_psc_v0.fit'],
    #                'extdir' : '/u/gl/mdwood/fermi/catalogs/Extended_archive_v16'})

    cat_3lac = Table.read('/u/gl/mdwood/fermi/catalogs/3lac_fx_fr_official.fits')
    if 'Fermi name' in cat_3lac.columns:
        cat_3lac['Fermi name'].name='Source_Name'
    cat_3fgl = Catalog3FGL('/u/gl/mdwood/fermi/catalogs/gll_psc_v16_ext.fit').table
    cat_3fgl = join(cat_3fgl, cat_3lac, join_type='outer', keys='Source_Name')
    cat_3fgl['Redshift'].fill_value = np.nan
    cat_3fgl['Radio flux(mJy)'].fill_value = np.nan
    cat_3fgl['X Flux(erg cm-2 s-1)'].fill_value = np.nan
    cat_3fgl['radio_flux'].fill_value = np.nan
    cat_3fgl['radio_freq'].fill_value = np.nan
    cat_3fgl['FX_FR'].fill_value = np.nan
    cat_3fgl['X Flux(erg cm-2 s-1)'].fill_value = np.nan
    cat_3fgl['Optical Class'].fill_value = ''
    cat_3fgl['SED Class'].fill_value = ''
    cat_3fgl['Log Sync.Z corr.(Hz)'][cat_3fgl['Log Sync.Z corr.(Hz)'] == np.inf] = np.nan
    cat_3fgl['radio_freq'][cat_3fgl['radio_freq'] == 0.0] = np.nan
    cat_3fgl['FX_FR'][cat_3fgl['FX_FR'] == 0.0] = np.nan
    cat_3fgl['SED Class'][cat_3fgl['SED Class'] == '-'] = ''
    cat_3fgl = cat_3fgl.filled()    
    
    cat_3fhl = Catalog3FGL('/u/gl/mdwood/fermi/catalogs/gll_psch_v11.fit').table
    skydir_3fgl = SkyCoord(cat_3fgl['RAJ2000'],cat_3fgl['DEJ2000'],unit='deg')
    skydir_3fhl = SkyCoord(cat_3fhl['RAJ2000'],cat_3fhl['DEJ2000'],unit='deg')

    model_3fgl = np.sqrt(np.array(cat_3fgl['Model_SemiMajor'])*np.array(cat_3fgl['Model_SemiMinor']))
    sigma95_3fgl = np.sqrt(np.array(cat_3fgl['Conf_95_SemiMajor'])*np.array(cat_3fgl['Conf_95_SemiMinor']))
    sigma95_3fgl[~np.isfinite(sigma95_3fgl)] = 0.1
    sigma95_3fgl[np.isfinite(model_3fgl)] += model_3fgl[np.isfinite(model_3fgl)]

    model_3fhl = np.sqrt(np.array(cat_3fhl['Model_SemiMajor'])*np.array(cat_3fhl['Model_SemiMinor']))
    sigma95_3fhl = np.sqrt(np.array(cat_3fhl['Conf_95_SemiMajor'])*np.array(cat_3fhl['Conf_95_SemiMinor']))
    sigma95_3fhl[~np.isfinite(sigma95_3fhl)] = 0.1
    sigma95_3fhl[np.isfinite(model_3fhl)] += model_3fhl[np.isfinite(model_3fhl)]
    
    for d in sorted(dirs):

        logfile0 = os.path.join(d,'run-region-analysis.log')
        logfile1 = os.path.join(d,'run-halo-analysis.log')

        #if not os.path.isfile(logfile0):
        #    continue

        if not os.path.isfile(logfile0) or not check_log(logfile0)=='Successful':
            print('skipping because region failed', d)
            continue

        if not os.path.isfile(logfile1) or not check_log(logfile1)=='Successful':
            print('skipping because halo failed', d)
            continue
        
        #if not os.path.isfile(logfile1) or check_log(logfile1)=='Exited':
        #    continue
        
        print d
        
        if not os.path.isfile(os.path.join(d,'new_source_data.npy')):
            print 'skipping'
            continue

        new_src_data = np.load(os.path.join(d,'new_source_data.npy'))
        fit_data = load_tables(d,'%s_roi'%suffix)
        roi_data = load_tables(d,'%s_roi'%suffix, hdu='ROI')

        tab_halo_data = load_tables(d,'%s_%s_data'%(suffix,halo_name))
        tab_new_src_data = [Table()] + load_tables(d,'%s_new_source_data'%(suffix))
        
        halo_data = load_npy_files(d,'%s_%s_data'%(suffix,halo_name))
        ext_gauss_data = load_tables(d,'%s_ext_gauss_ext'%suffix, True)
        ext_gauss_data_psfhi = load_tables(d,'%s_ext_gauss_ext_psfhi'%suffix, True)
        ext_gauss_data_psflo = load_tables(d,'%s_ext_gauss_ext_psflo'%suffix, True)
        ext_gauss_roi = load_tables(d,'%s_ext_gauss_roi'%suffix)
        ext_disk_data = load_tables(d,'%s_ext_disk_ext'%suffix, True)
        ext_disk_data_psfhi = load_tables(d,'%s_ext_disk_ext_psfhi'%suffix, True)
        ext_disk_data_psflo = load_tables(d,'%s_ext_disk_ext_psflo'%suffix, True)
        ext_disk_roi = load_tables(d,'%s_ext_disk_roi'%suffix)
        sed_data = load_npy_files(d,'%s_sed'%suffix)
        halo_data_sed = load_npy_files(d,'_cov05_*_halo_radialgaussian_sed',wildcard=True)
        
        if len(fit_data) == 0:
            print 'skipping'
            continue

        config = yaml.load(fit_data[0].meta['CONFIG'])        
        src_name = config['selection']['target']
        nebins = len(roi_data[0][0]['model_counts'])
        
        src_data = []
        new_srcs = []        
        for i, fd in enumerate(fit_data):            
            src_data += [fd[fd['name'] == src_name][0]]
        src = src_data[0]
            
        if tab is None:
            tab, tab_lnl, tab_sed = create_tables(nebins,
                                                  len(ext_gauss_data[0]['width']),
                                                  halo_scan_shape,
                                                  eflux_scan_pts)
            
        row_dict = {}
        for c in tab.columns:
            
            if tab.columns[c].dtype.kind in ['U','S']:

                if tab.columns[c].ndim == 2:
                    row_dict[c] = np.zeros(tab.columns[c].shape[1:],
                                           dtype=tab.columns[c].dtype)
                else:
                    row_dict[c] = ''
            elif tab.columns[c].dtype.kind in ['i']:
                row_dict[c] = np.zeros(tab.columns[c].shape[1:])
            else:
                row_dict[c] = np.ones(tab.columns[c].shape[1:])*np.nan
        
        halo_seds = []
        for hd in halo_data_sed:
            halo_seds += [extract_halo_sed(hd,halo_scan_shape)]

        ext_r68 = []
        for i, (ext,ext_roi) in enumerate(zip(ext_gauss_data,ext_gauss_roi)):
            ext_r68 += [ext['width']]

        extract_ext_data(row_dict, 'ext_gauss', src_name,
                         ext_gauss_data, ext_gauss_roi, ext_gauss_data_psfhi, ext_gauss_data_psflo)

        extract_ext_data(row_dict, 'ext_disk', src_name,
                         ext_disk_data, ext_disk_roi, ext_disk_data_psfhi, ext_disk_data_psflo)
                    
        halo_scan_r68 = np.array(tab_halo_data[0]['halo_width']).reshape(halo_scan_shape)
        halo_scan_index = np.array(tab_halo_data[0]['halo_index']).reshape(halo_scan_shape)
            
        for i in range(len(tab_halo_data)):

            halo_scan_ts = np.array(tab_halo_data[i]['ts']).reshape(halo_scan_shape)
            halo_scan_loglike = np.max(np.array(tab_halo_data[i]['loglike_scan']),axis=1).reshape(halo_scan_shape)
            halo_scan_eflux = np.array(tab_halo_data[i]['eflux']).reshape(halo_scan_shape)
            halo_scan_eflux_err = np.array(tab_halo_data[i]['eflux_err']).reshape(halo_scan_shape)
            halo_scan_eflux_ul95 = np.array(tab_halo_data[i]['eflux_ul95']).reshape(halo_scan_shape)
            halo_scan_flux = np.array(tab_halo_data[i]['flux']).reshape(halo_scan_shape)
            halo_scan_flux_err = np.array(tab_halo_data[i]['flux_err']).reshape(halo_scan_shape)
            halo_scan_flux_ul95 = np.array(tab_halo_data[i]['flux_ul95']).reshape(halo_scan_shape)
            halo_scan_idx = np.argmax(halo_scan_ts)
            halo_scan_max_ts = halo_scan_ts.flat[halo_scan_idx]
            halo_scan_max_loglike = halo_scan_loglike.flat[halo_scan_idx]
            halo_scan_max_eflux = halo_scan_eflux.flat[halo_scan_idx]
            halo_scan_max_eflux_ul95 = halo_scan_eflux_ul95.flat[halo_scan_idx]
            halo_scan_max_flux = halo_scan_flux.flat[halo_scan_idx]
            halo_scan_max_flux_ul95 = halo_scan_flux_ul95.flat[halo_scan_idx]
            halo_scan_max_r68 = tab_halo_data[0]['halo_width'].flat[halo_scan_idx]
            halo_scan_max_index = tab_halo_data[0]['halo_index'].flat[halo_scan_idx] 
            
            if halo_scan_max_ts > 2.0:
                
                ix, iy = np.unravel_index(np.argmax(halo_scan_ts),halo_scan_ts.shape)
                o = fit_parabola(halo_scan_ts,ix,iy,dpix=2)
                pix = np.array([[o['x0']],[o['y0']]])
                halo_fit_ts = map_coordinates(halo_scan_ts, pix, mode='nearest')
                halo_fit_loglike = map_coordinates(halo_scan_loglike, pix, mode='nearest')
                halo_fit_eflux = map_coordinates(halo_scan_eflux, pix, mode='nearest')
                halo_fit_eflux_err = map_coordinates(halo_scan_eflux_err,pix, mode='nearest')
                halo_fit_eflux_ul95 = map_coordinates(halo_scan_eflux_ul95,pix, mode='nearest')
                halo_fit_flux = map_coordinates(halo_scan_flux,pix, mode='nearest')
                halo_fit_flux_err = map_coordinates(halo_scan_flux_err,pix, mode='nearest')
                halo_fit_flux_ul95 = map_coordinates(halo_scan_flux_ul95,pix, mode='nearest')
                halo_fit_r68 = map_coordinates(halo_scan_r68,pix, mode='nearest')
                halo_fit_index = map_coordinates(halo_scan_index, pix, mode='nearest')

                # Interpolation
                x = np.linspace(0,halo_scan_shape[0]-1,halo_scan_shape[0])
                y = np.linspace(0,halo_scan_shape[1]-1,halo_scan_shape[1])
                interp_ts_x = map_coordinates(halo_scan_ts,np.array([x,o['y0']*np.ones_like(x)]))
                interp_ts_y = map_coordinates(halo_scan_ts,np.array([o['x0']*np.ones_like(y),y]))
                halo_r68_lims = get_parameter_limits(halo_scan_r68[:,iy], 0.5*interp_ts_x)
                halo_index_lims = get_parameter_limits(halo_scan_index[ix,:], 0.5*interp_ts_y)
                halo_fit_r68_err = halo_r68_lims['err']
                halo_fit_index_err = halo_index_lims['err']
                
                row_dict['fitn_halo_ts'][i] = halo_fit_ts
                row_dict['fitn_halo_loglike'][i] = halo_fit_loglike
                row_dict['fitn_halo_eflux'][i] = halo_fit_eflux
                row_dict['fitn_halo_eflux_err'][i] = halo_fit_eflux_err
                row_dict['fitn_halo_eflux_ul95'][i] = halo_fit_eflux_ul95
                row_dict['fitn_halo_flux'][i] = halo_fit_flux
                row_dict['fitn_halo_flux_err'][i] = halo_fit_flux_err
                row_dict['fitn_halo_flux_ul95'][i] = halo_fit_flux_ul95
                row_dict['fitn_halo_r68'][i] = halo_fit_r68
                row_dict['fitn_halo_r68_err'][i] = halo_fit_r68_err
                row_dict['fitn_halo_index'][i] = halo_fit_index
                row_dict['fitn_halo_index_err'][i] = halo_fit_index_err
            else:
                
                row_dict['fitn_halo_ts'][i] = halo_scan_max_ts
                row_dict['fitn_halo_loglike'][i] = halo_scan_max_loglike
                row_dict['fitn_halo_eflux'][i] = halo_scan_max_eflux
                row_dict['fitn_halo_eflux_ul95'][i] = halo_scan_max_eflux_ul95
                row_dict['fitn_halo_flux'][i] = halo_scan_max_flux
                row_dict['fitn_halo_flux_ul95'][i] = halo_scan_max_flux_ul95
                row_dict['fitn_halo_r68'][i] = halo_scan_max_r68
                row_dict['fitn_halo_index'][i] = halo_scan_max_index            

        for i in range(len(src_data)):
            row_dict['fitn_ts'][i] = src_data[i]['ts']
            row_dict['fitn_npred'][i] = src_data[i]['npred']
            row_dict['fitn_offset'][i] = src_data[i]['offset']
            row_dict['fitn_ra'][i] = src_data[i]['ra']
            row_dict['fitn_dec'][i] = src_data[i]['dec']
            row_dict['fitn_glon'][i] = src_data[i]['glon']
            row_dict['fitn_glat'][i] = src_data[i]['glat']
            row_dict['fitn_pos_err'][i] = src_data[i]['pos_err']
            row_dict['fitn_pos_r68'][i] = src_data[i]['pos_r68']
            row_dict['fitn_pos_r95'][i] = src_data[i]['pos_r95']
            row_dict['fitn_pos_r99'][i] = src_data[i]['pos_r99']
            index, index_err = extract_index(src_data[i], 1E3)
            row_dict['fitn_index'][i] = np.abs(index)
            row_dict['fitn_index_err'][i] = np.abs(index_err)

        idx = len(src_data)-1
        for colname in src_colnames:
            row_dict['fit_%s'%colname] = row_dict['fitn_%s'%colname][idx]
        
        # AIC = 2k - 2lnL        

        for i in range(len(fit_data)):

            loglike = roi_data[i][0]['loglike']
            loglike1 = roi_data[0][0]['loglike']            
            loglike_ext_gauss = row_dict['fitn_ext_gauss_loglike'][i]
            loglike_ext_disk = row_dict['fitn_ext_disk_loglike'][i]
            loglike_halo = row_dict['fitn_halo_loglike'][i] #loglike+row_dict['fitn_halo_ts'][i]/2.
            
            row_dict['fitn_dlike1'][i] = loglike - loglike1
            row_dict['fitn_dlike1_ext_gauss'][i] = loglike_ext_gauss - loglike1
            row_dict['fitn_dlike1_ext_disk'][i] = loglike_ext_disk - loglike1
            row_dict['fitn_dlike1_halo'][i] = loglike_halo - loglike1
            row_dict['fitn_nparam_ps'][i] = count_free_params(tab_new_src_data[i])
            row_dict['fitn_nparam_halo'][i] = row_dict['fitn_nparam_ps'][i] + 3
            row_dict['fitn_nparam_ext_gauss'][i] = row_dict['fitn_nparam_ps'][i] + 1
            row_dict['fitn_nparam_ext_disk'][i] = row_dict['fitn_nparam_ps'][i] + 1
            row_dict['fitn_aic_ps'][i] = -2.0*(row_dict['fitn_nparam_ps'][i] - row_dict['fitn_dlike1'][i])
            row_dict['fitn_aic_halo'][i] = -2.0*(row_dict['fitn_nparam_halo'][i] - row_dict['fitn_dlike1_halo'][i])
            row_dict['fitn_aic_ext_gauss'][i] = -2.0*(row_dict['fitn_nparam_ext_gauss'][i] - row_dict['fitn_dlike1_ext_gauss'][i])
            row_dict['fitn_aic_ext_disk'][i] = -2.0*(row_dict['fitn_nparam_ext_disk'][i] - row_dict['fitn_dlike1_ext_disk'][i])
            
        # PS+1 : x0, y0, norm0, index0, x1, y1, norm1, index1
        # Ext: R68, x0, y0, norm0, index0
        # Halo : R68, x0, y0, norm0, index0, norm1, index1

        row_dict['fitn_dlike'][1:] = row_dict['fitn_dlike1'][1:] - row_dict['fitn_dlike1'][:-1]
        row_dict['fitn_dlike_ext_gauss'][1:] = row_dict['fitn_dlike1_ext_gauss'][1:] - row_dict['fitn_dlike1_ext_gauss'][:-1]
        row_dict['fitn_dlike_ext_disk'][1:] = row_dict['fitn_dlike1_ext_disk'][1:] - row_dict['fitn_dlike1_ext_disk'][:-1]
        row_dict['fitn_dlike_halo'][1:] = row_dict['fitn_dlike1_halo'][1:] - row_dict['fitn_dlike1_halo'][:-1]
                
        # L_(N+1) - L_(ext,halo)
        row_dict['fitn_dlike_ps_ext_gauss'][:-1] = row_dict['fitn_dlike1'][1:] - row_dict['fitn_dlike1_ext_gauss'][:-1]
        row_dict['fitn_dlike_ps_ext_disk'][:-1] = row_dict['fitn_dlike1'][1:] - row_dict['fitn_dlike1_ext_disk'][:-1]
        row_dict['fitn_dlike_ps_halo'][:-1] = row_dict['fitn_dlike1'][1:] - row_dict['fitn_dlike1_halo'][:-1]
        row_dict['fitn_daic_ps_ext_gauss'][:-1] = row_dict['fitn_aic_ps'][1:] - row_dict['fitn_aic_ext_gauss'][:-1]
        row_dict['fitn_daic_ps_ext_disk'][:-1] = row_dict['fitn_aic_ps'][1:] - row_dict['fitn_aic_ext_disk'][:-1]
        row_dict['fitn_daic_ps_halo'][:-1] = row_dict['fitn_aic_ps'][1:] - row_dict['fitn_aic_halo'][:-1]

        # -(n-L)
        
        # -nparam_ps[i+1] + nparam_ext[i]
        

        
        print '-'*80
        print('%4s %10s %10s %10s %10s %10s %10s %4s'%('iter','dlike1','dlike1_ext','dlike','TS/2','dAIC','dlnL','N'))
        for i in range(len(fit_data)):
            print '%4i %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %4i %8.2f'%(i,
                                                             row_dict['fitn_dlike1'][i],
                                                             row_dict['fitn_dlike1_ext_gauss'][i],
                                                             row_dict['fitn_dlike'][i],
                                                             0.5*row_dict['fitn_ext_gauss_ts_ext'][i],
                                                             row_dict['fitn_daic_ps_ext_gauss'][i],
                                                             row_dict['fitn_dlike_ps_ext_gauss'][i],
                                                                   row_dict['fitn_nparam_ext_gauss'][i],
                                                                             row_dict['fitn_ext_gauss_r68'][i])

        print '+'*80
        for i in range(len(fit_data)):
            print '%4i %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %4i'%(i,
                                                             row_dict['fitn_dlike1'][i],
                                                             row_dict['fitn_dlike1_ext_disk'][i],
                                                             row_dict['fitn_dlike'][i],
                                                             0.5*row_dict['fitn_ext_disk_ts_ext'][i],
                                                             row_dict['fitn_daic_ps_ext_disk'][i],
                                                             row_dict['fitn_dlike_ps_ext_disk'][i],
                                                                       row_dict['fitn_nparam_ext_disk'][i])
            
        print '+'*80
        for i in range(len(fit_data)):            
            print '%4i %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f'%(i,
                                                             row_dict['fitn_dlike1'][i],
                                                             row_dict['fitn_dlike1_halo'][i],
                                                             row_dict['fitn_dlike'][i],
                                                             0.5*row_dict['fitn_halo_ts'][i],
                                                             row_dict['fitn_daic_ps_halo'][i],
                                                             row_dict['fitn_dlike_ps_halo'][i])

        fit_idx = len(fit_data)-1
        fit_idx_halo = find_best_model(row_dict,'halo','ts')
        fit_idx_ext_gauss = find_best_model(row_dict,'ext_gauss','ts_ext')
        fit_idx_ext_disk = find_best_model(row_dict,'ext_disk','ts_ext')
        
        # ---------------------------------------------------------
        # Halo Results
        for k in halo_colnames:
            row_dict['fit_halo_%s'%k] = row_dict['fitn_halo_%s'%k][fit_idx_halo]
        
        # ---------------------------------------------------------
        # Extension Results
        for k in ext_colnames:
            row_dict['fit_ext_gauss_%s'%k] = row_dict['fitn_ext_gauss_%s'%k][fit_idx_ext_gauss]
            row_dict['fit_ext_disk_%s'%k] = row_dict['fitn_ext_disk_%s'%k][fit_idx_ext_disk]

        for k in ['dlike_ps_halo','dlike1_halo','dlike_halo','daic_ps_halo']:
            row_dict['fit_%s'%k] = row_dict['fitn_%s'%k][fit_idx_halo]

        for k in ['dlike_ps','dlike1','dlike','daic_ps', 'aic']:
            row_dict['fit_%s_ext_gauss'%k] = row_dict['fitn_%s_ext_gauss'%k][fit_idx_ext_gauss]
            row_dict['fit_%s_ext_disk'%k] = row_dict['fitn_%s_ext_disk'%k][fit_idx_ext_disk]


        if row_dict['fit_ext_gauss_ts_ext'] > 16.:

            if fit_idx_ext_disk != fit_idx_ext_gauss:
                if fit_idx_ext_disk < fit_idx_ext_gauss:
                    fit_ext_model = 'ext_disk'
                    row_dict['fit_ext_model'] = 'disk'
                else:
                    fit_ext_model = 'ext_gauss'
                    row_dict['fit_ext_model'] = 'gauss'                
            elif row_dict['fit_aic_ext_disk'] < row_dict['fit_aic_ext_gauss']:
                fit_ext_model = 'ext_gauss'
                row_dict['fit_ext_model'] = 'gauss'
            else:
                fit_ext_model = 'ext_disk'
                row_dict['fit_ext_model'] = 'disk'            
        else:
            fit_ext_model = 'ext_gauss'
            row_dict['fit_ext_model'] = 'gauss'
                        
        #if ((row_dict['fit_ext_gauss_ts_ext'] < 16.) or
        #    (fit_idx_ext_disk != fit_idx_ext_gauss and fit_idx_ext_disk > fit_idx_ext_gauss) or
        #    (row_dict['fit_aic_ext_disk'] < row_dict['fit_aic_ext_gauss'])):
        #    fit_ext_model = 'ext_gauss'
        #    row_dict['fit_ext_model'] = 'gauss'
        #else:
        #    fit_ext_model = 'ext_disk'
        #    row_dict['fit_ext_model'] = 'disk'
            
            
        #if max(row_dict['fit_ext_disk_ts_ext'],row_dict['fit_ext_disk_ts_ext']) < 16.:
        #    fit_ext_model = 'ext_gauss'
        #    row_dict['fit_ext_model'] = 'gauss'
        #elif fit_idx_ext_disk != fit_idx_ext_gauss:
        #    if fit_idx_ext_disk < fit_idx_ext_gauss:
        #        fit_ext_model = 'ext_disk'
        #        row_dict['fit_ext_model'] = 'disk'
        #    else:
        #        fit_ext_model = 'ext_gauss'
        #        row_dict['fit_ext_model'] = 'gauss'
        #else:
        #    if row_dict['fit_aic_ext_disk'] > row_dict['fit_aic_ext_gauss']:
        #        fit_ext_model = 'ext_disk'
        #        row_dict['fit_ext_model'] = 'disk'
        #    else:
        #        fit_ext_model = 'ext_gauss'
        #        row_dict['fit_ext_model'] = 'gauss'

        for k in ext_colnames:
            row_dict['fit_ext_%s'%k] = row_dict['fit_%s_%s'%(fit_ext_model, k)]

        for k in ['dlike_ps','dlike1','dlike','daic_ps', 'aic']:
            row_dict['fit_%s_ext'%k] = row_dict['fit_%s_%s'%(k,fit_ext_model)]
            
        print('best halo', fit_idx_halo)
        print('best ext gauss', fit_idx_ext_gauss)
        print('best ext disk', fit_idx_ext_disk)
        print('best ext model', fit_ext_model)
        print('AIC gauss',row_dict['fit_aic_ext_gauss'])
        print('AIC disk',row_dict['fit_aic_ext_disk'])

        
        print(row_dict['assoc'])
        
        x = [src_data[-1]['offset_glon']]
        y = [src_data[-1]['offset_glat']]

        for j in range(len(new_src_data)):
            x += [new_src_data[j]['offset_glon']]
            y += [new_src_data[j]['offset_glat']]

        xc, yc = centroid(x,y)
        fit_mean_sep = mean_separation(x,y,xc,yc)
        fit_min_sep = min_separation(x,y)
        if fit_idx == 0:
            fit_mean_sep = np.nan
            fit_min_sep = np.nan

        #src_sed_tab = Table.read('fit%i_%s_sed.fits'%(fit_idx,codename))
            
        row_dict['fit_halo_scan_ts'] = np.nan*np.ones(halo_scan_shape)
        row_dict['fit_halo_scan_eflux_ul95'] = np.nan*np.ones(halo_scan_shape)
        row_dict['fit_halo_scan_ts'] = np.array(tab_halo_data[fit_idx_halo]['ts']).reshape(halo_scan_shape)
        row_dict['fit_halo_scan_eflux_ul95'] = np.array(tab_halo_data[fit_idx_halo]['eflux_ul95']).reshape(halo_scan_shape)

        fitn_ra = np.array([np.nan]*nfit)
        fitn_dec = np.array([np.nan]*nfit)
        fitn_ra[:len(src_data)] = np.array([sd['ra'] for sd in src_data])
        fitn_dec[:len(src_data)] = np.array([sd['dec'] for sd in src_data])

        #for k, v in src_data[0].items():
        #    if k in tab.columns:
        #        row_dict[k] = v

        if row_dict['fit_%s_ts_ext'%fit_ext_model] < 16.:
            c = SkyCoord(fitn_ra[fit_idx],
                         fitn_dec[fit_idx],unit='deg')

            for k in src_colnames:
                row_dict[k] = row_dict['fit_%s'%k]
            
            for x in ['','100','1000','10000']:
                row_dict['eflux%s'%x] = src_data[fit_idx]['eflux%s'%x]
                row_dict['eflux%s_err'%x] = src_data[fit_idx]['eflux%s_err'%x]
                row_dict['flux%s'%x] = src_data[fit_idx]['flux%s'%x]
                row_dict['flux%s_err'%x] = src_data[fit_idx]['flux%s_err'%x]
                
            row_dict['spatial_model'] = src_data[0]['SpatialModel']
            sigma95 = src_data[fit_idx]['pos_r95']
        else:
            c = SkyCoord(row_dict['fit_%s_ra'%fit_ext_model],
                         row_dict['fit_%s_dec'%fit_ext_model],
                         unit='deg')

            for k in src_colnames:
                row_dict[k] = row_dict['fit_%s_%s'%(fit_ext_model,k)]
                        
            for x in ['','100','1000','10000']:
                row_dict['eflux%s'%x] = row_dict['fit_%s_eflux%s'%(fit_ext_model,x)]
                row_dict['eflux%s_err'%x] = row_dict['fit_%s_eflux%s_err'%(fit_ext_model,x)]
                row_dict['flux%s'%x] = row_dict['fit_%s_flux%s'%(fit_ext_model,x)]
                row_dict['flux%s_err'%x] = row_dict['fit_%s_flux%s_err'%(fit_ext_model,x)]
            
            row_dict['spatial_model'] = 'RadialGaussian' if fit_ext_model == 'ext_gauss' else 'RadialDisk'
            sigma95 = row_dict['fit_%s_r68'%fit_ext_model]

        match_catalog(row_dict, cat_3fgl, sigma95, sigma95_3fgl, c, skydir_3fgl, '3fgl')
        match_catalog(row_dict, cat_3fhl, sigma95, sigma95_3fhl, c, skydir_3fhl, '3fhl')
        
        row_dict['name'] = create_source_name(c,prefix='FHES')
        codename = os.path.basename(d)
        linkname = '{%s}'%codename.replace('+','p').replace('.','_')
        
        row_dict['codename'] = codename
        row_dict['linkname'] = linkname

        if row_dict['assoc_3fgl_num'] > 0:
            m = (cat_3fgl['Source_Name'] == row_dict['name_3fgl'])
            row_3fgl = cat_3fgl[m][0]
            row_dict['assoc'] = row_3fgl['ASSOC1']
            row_dict['class'] = row_3fgl['CLASS1']
            row_dict['class_optical'] = row_3fgl['Optical Class']
            row_dict['class_sed'] = row_3fgl['SED Class']
            row_dict['var_index'] = row_3fgl['Variability_Index']
            if row_3fgl['Redshift'] > 0.0:
                row_dict['redshift'] = row_3fgl['Redshift']
            if row_3fgl['Log Sync.Z corr.(Hz)'] > 0.0 and row_3fgl['Log Sync.Z corr.(Hz)'] < 100:
                row_dict['nupeak'] = 10**row_3fgl['Log Sync.Z corr.(Hz)']
            if row_3fgl['Radio flux(mJy)'] > 0.0:
                row_dict['3lac_radio_flux'] = row_3fgl['Radio flux(mJy)']
            if row_3fgl['X Flux(erg cm-2 s-1)'] > 0.0:
                row_dict['3lac_xray_flux'] = row_3fgl['X Flux(erg cm-2 s-1)']

            row_dict['3lac_radio_nu'] = row_3fgl['radio_freq']
            row_dict['3lac_radio_nufnu'] = row_3fgl['radio_flux']
            row_dict['3lac_fx_fr'] = row_3fgl['FX_FR']
                
        elif row_dict['assoc_3fhl_num'] > 0:
            m = (cat_3fhl['Source_Name'] == row_dict['name_3fhl'])
            row_dict['assoc'] = cat_3fhl[m]['ASSOC1'][0]
            row_dict['class'] = cat_3fhl[m]['CLASS'][0]
            row_dict['redshift'] = cat_3fhl[m]['Redshift'][0]
            row_dict['nupeak'] = cat_3fhl[m]['NuPeak_obs'][0]
            if np.log10(row_dict['nupeak']) < 14:
                row_dict['class_sed'] = 'LSP'
            elif np.log10(row_dict['nupeak']) > 15:
                row_dict['class_sed'] = 'HSP'
            elif np.log10(row_dict['nupeak']) > 14 or np.log10(row_dict['nupeak']) < 15:
                row_dict['class_sed'] = 'ISP'                
        else:
            row_dict['assoc'] = ''
            row_dict['class'] = ''
            row_dict['redshift'] = np.nan
            row_dict['nupeak'] = np.nan
                                    
        row_dict['spectrum_type'] = src_data[0]['SpectrumType']            
        row_dict['fit_mean_sep'] = fit_mean_sep
        row_dict['fit_min_sep'] = fit_min_sep
        row_dict['fit_idx'] = fit_idx
        row_dict['fit_idx_halo'] = fit_idx_halo
        row_dict['fit_idx_ext_gauss'] = fit_idx_ext_gauss
        row_dict['fit_idx_ext_disk'] = fit_idx_ext_disk
        if row_dict['fit_ext_model'] == 'disk':
            row_dict['fit_idx_ext'] = fit_idx_ext_disk
        else:
            row_dict['fit_idx_ext'] = fit_idx_ext_gauss
                        
        from scipy.interpolate import UnivariateSpline

        fit_dlnl_interp = np.zeros(halo_scan_shape + (len(eflux_scan_pts),))
        fit_sed_dlnl_interp = np.zeros((halo_scan_shape[0],nebins,len(eflux_scan_pts)))
        fit_halo_dloglike = np.array(tab_halo_data[fit_idx_halo]['dloglike_scan']).reshape(halo_scan_shape + (-1,))
        fit_halo_eflux = np.array(tab_halo_data[fit_idx_halo]['eflux_scan']).reshape(halo_scan_shape + (-1,))
        fit_halo_dloglike -= fit_halo_dloglike[:,:,0][:,:,None]
        
        for i in range(halo_scan_shape[0]):

            if not len(halo_seds):
                continue
            
            for j in range(nebins):
                x = halo_seds[fit_idx_halo]['dlnl_eflux'][i,j]
                y = halo_seds[fit_idx_halo]['dlnl'][i,j]
                isort = np.argsort(x)
                sp = UnivariateSpline(x[isort], y[isort], k=2, s=0.001)
                ys = sp(eflux_scan_pts)
                fit_sed_dlnl_interp[i,j] = ys 

                if np.any(~np.isfinite(sp(eflux_scan_pts))):
                    print(i,j,fit_idx_halo,x, y)
                    print(i,j,sp(eflux_scan_pts))
                
            
            for j in range(halo_scan_shape[1]):                
                sp = UnivariateSpline(fit_halo_eflux[i,j],
                                      fit_halo_dloglike[i,j],k=2,s=0.001)
                fit_dlnl_interp[i,j] = sp(eflux_scan_pts)
                
        row_dict_lnl['fit_src_sed_scan_dlnl'] = sed_data[fit_idx_halo]['dloglike_scan']
        row_dict_lnl['fit_src_sed_scan_eflux'] = sed_data[fit_idx_halo]['norm_scan'].copy()
        row_dict_lnl['fit_src_sed_scan_eflux'] *= sed_data[fit_idx_halo]['ref_eflux'][:,None]
        
        row_dict_sed['name'] = row_dict['name']
        row_dict_sed['norm'] = sed_data[-1]['norm']
        row_dict_sed['norm_err'] = sed_data[-1]['norm_err']
        row_dict_sed['norm_errp'] = sed_data[-1]['norm_err_hi']
        row_dict_sed['norm_errn'] = sed_data[-1]['norm_err_lo']
        row_dict_sed['norm_ul'] = sed_data[-1]['norm_ul']
        row_dict_sed['norm_scan'] = sed_data[-1]['norm_scan']
        row_dict_sed['dloglike_scan'] = sed_data[-1]['dloglike_scan']
        row_dict_sed['ref_flux'] = sed_data[-1]['ref_flux']
        row_dict_sed['ref_eflux'] = sed_data[-1]['ref_eflux']
        row_dict_sed['ref_npred'] = sed_data[-1]['ref_npred']
        row_dict_sed['ref_dnde'] = sed_data[-1]['ref_dnde']
            
        row_dict_lnl['name'] = row_dict['name']        
        row_dict_lnl['fit_ext_gauss_scan_dlnl'] = ext_gauss_data[fit_idx_ext_gauss]['dloglike']
        row_dict_lnl['fit_ext_gauss_scan_psfhi_dlnl'] = ext_gauss_data_psfhi[fit_idx_ext_gauss]['dloglike']
        row_dict_lnl['fit_ext_gauss_scan_psflo_dlnl'] = ext_gauss_data_psflo[fit_idx_ext_gauss]['dloglike']
        row_dict_lnl['fit_ext_disk_scan_dlnl'] = ext_disk_data[fit_idx_ext_disk]['dloglike']
        row_dict_lnl['fit_ext_disk_scan_psfhi_dlnl'] = ext_disk_data_psfhi[fit_idx_ext_disk]['dloglike']
        row_dict_lnl['fit_ext_disk_scan_psflo_dlnl'] = ext_disk_data_psflo[fit_idx_ext_disk]['dloglike']
        row_dict_lnl['fit_halo_scan_dlnl'] = fit_dlnl_interp
        row_dict_lnl['fit_halo_sed_scan_dlnl'] = fit_sed_dlnl_interp
        
        tab.add_row([row_dict[k] for k in tab.columns])
        tab_lnl.add_row([row_dict_lnl[k] for k in tab_lnl.columns])
        tab_sed.add_row([row_dict_sed[k] for k in tab_sed.columns])
                

    m = tab['class']==''
    tab['class'][m] = 'unkn'

    output_cat = os.path.splitext(output)[0] + '_cat.fits'
    output_lnl = os.path.splitext(output)[0] + '_lnl.fits'
    output_sed = os.path.splitext(output)[0] + '_sed.fits'

    cols_dict = OrderedDict()
    if len(tab_halo_data):

        halo_r68 = np.array(tab_halo_data[0]['halo_width']).reshape(halo_scan_shape)
        halo_index = np.array(tab_halo_data[0]['halo_index']).reshape(halo_scan_shape)
        cols_dict['halo_scan_r68'] = dict(dtype='f8', format='%.3f',
                                            data=halo_r68[:,0][np.newaxis,:])
        cols_dict['halo_scan_index'] = dict(dtype='f8', format='%.3f',
                                            data=halo_index[0,:][np.newaxis,:])
        cols_dict['halo_scan_eflux'] = dict(dtype='f8', format='%.3f',
                                            data=eflux_scan_pts[np.newaxis,:])
    cols_dict['ext_r68'] = dict(dtype='f8', format='%.3f',
                                data=ext_r68[0][np.newaxis,:])
    
    tab_grid = Table([Column(name=k, **v) for k, v in cols_dict.items()])

    cols_dict = OrderedDict()
    cols_dict['e_min'] = dict(dtype='f8', format='%.3f',
                              data=sed_data[0]['e_min'],unit='MeV')
    cols_dict['e_max'] = dict(dtype='f8', format='%.3f',
                              data=sed_data[0]['e_max'],unit='MeV')
    cols_dict['e_ref'] = dict(dtype='f8', format='%.3f',
                              data=sed_data[0]['e_ctr'],unit='MeV')
    cols_dict['ref_flux'] = dict(dtype='f8', format='%.3f',
                                 data=sed_data[0]['ref_flux'],
                                 unit='ph / (cm2 s)')
    cols_dict['ref_eflux'] = dict(dtype='f8', format='%.3f',
                                  data=sed_data[0]['ref_eflux'],
                                  unit='MeV / (cm2 s)')
    cols_dict['ref_dnde'] = dict(dtype='f8', format='%.3f',
                                 data=sed_data[0]['ref_dnde'],
                                 unit='1 / (MeV cm2 s)')
        
    tab_ebounds = Table([Column(name=k, **v) for k, v in cols_dict.items()])
    
    hdulist = fits.HDUList()
    hdulist.append(fits.table_to_hdu(tab))
    hdulist.append(fits.table_to_hdu(tab_grid))
    hdulist.append(fits.table_to_hdu(tab_ebounds))
    hdulist[1].name = 'CATALOG'
    hdulist[2].name = 'SCAN_PARS'
    hdulist[3].name = 'EBOUNDS'    
    hdulist.writeto(output_cat,clobber=True)

    hdulist = fits.HDUList()
    hdulist.append(fits.table_to_hdu(tab))
    hdulist.append(fits.table_to_hdu(tab_sed))
    hdulist.append(fits.table_to_hdu(tab_lnl))
    hdulist.append(fits.table_to_hdu(tab_grid))
    hdulist.append(fits.table_to_hdu(tab_ebounds))
    hdulist[1].name = 'CATALOG'
    hdulist[2].name = 'SED'
    hdulist[3].name = 'LIKELIHOOD'
    hdulist[4].name = 'SCAN_PARS'
    hdulist[5].name = 'EBOUNDS'    
    hdulist.writeto(output_lnl,clobber=True)

    hdulist = fits.HDUList()
    hdulist.append(fits.table_to_hdu(tab_sed))
    hdulist.append(fits.table_to_hdu(tab_grid))
    hdulist.append(fits.table_to_hdu(tab_ebounds))
    hdulist[1].name = 'SED'
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
