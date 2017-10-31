import os
import sys
import copy
import pprint
import argparse
from collections import OrderedDict

import yaml

from scipy.ndimage.interpolation import map_coordinates
from scipy.ndimage.interpolation import shift
import numpy as np
from astropy.table import Table, Column, join
from astropy.io import fits
import matplotlib.pyplot as plt
import matplotlib

import fermipy.utils as utils
from haloanalysis.utils import create_mask

def extract_halo_prop(halo_scan_ts, halo_scan_eflux, halo_scan_eflux_err,
                      halo_scan_eflux_ul95, halo_scan_r68, halo_scan_index):

    shape = halo_scan_ts.shape
    halo_scan_r68 = halo_scan_r68[:,None]*np.ones(shape)
    halo_scan_index = halo_scan_index[None,:]*np.ones(shape)
    
    #halo_scan_ts = np.array(tab_halo_data[i]['ts']).reshape(shape)
    #halo_scan_eflux = np.array(tab_halo_data[i]['eflux']).reshape(shape)
    #halo_scan_eflux_err = np.array(tab_halo_data[i]['eflux_err']).reshape(shape)
    #halo_scan_eflux_ul95 = np.array(tab_halo_data[i]['eflux_ul95']).reshape(shape)
    #halo_scan_flux = np.array(tab_halo_data[i]['flux']).reshape(shape)
    #halo_scan_flux_err = np.array(tab_halo_data[i]['flux_err']).reshape(shape)
    #halo_scan_flux_ul95 = np.array(tab_halo_data[i]['flux_ul95']).reshape(shape)
    halo_scan_idx = np.unravel_index(np.argmax(halo_scan_ts),shape)
    halo_scan_max_ts = halo_scan_ts[halo_scan_idx]
    halo_scan_max_eflux = halo_scan_eflux[halo_scan_idx]
    halo_scan_max_eflux_ul95 = halo_scan_eflux_ul95[halo_scan_idx]
    #halo_scan_max_flux = halo_scan_flux[halo_scan_idx]
    #halo_scan_max_flux_ul95 = halo_scan_flux_ul95.flat[halo_scan_idx]
    halo_scan_max_r68 = halo_scan_r68[halo_scan_idx]
    halo_scan_max_index = halo_scan_index[halo_scan_idx] 

    o = {}
    
    if halo_scan_max_ts > 2.0:

        ix, iy = np.unravel_index(np.argmax(halo_scan_ts),halo_scan_ts.shape)
        pfit = utils.fit_parabola(halo_scan_ts,ix,iy,dpix=2)
        
        print(ix,iy)
        pprint.pprint(pfit)
        if pfit['fit_success']:
            x0, y0 = pfit['x0'],pfit['y0']
        else:
            x0, y0 = float(ix), float(iy)
        
        halo_fit_ts = map_coordinates(halo_scan_ts,np.array([[x0],[y0]]))
        halo_fit_eflux = map_coordinates(halo_scan_eflux,np.array([[x0],[y0]]))
        halo_fit_eflux_err = map_coordinates(halo_scan_eflux_err,np.array([[x0],[y0]]))
        halo_fit_eflux_ul95 = map_coordinates(halo_scan_eflux_ul95,np.array([[x0],[y0]]))
        #halo_fit_flux = map_coordinates(halo_scan_flux,np.array([[x0],[y0]]))
        #halo_fit_flux_err = map_coordinates(halo_scan_flux_err,np.array([[x0],[y0]]))
        #halo_fit_flux_ul95 = map_coordinates(halo_scan_flux_ul95,np.array([[x0],[y0]]))
        halo_fit_r68 = map_coordinates(halo_scan_r68,np.array([[x0],[y0]]))
        halo_fit_index = map_coordinates(halo_scan_index,np.array([[x0],[y0]]))

        # Interpolation
        x = np.linspace(0,shape[0]-1,shape[0])
        y = np.linspace(0,shape[1]-1,shape[1])
        interp_ts_x = map_coordinates(halo_scan_ts,np.array([x,y0*np.ones_like(x)]))
        interp_ts_y = map_coordinates(halo_scan_ts,np.array([x0*np.ones_like(y),y]))
        halo_r68_lims = utils.get_parameter_limits(halo_scan_r68[:,iy], 0.5*interp_ts_x)
        halo_index_lims = utils.get_parameter_limits(halo_scan_index[ix,:], 0.5*interp_ts_y)
        halo_fit_r68_err = halo_r68_lims['err']
        halo_fit_index_err = halo_index_lims['err']

        o['ts'] = halo_fit_ts
        o['eflux'] = halo_fit_eflux
        o['eflux_err'] = halo_fit_eflux_err
        o['eflux_ul95'] = halo_fit_eflux_ul95
        #o['flux'] = halo_fit_flux
        #o['flux_err'] = halo_fit_flux_err
        #o['flux_ul95'] = halo_fit_flux_ul95
        o['r68'] = halo_fit_r68
        o['r68_err'] = halo_fit_r68_err
        o['index'] = halo_fit_index
        o['index_err'] = halo_fit_index_err
    else:
        o['ts'] = halo_scan_max_ts
        o['eflux'] = halo_scan_max_eflux
        o['eflux_err'] = np.nan
        o['eflux_ul95'] = halo_scan_max_eflux_ul95
        #o['flux'] = halo_scan_max_flux
        #o['flux_ul95'] = halo_scan_max_flux_ul95
        o['r68'] = halo_scan_max_r68
        o['r68_err'] = np.nan
        o['index'] = halo_scan_max_index
        o['index_err'] = np.nan

    pprint.pprint(o)
    return o
        
def main():

    usage = "usage: %(prog)s"
    description = "Perform stacking analysis."
    parser = argparse.ArgumentParser(usage=usage,description=description)

    parser.add_argument('--output', default = 'stack.fits')
    parser.add_argument('--prefix', default = '')
    parser.add_argument('--ext_colname', default = 'fit_ext_gauss_scan_dlnl')
    parser.add_argument('--targets', default = None, required=True)
    parser.add_argument('--start_fraction', default = None, type=float)
    parser.add_argument('--stop_fraction', default = None, type=float)
    parser.add_argument('--make_plots', default = False, action='store_true')

    parser.add_argument('tables', nargs='+', default = None,
                        help='Run stacking analysis.')

    args = parser.parse_args()

    targets = yaml.load(open(args.targets))
    tab0 = Table.read(args.tables[0],'LIKELIHOOD')
    tab1 = Table.read(args.tables[0],'CATALOG')
    tab0 = join(tab0,tab1)
    #print(len(tab0))

    #tab0 = Table.read('tevsources_tab.fits')
    
    ext_scan_r68 = np.array(Table.read(args.tables[0],hdu='SCAN_PARS')['ext_r68'])[0]
    halo_scan_r68 = np.array(Table.read(args.tables[0],hdu='SCAN_PARS')['halo_scan_r68'])[0]
    halo_scan_index = np.array(Table.read(args.tables[0],hdu='SCAN_PARS')['halo_scan_index'])[0]    
    halo_scan_eflux = np.array(Table.read(args.tables[0],hdu='SCAN_PARS')['halo_scan_eflux'])[0]
    halo_scan_fhalo = halo_scan_eflux/1E-6
    halo_scan_shape = (len(halo_scan_r68),len(halo_scan_index))
    halo_const_ts = np.zeros(halo_scan_shape)
    halo_const_eflux = np.zeros(halo_scan_shape)
    halo_const_eflux_err = np.zeros(halo_scan_shape)
    halo_const_eflux_ul95 = np.zeros(halo_scan_shape)
    halo_fhalo_ts = np.zeros(halo_scan_shape)
    halo_fhalo_eflux = np.zeros(halo_scan_shape)
    halo_fhalo_eflux_err = np.zeros(halo_scan_shape)
    halo_fhalo_eflux_ul95 = np.zeros(halo_scan_shape)
    
    targets_mask = {}
    
    for k,v in targets.items():
        targets_mask[k] = create_mask(tab0,v)
    
    cols_dict = OrderedDict()
    cols_dict['name'] = dict(dtype='S32', format='%s',description='Target Set')
    cols_dict['nobj'] = dict(dtype='i8')
    
    cols_dict['halo_const_scan_ts'] = dict(dtype='f8', format='%.2f',shape=halo_scan_shape)
    #cols_dict['halo_const_scan_dlnl'] = dict(dtype='f8', format='%.2f',shape=halo_scan_shape+(61,))
    cols_dict['halo_const_scan_eflux'] = dict(dtype='f8', format='%.4g',shape=halo_scan_shape)
    cols_dict['halo_const_scan_eflux_err'] = dict(dtype='f8', format='%.4g',shape=halo_scan_shape)
    cols_dict['halo_const_scan_eflux_ul95'] = dict(dtype='f8', format='%.4g',shape=halo_scan_shape)
    cols_dict['halo_const_ts'] = dict(dtype='f8', format='%.2f')
    cols_dict['halo_const_norm'] = dict(dtype='f8', format='%.3f')
    cols_dict['halo_const_norm_err'] = dict(dtype='f8', format='%.3f')    
    cols_dict['halo_const_index'] = dict(dtype='f8', format='%.3f')
    cols_dict['halo_const_index_err'] = dict(dtype='f8', format='%.3f')
    cols_dict['halo_const_r68'] = dict(dtype='f8', format='%.3f')
    cols_dict['halo_const_r68_err'] = dict(dtype='f8', format='%.3f')
    
    cols_dict['halo_fhalo_scan_ts'] = dict(dtype='f8', format='%.2f',shape=halo_scan_shape)
    #cols_dict['halo_fhalo_scan_dlnl'] = dict(dtype='f8', format='%.2f',shape=halo_scan_shape+(61,))
    cols_dict['halo_fhalo_scan_eflux'] = dict(dtype='f8', format='%.4g',shape=halo_scan_shape)
    cols_dict['halo_fhalo_scan_eflux_err'] = dict(dtype='f8', format='%.4g',shape=halo_scan_shape)
    cols_dict['halo_fhalo_scan_eflux_ul95'] = dict(dtype='f8', format='%.4g',shape=halo_scan_shape)
    cols_dict['halo_fhalo_ts'] = dict(dtype='f8', format='%.2f')
    cols_dict['halo_fhalo_norm'] = dict(dtype='f8', format='%.3f')
    cols_dict['halo_fhalo_norm_err'] = dict(dtype='f8', format='%.3f')    
    cols_dict['halo_fhalo_index'] = dict(dtype='f8', format='%.3f')
    cols_dict['halo_fhalo_index_err'] = dict(dtype='f8', format='%.3f')
    cols_dict['halo_fhalo_r68'] = dict(dtype='f8', format='%.3f')
    cols_dict['halo_fhalo_r68_err'] = dict(dtype='f8', format='%.3f')

    cols_dict['ext_ts_ext'] = dict(dtype='f8', format='%.2f')
    cols_dict['ext_sys_ts_ext'] = dict(dtype='f8', format='%.2f')
    cols_dict['ext_r68'] = dict(dtype='f8', format='%.4f')
    cols_dict['ext_r68_err'] = dict(dtype='f8', format='%.4f')
    cols_dict['ext_r68_sys_err'] = dict(dtype='f8', format='%.4f')
    cols_dict['ext_r68_ul95'] = dict(dtype='f8', format='%.4f')
    cols_dict['ext_r68_sys_ul95'] = dict(dtype='f8', format='%.4f')
    cols_dict['ext_scan_dlnl'] = dict(dtype='f8', format='%.4f',shape=ext_scan_r68.shape)
    cols_dict['ext_scan_r68'] = dict(dtype='f8', format='%.4f',shape=ext_scan_r68.shape)
    
    cols_dict['eflux10000'] = dict(dtype='f8', format='%.4g')
    cols_dict['eflux1000'] = dict(dtype='f8', format='%.4g')
    cols_dict['flux10000'] = dict(dtype='f8', format='%.4g')
    cols_dict['flux1000'] = dict(dtype='f8', format='%.4g')
    
    cols = [Column(name=k, **v) for k,v in cols_dict.items()]
    tab_out = Table(cols)
       
    row_dict = {}
    object_ext_dlnl = {}
    object_halo_const_dlnl = {}
    object_halo_fhalo_dlnl = {}
    object_ext_ts_ext = {}
    object_halo_ts = {}
    object_name = {}
    object_halo_fhalo_dlnl = {}
    
    for name, v in sorted(targets.items()):
        
        m = targets_mask[name]

        m &= ((tab0['codename'] != '3fgl_j0534.5+2201') &
              (tab0['codename'] != '3fgl_j0007.0+7302') &
              (tab0['codename'] != 'fhes_j0001.3+6841') &
              (tab0['codename'] != 'fhes_j0007.1+7303e_off') &
              (tab0['codename'] != 'fhes_j0042.5+4117e') &
              (tab0['codename'] != 'fhes_j0238.0+5238e') &
              (tab0['codename'] != '3fgl_j0425.8+5600') &
              (tab0['codename'] != 'fhes_j0429.3+3529') &
              (tab0['codename'] != 'fhes_j0940.8-6103e') &
              (tab0['codename'] != '3fgl_j1209.1-5224') &
              (tab0['codename'] != 'fhes_j1741.4-3857e') &
              (tab0['codename'] != 'fhes_j2125.8+5833e') &
              (tab0['codename'] != 'fhes_j2208.4+6501') &
              (tab0['codename'] != 'fhes_j2309.0+5429e') 
              )
        m &= (tab0['fit_ext_ts_ext'] < 16.0)
        
        tab = tab0[m]

        if args.start_fraction is not None:
            istart = int(float(len(tab))*args.start_fraction)
            istop = int(float(len(tab))*args.stop_fraction)
            tab = tab[istart:istop]

        #tab.write(name + '_tab.fits',overwrite=True)
        print name, len(tab)
        
        # These should eventually be extracted from the FITS file
        halo_eflux_dlnl = np.array(tab['fit_halo_scan_dlnl'])        
        halo_eflux_dlnl = halo_eflux_dlnl - halo_eflux_dlnl[:,:,:,:1]
        
        object_ext_dlnl[name] = tab[args.ext_colname]
        object_halo_const_dlnl[name] = halo_eflux_dlnl
        object_ext_ts_ext[name] = np.array(tab['fit_ext_ts_ext'])
        object_halo_ts[name] = tab['fit_halo_ts']
        object_name[name] = tab['name']

        # Convert offset to pixel coordinates
        eflux_offset = -np.log10(tab['eflux1000']/1E-6)*10.
        halo_fhalo_dlnl = np.zeros(halo_eflux_dlnl.shape)
        for i in range(len(halo_eflux_dlnl)):
            halo_fhalo_dlnl[i] = shift(halo_eflux_dlnl[i],
                                           [0.0,0.0,eflux_offset[i]],
                                           mode='nearest',order=1)

        halo_fhalo_dlnl = halo_fhalo_dlnl - halo_eflux_dlnl[:,:,:,:1]

        object_halo_fhalo_dlnl[name] = halo_fhalo_dlnl
        
        # Sum over all sources
        halo_const_dlnl = np.sum(halo_eflux_dlnl,axis=0)
        halo_fhalo_dlnl = np.sum(halo_fhalo_dlnl,axis=0)
        ext_dlnl = np.sum(tab['fit_ext_gauss_scan_dlnl'],axis=0)
        ext_dlnl -= ext_dlnl[0]
        ext_psfhi_dlnl = np.sum(tab['fit_ext_gauss_scan_psfhi_dlnl'],axis=0)
        ext_psfhi_dlnl -= ext_psfhi_dlnl[0]
        ext_psflo_dlnl = np.sum(tab['fit_ext_gauss_scan_psflo_dlnl'],axis=0)
        ext_psflo_dlnl -= ext_psflo_dlnl[0]
                
        for i in range(len(halo_scan_r68)):
            for j in range(len(halo_scan_index)):

                lims = utils.get_parameter_limits(halo_scan_eflux,halo_const_dlnl[i,j])
                halo_const_ts[i,j] = 2.0*lims['lnlmax']
                halo_const_eflux_ul95[i,j] = lims['ul']
                halo_const_eflux[i,j] = lims['x0']
                halo_const_eflux_err[i,j] = lims['err']

                lims = utils.get_parameter_limits(halo_scan_fhalo,halo_fhalo_dlnl[i,j])
                halo_fhalo_ts[i,j] = 2.0*lims['lnlmax']
                halo_fhalo_eflux_ul95[i,j] = lims['ul']
                halo_fhalo_eflux[i,j] = lims['x0']
                halo_fhalo_eflux_err[i,j] = lims['err']
                
        #halo_ts = 2.0*np.max(halo_eflux_dlnl,axis=2)

        halo_const_fit = extract_halo_prop(halo_const_ts, halo_const_eflux, halo_const_eflux_err,
                                           halo_const_eflux_ul95, halo_scan_r68, halo_scan_index)
        halo_fhalo_fit = extract_halo_prop(halo_fhalo_ts, halo_fhalo_eflux, halo_fhalo_eflux_err,
                                           halo_fhalo_eflux_ul95, halo_scan_r68, halo_scan_index)

        #print(halo_const_fit)
        #print(halo_fhalo_fit)
        
        halo_const_max_idx = np.unravel_index(np.argmax(halo_const_ts),halo_scan_shape)
        halo_fhalo_max_idx = np.unravel_index(np.argmax(halo_fhalo_ts),halo_scan_shape)
                
        row_dict['name'] = name
        row_dict['nobj'] = len(tab)

        row_dict['flux10000'] = np.sum(tab['flux10000'])
        row_dict['flux1000'] = np.sum(tab['flux1000'])
        row_dict['eflux10000'] = np.sum(tab['eflux10000'])
        row_dict['eflux1000'] = np.sum(tab['eflux1000'])

        row_dict['halo_const_scan_dlnl'] = halo_const_dlnl
        row_dict['halo_const_scan_ts'] = halo_const_ts
        row_dict['halo_const_scan_eflux'] = halo_const_eflux
        row_dict['halo_const_scan_eflux_err'] = halo_const_eflux_err
        row_dict['halo_const_scan_eflux_ul95'] = halo_const_eflux_ul95
        row_dict['halo_const_ts'] = halo_const_fit['ts']
        row_dict['halo_const_norm'] = halo_const_fit['eflux']
        row_dict['halo_const_norm_err'] = halo_const_fit['eflux_err']
        row_dict['halo_const_index'] = halo_const_fit['index']
        row_dict['halo_const_index_err'] = halo_const_fit['index_err']
        row_dict['halo_const_r68'] = halo_const_fit['r68']
        row_dict['halo_const_r68_err'] = halo_const_fit['r68_err']

        row_dict['halo_fhalo_scan_dlnl'] = halo_fhalo_dlnl
        row_dict['halo_fhalo_scan_ts'] = halo_fhalo_ts
        row_dict['halo_fhalo_scan_eflux'] = halo_fhalo_eflux
        row_dict['halo_fhalo_scan_eflux_err'] = halo_fhalo_eflux_err
        row_dict['halo_fhalo_scan_eflux_ul95'] = halo_fhalo_eflux_ul95
        row_dict['halo_fhalo_ts'] = halo_fhalo_fit['ts']
        row_dict['halo_fhalo_norm'] = halo_fhalo_fit['eflux']
        row_dict['halo_fhalo_norm_err'] = halo_fhalo_fit['eflux_err']
        row_dict['halo_fhalo_index'] = halo_fhalo_fit['index']
        row_dict['halo_fhalo_index_err'] = halo_fhalo_fit['index_err']
        row_dict['halo_fhalo_r68'] = halo_fhalo_fit['r68']
        row_dict['halo_fhalo_r68_err'] = halo_fhalo_fit['r68_err']
        
        row_dict['ext_scan_dlnl'] = ext_dlnl
        row_dict['ext_scan_r68'] = ext_scan_r68

        ext_lims = utils.get_parameter_limits(ext_scan_r68,ext_dlnl)
        ext_psfhi_lims = utils.get_parameter_limits(ext_scan_r68,ext_psfhi_dlnl)
        ext_psflo_lims = utils.get_parameter_limits(ext_scan_r68,ext_psflo_dlnl)

        #print(ext_lims)
        #print(ext_psfhi_lims)
        #print(ext_psflo_lims)
        
        row_dict['ext_ts_ext'] = 2.0*(ext_lims['lnlmax'])
        row_dict['ext_sys_ts_ext'] = min(2.0*(ext_lims['lnlmax']),
                                         2.0*(ext_psfhi_lims['lnlmax']),
                                         2.0*(ext_psflo_lims['lnlmax']))                                         
        row_dict['ext_r68'] = ext_lims['x0']
        row_dict['ext_r68_ul95'] = ext_lims['ul']
        row_dict['ext_r68_err'] = ext_lims['err']
        row_dict['ext_r68_sys_err'] = 0.5*np.abs(ext_psfhi_lims['x0']-ext_psflo_lims['x0'])
        row_dict['ext_r68_sys_ul95'] = max(ext_lims['ul'],
                                           ext_psfhi_lims['ul'],
                                           ext_psflo_lims['ul'])
        
        tab_out.add_row([row_dict[k] for k in tab_out.columns])

#    tab_out.write(args.output,format='fits',overwrite=True)

    hdulist = fits.HDUList()
    hdulist.append(fits.table_to_hdu(tab_out))
    hdulist.append(fits.table_to_hdu(Table.read(args.tables[0],hdu='SCAN_PARS')))
    hdulist[1].name = 'SCAN_DATA'    
    hdulist.writeto(args.output, clobber=True)
        
    np.save(os.path.splitext(args.output)[0] + '.npy',
            {'object_ext_dlnl' : object_ext_dlnl,
             'object_ext_ts_ext' : object_ext_ts_ext,
             'object_name' : object_name,
             'object_halo_const_dlnl' : object_halo_const_dlnl,
             'object_halo_fhalo_dlnl' : object_halo_fhalo_dlnl,
             })
    
    if not args.make_plots:
        sys.exit(0)

    for row in tab_out:
        continue
        
        name = row['name']
        
        plt.figure()
        plt.gca().set_title(name)
        for i, idx in enumerate(halo_scan_index):
            plt.plot(halo_scan_r68,row['fit_halo_scan_ts'][:,i],
                     label='%.2f'%idx,
                     color=matplotlib.cm.spectral(float(i)/len(halo_scan_index)))

        plt.gca().set_xscale('log')
        plt.gca().legend(frameon=False)
        plt.gca().set_xlabel('R$_{68}$ [deg]')
        plt.gca().set_ylabel('TS')
        plt.gca().set_ylim(0,16)
        plt.savefig('composite_halo_scan_ts_%s.png'%name)
                     
#        im = plt.imshow(row['fit_halo_scan_ts'].T,origin='lower',interpolation='nearest')
#        plt.colorbar(im)
#        plt.savefig('composite_halo_scan_ts_%s.png'%name)
#        plt.figure()
#        plt.gca().set_title(name)
#        im = plt.imshow(row['fit_halo_scan_eflux_ul95'].T,origin='lower',interpolation='nearest')
#        plt.colorbar(im)
#        plt.savefig('composite_halo_scan_eflux_%s.png'%name)


        
    #fig_ext = plt.figure()
    #for row in tab_out:

    #    plt.figure(fig_ext.number)
    #    plt.plot(ext_r68,row['fit_ext_dlnl'],label=row['name'])

    #plt.gca().legend(frameon=False)
    #plt.gca().set_ylim(-10,10)
    #plt.gca().set_xscale('log')


    
    for row in tab_out:

        name = row['name']

        plt.figure()
        plt.plot(ext_scan_r68,2.0*row['ext_scan_dlnl'],label=name,linewidth=1.5)
        plt.axvline(row['ext_r68_ul95'],label=name + ' 95% UL',linewidth=1.0,color='b',
                    linestyle='--')
        for i in range(object_ext_dlnl[name].shape[0]):
            plt.plot(ext_scan_r68,2.0*object_ext_dlnl[name][i],color='grey',alpha=0.3)

        i = np.argmax(object_ext_ts_ext[name])
        plt.plot(ext_scan_r68,2.0*object_ext_dlnl[name][i],
                 label=object_name[name][i],
                 linewidth=1.5)

        plt.axhline(0.0,color='k',linestyle='--')
        plt.gca().legend(frameon=False)
        plt.gca().set_ylim(-10,30)
        plt.gca().set_xlim(0.00316,2.0)
        plt.gca().set_xscale('log')
        plt.gca().set_xlabel('R$_{68}$ [deg]')
        plt.gca().set_ylabel('2 $\\times$ Delta LogLikelihood')
        plt.savefig('%scomposite_ext_%s.png'%(args.prefix,name))

#        1/0
        
        plt.figure()
        plt.plot(halo_scan_eflux,2.0*row['halo_const_scan_dlnl'][4,2],label=name,linewidth=1.5)
        plt.axvline(row['halo_const_scan_eflux_ul95'][4,2],label=name + ' 95% UL',linewidth=1.0,color='b',
                    linestyle='--')
        for i in range(object_ext_dlnl[name].shape[0]):
            plt.plot(halo_scan_eflux,2.0*object_halo_const_dlnl[name][i,4,2],color='grey',alpha=0.3)

        i = np.argmax(np.max(object_halo_const_dlnl[name][:,4,2,:],axis=1))
        plt.plot(halo_scan_eflux,2.0*object_halo_const_dlnl[name][i,4,2],
                 label=object_name[name][i],
                 linewidth=1.5)

        plt.gca().text(0.05,0.8,
                       'Halo Model\nR68 = %.3f deg\nIndex = %.2f'%(halo_scan_r68[4],halo_scan_index[2]),
                       transform=plt.gca().transAxes)
        
        plt.axhline(0.0,color='k',linestyle='--')
        plt.gca().legend(frameon=False)
        plt.gca().set_ylim(-10,20)
        plt.gca().set_xscale('log')
        plt.gca().set_xlabel('Energy Flux (E = 1-316 GeV) [MeV cm$^{-2}$ s$^{-1}$]')
        plt.gca().set_ylabel('2 $\\times$ Delta LogLikelihood')

        plt.savefig('%scomposite_halo_const_%s.png'%(args.prefix,name))

        plt.figure()
        plt.plot(halo_scan_fhalo,2.0*row['halo_fhalo_scan_dlnl'][4,2],label=name,linewidth=1.5)
        plt.axvline(row['halo_fhalo_scan_eflux_ul95'][4,2],label=name + ' 95% UL',linewidth=1.0,color='b',
                    linestyle='--')
        for i in range(object_ext_dlnl[name].shape[0]):
            plt.plot(halo_scan_fhalo,2.0*object_halo_fhalo_dlnl[name][i,4,2],color='grey',alpha=0.3)

        i = np.argmax(np.max(object_halo_fhalo_dlnl[name][:,4,2,:],axis=1))
        plt.plot(halo_scan_fhalo,2.0*object_halo_fhalo_dlnl[name][i,4,2],
                 label=object_name[name][i],
                 linewidth=1.5)

        plt.gca().text(0.05,0.8,
                       'Halo Model\nR68 = %.3f deg\nIndex = %.2f'%(halo_scan_r68[4],halo_scan_index[2]),
                       transform=plt.gca().transAxes)
        
        plt.axhline(0.0,color='k',linestyle='--')
        plt.gca().legend(frameon=False)
        plt.gca().set_ylim(-10,20)
        plt.gca().set_xscale('log')
        plt.gca().set_xlabel('Halo Fraction')
        plt.gca().set_ylabel('2 $\\times$ Delta LogLikelihood')

        plt.savefig('%scomposite_halo_fhalo_%s.png'%(args.prefix,name))
        
if __name__ == "__main__":
    main()
