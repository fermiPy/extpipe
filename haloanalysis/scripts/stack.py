import os
import sys
import copy
import argparse
from collections import OrderedDict

import yaml

from scipy.ndimage.interpolation import shift
import numpy as np
from astropy.table import Table, Column, join
from astropy.io import fits
import matplotlib.pyplot as plt
import matplotlib

import fermipy.utils as utils
from haloanalysis.utils import create_mask

def main():

    usage = "usage: %(prog)s"
    description = "Perform stacking analysis."
    parser = argparse.ArgumentParser(usage=usage,description=description)

    parser.add_argument('--output', default = 'stack.fits')
    parser.add_argument('--targets', default = None, required=True)
    parser.add_argument('--make_plots', default = False, action='store_true')

    parser.add_argument('tables', nargs='+', default = None,
                        help='Run stacking analysis.')

    args = parser.parse_args()

    targets = yaml.load(open(args.targets))
    tab0 = Table.read(args.tables[0])
    tab1 = Table.read(args.tables[1])
    tab0 = join(tab0,tab1)

#    tab0 = tab0[:1200]
    
    ext_width = np.array(Table.read(args.tables[0],hdu='SCAN_PARS')['ext_width'])[0]
    halo_scan_width = np.array(Table.read(args.tables[0],hdu='SCAN_PARS')['halo_scan_width'])[0]
    halo_scan_index = np.array(Table.read(args.tables[0],hdu='SCAN_PARS')['halo_scan_index'])[0]    
    halo_scan_eflux = np.array(Table.read(args.tables[0],hdu='SCAN_PARS')['halo_scan_eflux'])[0]
    halo_scan_fhalo = halo_scan_eflux/1E-6
    halo_scan_shape = (len(halo_scan_width),len(halo_scan_index))
    halo_const_ts = np.zeros(halo_scan_shape)
    halo_const_eflux_ul95 = np.zeros(halo_scan_shape)
    halo_fhalo_ts = np.zeros(halo_scan_shape)
    halo_fhalo_eflux_ul95 = np.zeros(halo_scan_shape)
    
    targets_mask = {}
    
    for k,v in targets.items():
        targets_mask[k] = create_mask(tab0,v)
    
    cols_dict = OrderedDict()
    cols_dict['name'] = dict(dtype='S20', format='%s',description='Target Set')        
    cols_dict['halo_const_scan_ts'] = dict(dtype='f8', format='%.2f',shape=halo_scan_shape)
    cols_dict['halo_const_scan_dlnl'] = dict(dtype='f8', format='%.2f',shape=halo_scan_shape+(41,))
    cols_dict['halo_const_scan_eflux_ul95'] = dict(dtype='f8', format='%.4g',shape=halo_scan_shape)
    cols_dict['halo_const_ts'] = dict(dtype='f8', format='%.2f')
    cols_dict['halo_const_index'] = dict(dtype='f8', format='%.3f')
    cols_dict['halo_const_width'] = dict(dtype='f8', format='%.3f')
    cols_dict['halo_fhalo_scan_ts'] = dict(dtype='f8', format='%.2f',shape=halo_scan_shape)
    cols_dict['halo_fhalo_scan_dlnl'] = dict(dtype='f8', format='%.2f',shape=halo_scan_shape+(41,))
    cols_dict['halo_fhalo_scan_eflux_ul95'] = dict(dtype='f8', format='%.4g',shape=halo_scan_shape)
    cols_dict['halo_fhalo_ts'] = dict(dtype='f8', format='%.2f')
    cols_dict['halo_fhalo_index'] = dict(dtype='f8', format='%.3f')
    cols_dict['halo_fhalo_width'] = dict(dtype='f8', format='%.3f')

    cols_dict['ext_width'] = dict(dtype='f8', format='%.4f',shape=(len(ext_width),))
    
    cols_dict['ext_ts'] = dict(dtype='f8', format='%.2f')
    cols_dict['ext_mle'] = dict(dtype='f8', format='%.4f')
    cols_dict['ext_ul95'] = dict(dtype='f8', format='%.4f')
    cols_dict['ext_dlnl'] = dict(dtype='f8', format='%.4f',shape=(len(ext_width),))

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
    object_ext_ts = {}
    object_halo_ts = {}
    object_name = {}
    object_halo_fhalo_dlnl = {}
    
    for name, v in sorted(targets.items()):
        
        m = targets_mask[name]

        m &= (tab0['name'] != '3FGL J0534.5+2201')
        
        tab = tab0[m]

        print name, len(tab)
        
        # These should eventually be extracted from the FITS file
        halo_eflux_dlnl = np.array(tab['fit_halo_scan_dlnl'])        
        halo_eflux_dlnl = halo_eflux_dlnl - halo_eflux_dlnl[:,:,:,:1]
        
        object_ext_dlnl[name] = tab['fit_ext_scan_dlnl']
        object_halo_const_dlnl[name] = halo_eflux_dlnl
        object_ext_ts[name] = tab['fit_ext_ts']
        object_halo_ts[name] = tab['fit_halo_ts']
        object_name[name] = tab['name']

        eflux_offset = -np.log10(tab['eflux1000']/1E-6)*10.

        halo_fhalo_dlnl = np.zeros(halo_eflux_dlnl.shape)
        for i in range(len(halo_eflux_dlnl)):
#            print i, eflux_offset[i]
            halo_fhalo_dlnl[i] = shift(halo_eflux_dlnl[i],
                                           [0.0,0.0,eflux_offset[i]],
                                           mode='nearest',order=1)

        halo_fhalo_dlnl = halo_fhalo_dlnl - halo_eflux_dlnl[:,:,:,:1]

        object_halo_fhalo_dlnl[name] = halo_fhalo_dlnl
        
        # Sum over all sources
        halo_const_dlnl = np.sum(halo_eflux_dlnl,axis=0)
        halo_fhalo_dlnl = np.sum(halo_fhalo_dlnl,axis=0)
        ext_dlnl = np.sum(tab['fit_ext_scan_dlnl'],axis=0)
        ext_dlnl -= ext_dlnl[0]
                
        for i in range(len(halo_scan_width)):
            for j in range(len(halo_scan_index)):

                lims = utils.get_parameter_limits(halo_scan_eflux,halo_const_dlnl[i,j])
                halo_const_ts[i,j] = 2.0*lims['lnlmax']
                halo_const_eflux_ul95[i,j] = lims['ul']

                lims = utils.get_parameter_limits(halo_scan_fhalo,halo_fhalo_dlnl[i,j])
                halo_fhalo_ts[i,j] = 2.0*lims['lnlmax']
                halo_fhalo_eflux_ul95[i,j] = lims['ul']
                
        #halo_ts = 2.0*np.max(halo_eflux_dlnl,axis=2)

        halo_const_max_idx = np.unravel_index(np.argmax(halo_const_ts),halo_scan_shape)
        halo_fhalo_max_idx = np.unravel_index(np.argmax(halo_fhalo_ts),halo_scan_shape)
                
        row_dict['name'] = name

        row_dict['flux10000'] = np.sum(tab['flux10000'])
        row_dict['flux1000'] = np.sum(tab['flux1000'])
        row_dict['eflux10000'] = np.sum(tab['eflux10000'])
        row_dict['eflux1000'] = np.sum(tab['eflux1000'])


        
        row_dict['halo_const_scan_dlnl'] = halo_const_dlnl
        row_dict['halo_const_scan_ts'] = halo_const_ts
        row_dict['halo_const_scan_eflux_ul95'] = halo_const_eflux_ul95
        row_dict['halo_const_ts'] = halo_const_ts[halo_const_max_idx]
        row_dict['halo_const_index'] = halo_scan_index[halo_const_max_idx[1]]
        row_dict['halo_const_width'] = halo_scan_width[halo_const_max_idx[0]]

        row_dict['halo_fhalo_scan_dlnl'] = halo_fhalo_dlnl
        row_dict['halo_fhalo_scan_ts'] = halo_fhalo_ts
        row_dict['halo_fhalo_scan_eflux_ul95'] = halo_fhalo_eflux_ul95
        row_dict['halo_fhalo_ts'] = halo_fhalo_ts[halo_fhalo_max_idx]
        row_dict['halo_fhalo_index'] = halo_scan_index[halo_fhalo_max_idx[1]]
        row_dict['halo_fhalo_width'] = halo_scan_width[halo_fhalo_max_idx[0]]
        row_dict['ext_dlnl'] = ext_dlnl

        ext_lims = utils.get_parameter_limits(ext_width,ext_dlnl)

        row_dict['ext_ts'] = 2.0*(ext_lims['lnlmax'])
        row_dict['ext_ul95'] = ext_lims['ul']
        row_dict['ext_mle'] = ext_lims['x0']
        row_dict['ext_width'] = ext_width
        
        tab_out.add_row([row_dict[k] for k in tab_out.columns])

#    tab_out.write(args.output,format='fits',overwrite=True)

    hdulist = fits.HDUList()
    hdulist.append(fits.table_to_hdu(tab_out))
    hdulist.append(fits.table_to_hdu(Table.read(args.tables[0],hdu='SCAN_PARS')))
    hdulist[1].name = 'SCAN_DATA'    
    hdulist.writeto(args.output, clobber=True)
        
    np.save(os.path.splitext(args.output)[0] + '.npy',
            {'object_ext_dlnl' : object_ext_dlnl,
             'object_ext_ts' : object_ext_ts,
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
            plt.plot(halo_scan_width,row['fit_halo_scan_ts'][:,i],
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
    #    plt.plot(ext_width,row['fit_ext_dlnl'],label=row['name'])

    #plt.gca().legend(frameon=False)
    #plt.gca().set_ylim(-10,10)
    #plt.gca().set_xscale('log')


    
    for row in tab_out:

        name = row['name']

        plt.figure()
        plt.plot(ext_width,2.0*row['ext_dlnl'],label=name,linewidth=1.5)
        plt.axvline(row['ext_ul95'],label=name + ' 95% UL',linewidth=1.0,color='b',
                    linestyle='--')
        for i in range(object_ext_dlnl[name].shape[0]):
            plt.plot(ext_width,2.0*object_ext_dlnl[name][i],color='grey',alpha=0.3)

        i = np.argmax(object_ext_ts[name])
        plt.plot(ext_width,2.0*object_ext_dlnl[name][i],
                 label=object_name[name][i],
                 linewidth=1.5)

        plt.axhline(0.0,color='k',linestyle='--')
        plt.gca().legend(frameon=False)
        plt.gca().set_ylim(-10,100)
        plt.gca().set_xlim(0.00316,1.0)
        plt.gca().set_xscale('log')
        plt.gca().set_xlabel('R$_{68}$ [deg]')
        plt.gca().set_ylabel('2 $\\times$ Delta LogLikelihood')
        plt.savefig('composite_ext_%s.png'%name)
        
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
                       'Halo Model\nR68 = %.3f deg\nIndex = %.2f'%(halo_scan_width[4],halo_scan_index[2]),
                       transform=plt.gca().transAxes)
        
        plt.axhline(0.0,color='k',linestyle='--')
        plt.gca().legend(frameon=False)
        plt.gca().set_ylim(-10,20)
        plt.gca().set_xscale('log')
        plt.gca().set_xlabel('Energy Flux (E = 1-316 GeV) [MeV cm$^{-2}$ s$^{-1}$]')
        plt.gca().set_ylabel('2 $\\times$ Delta LogLikelihood')

        plt.savefig('composite_halo_const_%s.png'%name)

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
                       'Halo Model\nR68 = %.3f deg\nIndex = %.2f'%(halo_scan_width[4],halo_scan_index[2]),
                       transform=plt.gca().transAxes)
        
        plt.axhline(0.0,color='k',linestyle='--')
        plt.gca().legend(frameon=False)
        plt.gca().set_ylim(-10,20)
        plt.gca().set_xscale('log')
        plt.gca().set_xlabel('Halo Fraction')
        plt.gca().set_ylabel('2 $\\times$ Delta LogLikelihood')

        plt.savefig('composite_halo_fhalo_%s.png'%name)
        
if __name__ == "__main__":
    main()
