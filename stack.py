import sys
import argparse
from collections import OrderedDict

import yaml

import numpy as np
from astropy.table import Table, Column, join
import matplotlib.pyplot as plt
import matplotlib

import fermipy.utils as utils


def create_mask(tab,target_def):
    """Create a table mask from a target definition."""

    m = np.empty(len(tab),dtype=bool); m.fill(True)

    for k,v in target_def.items():
        
        if isinstance(v,list):

            m0 = np.zeros(len(tab0),dtype=bool)

            for t in v:
                m0 |= (tab0[k] == t)

            m &= m0
            
        elif isinstance(v,dict):

            m0 = np.empty(len(tab0),dtype=bool)
            m0.fill(True)
            if 'min' in v:
                m0 &= (tab0[k] >= v['min'])

            if 'max' in v:
                m0 &= (tab0[k] <= v['max'])
            m &= m0
            
    return m

    

if __name__ == '__main__':

    usage = "usage: %(prog)s"
    description = "Perform stacking analysis."
    parser = argparse.ArgumentParser(usage=usage,description=description)

    parser.add_argument('--output', default = 'stack.fits')
    parser.add_argument('--targets', default = None, required=True)
    parser.add_argument('--make_plots', default = False, action='store_true')

    parser.add_argument('tables', nargs='+', default = None,
                        help='Run analyses in all subdirectories of this '
                        'directory.')

    args = parser.parse_args()

    targets = yaml.load(open(args.targets))
    tab0 = Table.read(args.tables[0])
    tab1 = Table.read(args.tables[1])
    tab0 = join(tab0,tab1)

#    tab0 = tab0[:1200]
    
    ext_width = np.array(Table.read(args.tables[0],hdu='EXT_WIDTH')['ext_width'])
    halo_scan_width = np.array(Table.read(args.tables[0],hdu='HALO_SCAN_WIDTH')['halo_scan_width'])
    halo_scan_index = np.array(Table.read(args.tables[0],hdu='HALO_SCAN_INDEX')['halo_scan_index'])
    halo_scan_shape = (len(halo_scan_width),len(halo_scan_index))
    eflux = np.logspace(-9,-5,41)
    halo_ts = np.zeros(halo_scan_shape)
    halo_eflux_ul95 = np.zeros(halo_scan_shape)
    
    targets_mask = {}
    
    for k,v in targets.items():
        targets_mask[k] = create_mask(tab0,v)
    
    cols_dict = OrderedDict()
    cols_dict['name'] = dict(dtype='S20', format='%s',description='Target Set')        
    cols_dict['fit_halo_scan_ts'] = dict(dtype='f8', format='%.2f',shape=(9,7))
    cols_dict['fit_halo_scan_dlnl'] = dict(dtype='f8', format='%.2f',shape=(9,7,41))
    cols_dict['fit_halo_scan_eflux_ul95'] = dict(dtype='f8', format='%.4g',shape=(9,7))

    cols_dict['fit_ext_ts'] = dict(dtype='f8', format='%.2f')
    cols_dict['fit_ext_ul95'] = dict(dtype='f8', format='%.2f')
    cols_dict['fit_ext_dlnl'] = dict(dtype='f8', format='%.2f',shape=(26))
    
    cols = [Column(name=k, **v) for k,v in cols_dict.items()]
    tab_out = Table(cols)
       
    row_dict = {}
    object_fit_ext_dlnl = {}
    object_fit_halo_dlnl = {}
    object_fit_ext_ts = {}
    object_fit_halo_ts = {}
    object_name = {}
    
    for name, v in sorted(targets.items()):
        
        m = targets_mask[name]
        tab = tab0[m]

        print name, len(tab)
        
        # These should eventually be extracted from the FITS file
        fit_halo_scan_dlnl = np.array(tab['fit_halo_scan_dlnl'])        
        fit_halo_scan_dlnl = fit_halo_scan_dlnl - fit_halo_scan_dlnl[:,:,:,:1]
        
        object_fit_ext_dlnl[name] = tab['fit_ext_scan_dlnl']
        object_fit_halo_dlnl[name] = fit_halo_scan_dlnl
        object_fit_ext_ts[name] = tab['fit_ext_ts']
        object_fit_halo_ts[name] = tab['fit_halo_ts']
        object_name[name] = tab['name']
        
        # Sum over all sources
        fit_halo_dlnl = np.sum(fit_halo_scan_dlnl,axis=0)
        fit_ext_dlnl =  np.sum(tab['fit_ext_scan_dlnl'],axis=0)
        fit_ext_dlnl -= fit_ext_dlnl[0]
                
        for i in range(9):
            for j in range(7):

                fit_halo_dlnl_sum = fit_halo_dlnl[i,j]
                lims = utils.get_parameter_limits(eflux,fit_halo_dlnl_sum)
                halo_ts[i,j] = 2.0*lims['lnlmax']
                halo_eflux_ul95[i,j] = lims['ul']

        #halo_ts = 2.0*np.max(fit_halo_dlnl,axis=2)
                
        row_dict['name'] = name
        row_dict['fit_halo_scan_dlnl'] = fit_halo_dlnl
        row_dict['fit_halo_scan_ts'] = halo_ts
        row_dict['fit_halo_scan_eflux_ul95'] = halo_eflux_ul95
        row_dict['fit_ext_dlnl'] = fit_ext_dlnl

        ext_lims = utils.get_parameter_limits(ext_width,fit_ext_dlnl)

        row_dict['fit_ext_ts'] = 2.0*(ext_lims['lnlmax'])
        row_dict['fit_ext_ul95'] = ext_lims['ul']
        
        tab_out.add_row([row_dict[k] for k in cols_dict.keys()])

    tab_out.write(args.output,format='fits',overwrite=True)

    if not args.make_plots:
        sys.exit(0)

    for row in tab_out:

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

    sys.exit(0)
        
    fig_ext = plt.figure()
    for row in tab_out:

        plt.figure(fig_ext.number)
        plt.plot(ext_width,row['fit_ext_dlnl'],label=row['name'])

    plt.gca().legend(frameon=False)
    plt.gca().set_ylim(-10,10)
    plt.gca().set_xscale('log')

    for row in tab_out:

        name = row['name']

        plt.figure()
        plt.plot(ext_width,2.0*row['fit_ext_dlnl'],label=name,linewidth=1.5)
        plt.axvline(row['fit_ext_ul95'],label=name + ' 95% UL',linewidth=1.0,color='b',
                    linestyle='--')
        for i in range(object_fit_ext_dlnl[name].shape[0]):
            plt.plot(ext_width,2.0*object_fit_ext_dlnl[name][i],color='grey',alpha=0.3)

        i = np.argmax(object_fit_ext_ts[name])
        plt.plot(ext_width,2.0*object_fit_ext_dlnl[name][i],
                 label=object_name[name][i],
                 linewidth=1.5)

        plt.axhline(0.0,color='k',linestyle='--')
        plt.gca().legend(frameon=False)
        plt.gca().set_ylim(-10,100)
        plt.gca().set_xlim(0.01,1.0)
        plt.gca().set_xscale('log')
        plt.gca().set_xlabel('R$_{68}$ [deg]')
        plt.gca().set_ylabel('2 $\\times$ Delta LogLikelihood')
        plt.savefig('composite_ext_%s.png'%name)

        plt.figure()
        plt.plot(eflux,2.0*row['fit_halo_scan_dlnl'][4,2],label=name,linewidth=1.5)
        plt.axvline(row['fit_halo_scan_eflux_ul95'][4,2],label=name + ' 95% UL',linewidth=1.0,color='b',
                    linestyle='--')
        for i in range(object_fit_ext_dlnl[name].shape[0]):
            plt.plot(eflux,2.0*object_fit_halo_dlnl[name][i,4,2],color='grey',alpha=0.3)

        i = np.argmax(np.max(object_fit_halo_dlnl[name][:,4,2,:],axis=1))
#        i = np.argmax(object_fit_halo_ts[name])
        plt.plot(eflux,2.0*object_fit_halo_dlnl[name][i,4,2],
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

        plt.savefig('composite_halo_%s.png'%name)
        
