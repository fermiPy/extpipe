import argparse
from collections import OrderedDict

import yaml

import numpy as np
from astropy.table import Table, Column, join
import matplotlib.pyplot as plt


import fermipy.utils as utils


def create_mask(tab,target_def):
    """Create a table mask from a target definition."""

    m = np.empty(len(tab),dtype=bool); m.fill(True)

    for k,v in target_def.items():

        print k, v, type(v)
        
        if isinstance(v,list):

            m0 = np.empty(len(tab0),dtype=bool)

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

            print m0
                
            m &= m0
            
    
    return m

    

if __name__ == '__main__':

    usage = "usage: %(prog)s"
    description = "Perform stacking analysis."
    parser = argparse.ArgumentParser(usage=usage,description=description)

    parser.add_argument('--output', default = 'stack.fits')
    parser.add_argument('--targets', default = None, required=True)

    parser.add_argument('tables', nargs='+', default = None,
                        help='Run analyses in all subdirectories of this '
                        'directory.')

    args = parser.parse_args()

    targets = yaml.load(open(args.targets))
    tab0 = Table.read(args.tables[0])
    tab1 = Table.read(args.tables[1])
    tab0 = join(tab0,tab1)

    targets_mask = {}
    
    for k,v in targets.items():
        targets_mask[k] = create_mask(tab0,v)
    

    cols_dict = OrderedDict()
    cols_dict['name'] = dict(dtype='S20', format='%s',description='Target Set')        
    cols_dict['fit_halo_scan_ts'] = dict(dtype='f8', format='%.2f',shape=(9,7))
    cols_dict['fit_halo_scan_dlnl'] = dict(dtype='f8', format='%.2f',shape=(9,7,31))
    cols_dict['fit_halo_scan_eflux_ul95'] = dict(dtype='f8', format='%.4g',shape=(9,7))

    cols_dict['fit_ext_ts'] = dict(dtype='f8', format='%.2f')
    cols_dict['fit_ext_dlnl'] = dict(dtype='f8', format='%.2f',shape=(26))

    
    cols = [Column(name=k, **v) for k,v in cols_dict.items()]

    tab_out = Table(cols)
    
    row_dict = {}
    for name, v in targets.items():
        
        m = targets_mask[name]
        tab = tab0[m]

        # These should eventually be extracted from the FITS file
        ext_width = np.logspace(-2.5,0,26)
        eflux = np.logspace(-8,-5,31)
        halo_ts = np.zeros((9,7))
        halo_eflux_ul95 = np.zeros((9,7))

        # Sum over all sources
        fit_halo_dlnl = np.sum(tab['fit_halo_scan_dlnl'],axis=0)
        fit_ext_dlnl =  np.sum(tab['fit_ext_scan_dlnl'],axis=0)

        fit_ext_dlnl -= fit_ext_dlnl[0]
        
        for i in range(9):
            for j in range(7):

                fit_halo_dlnl_sum = fit_halo_dlnl[i,j]
                fit_halo_dlnl_sum -= fit_halo_dlnl_sum[0]                
                lims = utils.get_parameter_limits(eflux,fit_halo_dlnl_sum)
                halo_ts[i,j] = 2.0*lims['lnlmax']
                halo_eflux_ul95[i,j] = lims['ul']

        row_dict['name'] = name
        row_dict['fit_halo_scan_dlnl'] = fit_halo_dlnl
        row_dict['fit_halo_scan_ts'] = halo_ts
        row_dict['fit_halo_scan_eflux_ul95'] = halo_eflux_ul95
        row_dict['fit_ext_dlnl'] = fit_ext_dlnl


        ext_lims = utils.get_parameter_limits(ext_width,fit_ext_dlnl)

        row_dict['fit_ext_ts'] = 2.0*(ext_lims['lnlmax'])        
        
        tab_out.add_row([row_dict[k] for k in cols_dict.keys()])


        xedge = np.logspace(-1,0,9)
        yedge = np.linspace(1.5,3.0,7)
        
        # Make diagnostics
        plt.figure()
        im = plt.pcolormesh(xedge,yedge,halo_ts.T)
        plt.colorbar(im)
        
        plt.gca().set_xscale('log')
        plt.gca().set_ylim(1.5,3.0)
        plt.gca().set_title(name)
        
#        im = plt.pcolormesh(xedge,yedge,halo_ts.T,origin='lower',interpolation='bicubic')
#        plt.colorbar(im)
#        plt.contour(halo_ts.T)
        

    tab_out.write(args.output,format='fits',overwrite=True)




        
    
    
    sys.exit(0)
        

    

    plt.figure()

    for i in range(len(tab)):
        plt.plot(eflux,fit_dlnl[i],color='grey')
    plt.plot(eflux,fit_dlnl_sum,color='k',linewidth=2)
    plt.gca().set_ylim(-5,5)
    plt.gca().set_xscale('log')

    plt.figure()

    for i in range(len(tab)):
        plt.plot(width,fit_ext_dlnl[i],color='grey')
    plt.plot(width,fit_ext_dlnl_sum,color='k',linewidth=2)
    plt.gca().set_ylim(-5,5)
    plt.gca().set_xscale('log')
