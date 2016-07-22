import os
import sys
import argparse
from collections import OrderedDict

import yaml

import numpy as np
from astropy.io import fits
from astropy.table import Table, Column, join
import matplotlib.pyplot as plt
import matplotlib

import fermipy.utils as utils
from fermipy.spectrum import *
from fermipy.castro import CastroData

from haloanalysis.utils import create_mask
from haloanalysis.sed import HaloSED
from haloanalysis.model import CascModel, CascLike


def load_source_table(tab_ext_filename, tab_lnl_filename, name = '1ES 0229+200'):

    tab_ext = Table.read(tab_ext_filename)
    tab_lnl = Table.read(tab_lnl_filename)
    tab = join(tab_ext, tab_lnl)
    mask = create_mask(tab, {'assoc': [name]})
    return tab[mask]


def load_cache(tab_ext_filename, tab_lnl_filename, cachefile='cache.fits'):

    if os.path.isfile(cachefile):
        return Table.read(cachefile)

    tab = load_source_table(tab_ext_filename, tab_lnl_filename)
    tab.write(cachefile)

    return tab


def scan_igmf_likelihood(tab, modelfile, sedfile, outputfile, nstep):
    """This function uses the primary and cascade SEDs to scan the
    likelihood space in the IGMF parameter (B,L_coh).

    """
    
    fn = PLExpCutoff([1E-11,-1.5,1E6],scale=1E3)

    lcoh_scan = np.linspace(-4,4,nstep)
    igmf_scan = np.linspace(-20,-12,nstep)    
    bpars = np.meshgrid(lcoh_scan, igmf_scan)

    sed_prim = CastroData.create_from_flux_points(sedfile)
    sed_casc = HaloSED.create_from_fits(tab[0])
    hmm = CascModel.create_from_fits(modelfile)
    hl = CascLike(hmm, fn, sed_casc, sed_prim)

    model_lnl = np.zeros(bpars[0].shape)*np.nan
    p1 = fn.params

    for idx, x in np.ndenumerate(bpars[0]):

        p0 = [bpars[0][idx], bpars[1][idx]]
        lnl, p1 = hl.fit(p0,p1,method='SLSQP')
        model_lnl[idx] = lnl
        print(idx, lnl, p0, p1)
        #print bpars, lnl, p1

    cols_dict = OrderedDict()
    cols_dict['name'] = dict(dtype='S32', format='%s', description='name')
    cols_dict['assoc'] = dict(dtype='S32', format='%s', description='assoc')
    cols_dict['igmf_scan_dlnl'] = dict(dtype='f8', format='%.3f',
                                       shape=model_lnl.shape)
    
    tab_scan = Table([Column(name=k, **v) for k, v in cols_dict.items()])
    row_dict = {}
    row_dict['name'] = tab['name']
    row_dict['assoc'] = tab['assoc']
    row_dict['igmf_scan_dlnl'] = model_lnl
    tab_scan.add_row([row_dict[k] for k in cols_dict.keys()])

    cols_dict = OrderedDict()
    cols_dict['lcoh'] = dict(dtype='f8', format='%.3f', shape=bpars[0].shape,
                             data=bpars[0])
    cols_dict['igmf'] = dict(dtype='f8', format='%.3f', shape=bpars[1].shape,
                             data=bpars[1])
    tab_scan_grid = Table([Column(name=k, **v) for k, v in cols_dict.items()])
    
    hdulist = fits.HDUList()
    hdulist.append(fits.table_to_hdu(tab_scan))
    hdulist.append(fits.table_to_hdu(tab_scan_grid))

    hdulist[1].name = 'SCAN_DATA'
    hdulist[2].name = 'SCAN_GRID'
    
    hdulist.writeto(outputfile,clobber=True)

def main():

    usage = "usage: %(prog)s"
    description = "Fit SED for cascade component."
    parser = argparse.ArgumentParser(usage=usage,description=description)

    parser.add_argument('--output', default = 'igmf_casc_lnl.fits')
    parser.add_argument('--make_plots', default = False, action='store_true')
    parser.add_argument('--cache', default = False, action='store_true')
    parser.add_argument('--modelfile', default = False, required=True)
    parser.add_argument('--sedfile', default = False, required=True)
    parser.add_argument('--nstep', default = 5)
    parser.add_argument('tables', nargs='+', default = None,
                        help='Extension and likelihood tables.')

    args = parser.parse_args()

    # Use cached fits file
    if args.cache:
        tab = load_cache(args.tables[0], args.tables[1])
    else:
        tab = load_source_table(args.tables[0], args.tables[1])

    scan_igmf_likelihood(tab, args.modelfile, args.sedfile, args.output, args.nstep)


if __name__ == "__main__":
    main()
