import os
import sys
import argparse
from collections import OrderedDict

import yaml

import numpy as np
from astropy.io import fits
from astropy.table import Table, Column, join, hstack
import matplotlib.pyplot as plt
import matplotlib

import fermipy.utils as utils
from fermipy.spectrum import *
from fermipy.castro import CastroData

from haloanalysis.utils import create_mask, load_source_rows
from haloanalysis.sed import HaloSED
from haloanalysis.model import CascModel, CascLike
from haloanalysis.model import scan_igmf_likelihood


def load_cache(tab_ext_filename, tab_lnl_filename, names,
               cachefile='cache.fits'):

    if os.path.isfile(cachefile):
        return Table.read(cachefile)

    tab0 = Table.read(tab_ext_filename)
    tab1 = Table.read(tab_lnl_filename)
    tab = join(tab0,tab1)
    tab = load_source_rows(tab, names)
    tab.write(cachefile)
    return tab


def main():

    usage = "usage: %(prog)s"
    description = "Fit SED for cascade component."
    parser = argparse.ArgumentParser(usage=usage,description=description)

    parser.add_argument('--output', default = 'igmf_casc_lnl.fits')
    parser.add_argument('--make_plots', default = False, action='store_true')
    parser.add_argument('--cache', default = False, action='store_true')
    parser.add_argument('--modelfile', default = False, required=True)
    parser.add_argument('--sedfile', default = False, required=True,
                        help='FITS file containing the TeV SEDs.')
    parser.add_argument('--nstep', default = 5)
    parser.add_argument('--name', default = '1es0229+200')
    parser.add_argument('tables', nargs='+', default = None,
                        help='Extension and likelihood tables.')

    args = parser.parse_args()

    # list of sources
    src_names = ['1es0229+200', '1ES0347-121']

    # Use cached fits file
    if args.cache:
        tab_casc = load_cache(args.tables[0], args.tables[1], src_names)
    else:
        tab0 = Table.read(args.table[0])
        tab1 = Table.read(args.table[1])
        tab_casc = join(tab0, tab1)
        tab_casc = load_source_rows(tab_casc, src_names)

    tab_tev = Table.read(args.sedfile)
        
    for row0 in tab_casc:
        row1 = load_source_rows(tab_tev, [row0['assoc']], key='SOURCE')        
        scan_igmf_likelihood(row, args.modelfile, args.sedfile, args.output, args.nstep)


if __name__ == "__main__":
    main()
