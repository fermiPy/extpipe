import os
import sys
import argparse
from collections import OrderedDict

import yaml

import numpy as np
from astropy.io import fits
from astropy.table import Table, Column, join, hstack, vstack
import matplotlib.pyplot as plt
import matplotlib

import fermipy.utils as utils
from fermipy.spectrum import *
from fermipy.castro import CastroData

from haloanalysis.utils import create_mask, load_source_rows
from haloanalysis.sed import HaloSED
from haloanalysis.model import CascModel, CascLike
from haloanalysis.model import scan_igmf_likelihood
from haloanalysis.model import LogParabolaExpCutoff
from fermipy.spectrum import PLExpCutoff



def load_cache(tables, names, cachefile='cache.fits'):

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
    parser.add_argument('--modelfile', default = False, required=True,
                        help='FITS file containing the IGMF models')
    parser.add_argument('--sedfile', default = False, required=True,
                        help='FITS file containing the TeV SEDs.')
    parser.add_argument('--nstep', default = 5, type=int)
    parser.add_argument('--casc_scale', default = 1.0, type=float)
    parser.add_argument('--casc_r68_scale', default = 1.0, type=float)
    parser.add_argument('--name', default = [], action='append')
    parser.add_argument('--key', default = 'SOURCE_FULL')
    parser.add_argument('--int_spec', default = 'PLExpCutOff',
    			help='intrinsic source spec. as in fermipy.spectrum',
			choices=['LogParabolaExpCutoff','PLExpCutoff'])
    parser.add_argument('tables', nargs='+', default = None,
                        help='Extension and likelihood tables.')

    args = parser.parse_args()

    # list of sources
    src_names = args.name
    
    casc_model = CascModel.create_from_fits(args.modelfile)
    casc_model.set_eblmodel(eblmodel = 'dominguez')

    tab_pars = Table.read(args.tables[0],hdu='SCAN_PARS')
    tab_ebounds = Table.read(args.tables[0],hdu='EBOUNDS')

    # intrinsic spec
    if args.int_spec == 'PLExpCutOff':
	fint = sp.PLExpCutoff  # intrinsic spectrum
	p0 = [1E-13,-1.5,1E7]  # initial parameters
	scale = 1E3            # pivot energy
    elif args.int_spec == 'LogParabolaExpCutoff':
	fint = LogParabolaExpCutoff  # intrinsic spectrum
	p0 = [1E-13,-1.5,0.,1E7]  # initial parameters
	scale = 1E3            # pivot energy
    
    # Use cached fits file
    if args.cache:
        tab_casc = load_cache(args.tables, src_names)
    else:
        tables = [Table.read(args.tables[0],'CATALOG'),
                  Table.read(args.tables[0],'LIKELIHOOD'),
                  Table.read(args.tables[0],'SED')]

        for i, t in enumerate(tables):
            if 'NAME' in t.columns:
                t['name'] = t['NAME']

        tab_casc = join(tables[0],tables[1])
        tab_casc = join(tab_casc,tables[2])
        #if src_names:
        #    tab_casc = load_source_rows(tab_casc, src_names)
        
    tab_sed_tev = Table.read(args.sedfile)

    tab_igmf = []

    for name in src_names:

        rows_sed_tev = load_source_rows(tab_sed_tev, [name], key=args.key)
        #cat_names = [ '3FGL %s'%row['3FGL_NAME'] for row in rows_sed_tev ]
        cat_names = [ row['3FGL_NAME'] for row in rows_sed_tev ]
        cat_names = np.unique(np.array(cat_names))
        #rows_sed_gev = load_source_rows(tab_casc, cat_names, key='name')
        #rows_casc = load_source_rows(tab_casc, cat_names, key='name')
        rows_sed_gev = load_source_rows(tab_casc, cat_names, key='name_3fgl')
        rows_casc = load_source_rows(tab_casc, cat_names, key='name_3fgl')
        #tab = scan_igmf_likelihood(casc_model, rows_sed_tev, rows_sed_gev,
                                   #rows_casc, tab_pars, tab_ebounds, args.nstep,
                                   #args.casc_scale, args.casc_r68_scale)
	tab = scan_igmf_likelihood(casc_model, rows_sed_tev, rows_sed_gev,
				    rows_casc, tab_pars, tab_ebounds, args.nstep,
				    casc_scale = args.casc_scale, 
				    casc_r68_scale = args.casc_r68_scale,
				    p0 = p0,
				    scale = scale,
				    fint = fint)
        tab_igmf += [tab]

    tab = vstack(tab_igmf)
        
    hdulist = fits.HDUList()
    hdulist.append(fits.table_to_hdu(tab))
    hdulist[1].name = 'SCAN_DATA'    
    hdulist.writeto(args.output, clobber=True)
        

if __name__ == "__main__":
    main()
