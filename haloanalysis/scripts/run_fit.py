import os
from os import path
import sys

import yaml
import logging

import numpy as np
from astropy.io import fits
from astropy.table import Table, join

from fermipy.spectrum import *

from haloanalysis.utils import create_mask, load_source_rows
from haloanalysis.model import CascModel, CascLike, LogParabolaExpCutoff
from haloanalysis.model import scan_igmf_likelihood

import argparse
from haloanalysis.batchfarm import utils,lsf

if __name__ == '__main__':
    usage = "usage: %(prog)s"
    description = "Fit SED for cascade component."
    parser = argparse.ArgumentParser(usage=usage,description=description)
    parser.add_argument('-c', default = False, required=True,
                        help='yaml config file')
    parser.add_argument('-i', required=False, default = 0, 
                        help='TeV SEDs id.', type=int)
    args = parser.parse_args()

    utils.init_logging('INFO')
    config = yaml.load(open(args.c))
    tmpdir, job_id = lsf.init_lsf(local_id = args.i)
    os.chdir(tmpdir)    # go to tmp directory
    logging.info('Entering directory {0:s}'.format(tmpdir))
    logging.info('PWD is {0:s}'.format(os.environ["PWD"]))

# --- get the extension catalog ------------------- #
    cat = utils.copy2scratch(config['cat'], tmpdir)
    tab_pars = Table.read(cat,hdu='SCAN_PARS')
    tab_ebounds = Table.read(cat,hdu='EBOUNDS')

    tables = [Table.read(cat,'CATALOG'),
	    Table.read(cat,'LIKELIHOOD'),
	    Table.read(cat,'SED')]

    tab_casc = join(tables[0],tables[1])
    tab_casc = join(tab_casc,tables[2])
# --- Load the TeV SEDs --------------------------- #
    cat_tev = utils.copy2scratch(config['cat_tev'], tmpdir)
    tab_sed_tev = Table.read(cat_tev)

    src_name = tab_sed_tev['SOURCE_FULL'][job_id - 1]

# --- Load the Cascade model ---------------------- #
    basedir = config['basedir']
    tmax_str = '{0:.0e}'.format(config['tmax']).replace('+','').replace('0','')
    if config['kind'] == 'elmag':
        filename = path.join(config['basedir'],
	                'th_jet{0[th_jet]:.2f}/gam-2.00/results_merged_z_th{0[th_jet]:.0f}d_t{1:s}.fits'.format(
			config, tmax_str))
    elif kind == 'analytic':
	filename = path.join(basedir,
			'th_jet{0[th_jet]:.2f}/gam-2.00/results_analytical_merged_z_th{0[th_jet]:.0f}d_t{1:s}.fits'.format(
			config, tmax_str))
	      
    filename = utils.copy2scratch(filename, tmpdir)
    casc_model = CascModel.create_from_fits(filename)
    casc_model.set_eblmodel(eblmodel = config['ebl_model'])
# --- Perform the analysis ------------------------- #
    # intrinsic spec
    if config['int_spec'] == 'PLExpCutOff':
	fint = sp.PLExpCutoff  # intrinsic spectrum
	p0 = [1E-13,-1.5,1E7]  # initial parameters
	scale = 1E3            # pivot energy
    elif config['int_spec'] == 'LogParabolaExpCutoff':
	fint = LogParabolaExpCutoff  # intrinsic spectrum
	p0 = [1E-13,-1.5,0.,1E7]  # initial parameters
	scale = 1E3            # pivot energy

    rows_sed_tev = load_source_rows(tab_sed_tev, [src_name], key='SOURCE_FULL')
    #cat_names = [ '3FGL %s'%row['3FGL_NAME'] for row in rows_sed_tev ]
    cat_names = [ row['3FGL_NAME'] for row in rows_sed_tev ]
		    
    cat_names = np.unique(np.array(cat_names))
    rows_sed_gev = load_source_rows(tab_casc, cat_names, key='name_3fgl')
    rows_casc = load_source_rows(tab_casc, cat_names, key='name_3fgl')
						        
    logging.info("Starting analysis for {0:s}".format(src_name))
    tab = scan_igmf_likelihood(casc_model, rows_sed_tev, rows_sed_gev,
		    rows_casc, tab_pars, tab_ebounds, config['nstep'],
		    casc_scale = config['casc_scale'], 
		    casc_r68_scale = config['casc_r68_scale'],
		    p0 = p0,
		    scale = scale,
		    fint = fint)

    hdulist = fits.HDUList()
    hdulist.append(fits.table_to_hdu(tab))
    hdulist[1].name = 'SCAN_DATA'    
    outfile = 'fit_igmf_th_jet{0[th_jet]:.0f}_tmax{0[tmax]:.0e}_lp_{0[kind]:s}_{1:04n}.fits'.format(config,job_id)
    hdulist.writeto(path.join(tmpdir,outfile), clobber=True)
    logging.info("Saved output to {0:s}".format(outfile))

    utils.copy2scratch(outfile,config['outdir'])

