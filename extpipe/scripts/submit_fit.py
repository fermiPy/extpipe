from extpipe.batchfarm import utils,lsf
from extpipe.batchfarm.lsf import lsfDefaults
import extpipe
import yaml
import argparse
from astropy.table import Table, vstack
from os import path
import os
import logging
from glob import glob
from astropy.io import fits
import numpy as np

if __name__ == '__main__':
    usage = "usage: %(prog)s"
    description = "Fit SED for cascade component."
    parser = argparse.ArgumentParser(usage=usage,description=description)
    parser.add_argument('-c', default = False, required=True,
                        help='yaml config file')
    parser.add_argument('--combined', default = 1, type = int, 
                        help='perform a fit to combines TeV spectra of 1ES0229, 1ES1218, and H2356')
    args = parser.parse_args()
    utils.init_logging('INFO', color = True)
    config = yaml.load(open(args.c))
    config['configname'] = 'fit_igmf_th_jet{0[th_jet]:.0f}_tmax{0[tmax]:.0e}_lp_{0[kind]:s}_ts{0[tev_scale]}'.format(config)
    config['combined'] = args.combined

    tab_sed_tev = Table.read(config['cat_tev'])
    njobs = tab_sed_tev['SOURCE_FULL'].data.shape[0]
    if args.combined:
	jobid = []
	jobsrc = []
	cs = ['1ES0229+200', '1ES1218+304', 'H2356-309']
	for s in cs:
	    njobs -= np.sum(tab_sed_tev['SOURCE'] == s)
	    njobs += 1
	njobs = int(njobs)
	for i,sf in enumerate(tab_sed_tev['SOURCE']):
	    if sf in cs and (not sf in jobid):
		jobid.append(sf)
		jobsrc.append(tab_sed_tev['SOURCE_FULL'][i])
	    elif not sf in cs:
		jobid.append(sf)
		jobsrc.append(tab_sed_tev['SOURCE_FULL'][i])
	config['jobsrc'] = jobsrc
    else:
	config['jobsrc'] = list(tab_sed_tev['SOURCE_FULL'].data)

    if config['tev_scale'] == 1.:
	outfile = 'fit_igmf_th_jet{0[th_jet]:.0f}_tmax{0[tmax]:.0e}_lp_{0[kind]:s}_com{0[combined]}_0*.fits'.format(config)
    else:
	outfile = 'fit_igmf_th_jet{0[th_jet]:.0f}_tmax{0[tmax]:.0e}_lp_{0[kind]:s}_com{0[combined]}_ts{0[tev_scale]}_0*.fits'.format(config)
    outfile = path.join(config['outdir'],outfile)
    missing = utils.missing_files(outfile, njobs, num = 4, split = '.fits')

    if len(missing) < njobs:
	njobs = missing
	logging.info('there are {0:n} files missing in {1:s}'.format(len(missing),
	    outfile ))

    if len(missing):
	script = path.join(path.dirname(extpipe.__file__), 'scripts/run_fit.py')
	lsf.submit_lsf(script,
	    config,'',njobs, jname = 'th{0[th_jet]:.0f}t{0[tmax]:.0e}'.format(config),
	       logdir = path.join(config['outdir'],'log/'))
    else:
	logging.info("All files present.")
	files = glob(outfile)
	if config['tev_scale'] == 1.:
	    files = sorted(files, key = lambda f : int(path.basename(f).split('.')[0][-4:]))
	else:
	    files = sorted(files, key = lambda f : int(path.basename(f).split('.')[1][-4:]))
	t = []
	for f in files:
	    t.append(Table.read(f)) 
	    if config['combined']:
		t[-1].remove_column('loglike_comp')
	tab = vstack(t)
	hdulist = fits.HDUList()
	hdulist.append(fits.table_to_hdu(tab))
	hdulist[1].name = 'SCAN_DATA'    

	hdulist.writeto(outfile.split('0*')[0] + 'combined.fits', clobber=True)
	logging.info("Saved combined result table to {0:s}".format(outfile.split('0*')[0] + 'combined.fits'))
