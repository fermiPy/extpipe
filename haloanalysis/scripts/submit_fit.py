from haloanalysis.batchfarm import utils,lsf
from haloanalysis.batchfarm.lsf import lsfDefaults
import haloanalysis
import yaml
import argparse
from astropy.table import Table, vstack
from os import path
import os
import logging
from glob import glob
from astropy.io import fits

if __name__ == '__main__':
    usage = "usage: %(prog)s"
    description = "Fit SED for cascade component."
    parser = argparse.ArgumentParser(usage=usage,description=description)
    parser.add_argument('-c', default = False, required=True,
                        help='yaml config file')
    args = parser.parse_args()
    utils.init_logging('INFO', color = True)
    config = yaml.load(open(args.c))
    config['configname'] = 'fit_igmf_th_jet{0[th_jet]:.0f}_tmax{0[tmax]:.0e}_lp_{0[kind]:s}'.format(config)

    tab_sed_tev = Table.read(config['cat_tev'])
    njobs = tab_sed_tev['SOURCE_FULL'].data.shape[0]

    outfile = 'fit_igmf_th_jet{0[th_jet]:.0f}_tmax{0[tmax]:.0e}_lp_{0[kind]:s}_0*.fits'.format(config)
    outfile = path.join(config['outdir'],outfile)
    missing = utils.missing_files(outfile, njobs, num = 4, split = '.fits')

    if len(missing) < njobs:
	njobs = missing
	logging.info('there are {0:n} files missing in {1:s}'.format(len(missing),
	    outfile ))

    if len(missing):
	script = path.join(path.dirname(haloanalysis.__file__), 'scripts/run_fit.py')
	lsf.submit_lsf(script,
	    config,'',njobs, jname = 'fit', logdir = path.join(config['outdir'],'log/'))
    else:
	logging.info("All files present.")
	files = glob(outfile)
	files = sorted(files, key = lambda f : int(path.basename(f).split('.')[0][-4:]))
	t = []
	for f in files:
	    t.append(Table.read(f)) 
	tab = vstack(t)
	hdulist = fits.HDUList()
	hdulist.append(fits.table_to_hdu(tab))
	hdulist[1].name = 'SCAN_DATA'    

	hdulist.writeto(outfile.split('0*')[0] + 'combined.fits', clobber=True)
	logging.info("Saved combined result table to {0:s}".format(outfile.split('0*')[0] + 'combined.fits'))
