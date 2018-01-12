import sys
import yaml
import numpy as np
from astropy.table import Table, vstack
from astropy.io import fits
import argparse
from haloanalysis.utils import create_mask

def main():

    usage = "usage: %(prog)s [config files]"
    description = "Merge FITS tables."
    parser = argparse.ArgumentParser(usage=usage,description=description)

    parser.add_argument('--output', default = None, required=True)
    parser.add_argument('--filter', default = None, type=str)
    parser.add_argument('--ts_threshold', default = None, type=float)   
    parser.add_argument('files', nargs='*', default = None,
                        help='One or more FITS files containing BINTABLEs.')

    args = parser.parse_args()
    h0 = fits.open(args.files[0])

    hdu_names = [t.name for t in h0]
    
    tables = {}

    for f in sorted(args.files):

        for t in hdu_names[1:]:
            
            tables.setdefault(t,[])
            tables[t] += [Table.read(f,t)]
        
    for k, v in tables.items():

        if k not in ['CATALOG','SED','LIKELIHOOD']:
            tables[k] = v[0]
        else:
            tables[k] = vstack(v)
            
    cat = tables['CATALOG']
    m = np.ones(len(cat),dtype=bool)
    
    if args.filter:
        cuts = yaml.load(open(args.filter))

        if 'inclusive' in cuts:
            m &= create_mask(tables['CATALOG'],cuts['inclusive'])

        if 'exclusive' in cuts:
            m &= ~create_mask(tables['CATALOG'],cuts['exclusive']) 

    if args.ts_threshold is not None:
        m &= (cat['ts'] > args.ts_threshold)
    
    tab_hdus = []
    for t in hdu_names[1:]:

        tab = tables[t]
        if t in ['CATALOG','SED','LIKELIHOOD']:
            tab = tab[m]
        
        tab_hdus += [fits.table_to_hdu(tab)]
            
    hdulist = fits.HDUList([fits.PrimaryHDU()] + tab_hdus)
    hdulist.writeto(args.output,clobber=True)
    

if __name__ == "__main__":
    main()
