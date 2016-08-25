import sys
from astropy.table import Table, vstack
import astropy.io.fits as pyfits
import argparse

def main():

    usage = "usage: %(prog)s [config files]"
    description = "Merge tables."
    parser = argparse.ArgumentParser(usage=usage,description=description)

    parser.add_argument('--output', default = None, required=True)   
    parser.add_argument('files', nargs='*', default = None,
                        help='One or more FITS files containing BINTABLEs.')

    args = parser.parse_args()
    
    h0 = pyfits.open(args.files[0])

    tables = []

    for f in sorted(args.files):

        print f

        tables += [Table.read(f)]


    tab = vstack(tables)
    tab.write(args.output,format='fits',overwrite=True)

    h = pyfits.open(args.output)
    h0[1] = h[1]
    h0.writeto(args.output,clobber=True)

if __name__ == "__main__":
    main()
