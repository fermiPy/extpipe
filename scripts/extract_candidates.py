from astropy.table import Table
import sys
import argparse
import numpy as np
import yaml
from fermipy.utils import tolist

exclude = ['3fgl_j0059.0-7242e', # SMC
           '3fgl_j1324.0-4330e', # Cen A Lobes
           '3fgl_j2051.0+3040e', # Cyg Loop
           '3fgl_j0526.6-6825e', # LMC
           '3fgl_j0524.5-6937',  # LMC
           '3fgl_j0525.2-6614',  # LMC
           '3fgl_j0456.2-6924',  # LMC
           '3fhl_j0537.9-6909',  # LMC
           '3fgl_j0534.5+2201',  # Crab Pulsar
           '3fgl_j0534.5+2201s', # Crab Sync
           '3fgl_j0004.2+6757',  # fhes_j0001.3+6841
           '3fgl_j0008.5+6853',  # fhes_j0001.3+6841
           '3fgl_j2356.9+6812',  # fhes_j0001.3+6841
           '3fgl_j2355.4+6939',  # fhes_j0001.3+6841
           '3fgl_j1334.3-4152',  # fhes_j1334.8-4124
           '3fgl_j1335.2-4056',  # fhes_j1334.8-4124
           '3fgl_j1626.2-2428c', # fhes_j1627.3-2430
           '3fgl_j1628.2-2431c', # fhes_j1627.3-2430
           '3fgl_j0431.7+3503',  # fhes_j0429.3+3529
           '3fgl_j0426.3+3510',  # fhes_j0429.3+3529
           '3fgl_j0429.8+3611c', # fhes_j0429.3+3529
           '3fgl_j1748.5-3912',  # fhes_j1748.5-3912e
           #'3fgl_j0127.9+2551',  # spurious source due to localization issues
           #'3fgl_j0226.7-4747',  # spurious source due to localization issues
           #'3fgl_j1948.1-7059',  # spurious source due to localization issues
           #'3fgl_j1729.0+6049',  # spurious source due to localization issues
           '3fgl_j2125.8+5832',  # fhes_j2125.8+5833e
           '3fgl_j2206.5+6451',  # fhes_j2208.4+6501
           '3fgl_j2210.2+6509',  # fhes_j2208.4+6501
           '3fgl_j1720.3-0428',  #
           #'fhes_j1820.4-3217e', # spurious source?
           #'fhes_j2125.8+5833',
           '3fgl_j2309.0+5428', # fhes_j2309.0+5429e
           'lmc_p1',
           'lmc_p2',
           'lmc_p2',
           'lmc_p4',
           '',
]

def main():

    usage = "usage: %(prog)s [config files]"
    description = "Merge FITS tables."
    parser = argparse.ArgumentParser(usage=usage,description=description)

    parser.add_argument('--output', default = None)
    parser.add_argument('--filter', default = None, type=str)
    parser.add_argument('--type', default = '', type=str)
    parser.add_argument('--ts_threshold', default = None, type=float)
    parser.add_argument('--glat', default = 0.0, type=float)   
    parser.add_argument('files', nargs='*', default = None,
                        help='One or more FITS files containing BINTABLEs.')

    args = parser.parse_args()

    tab = Table.read(sys.argv[1])

    tab.sort('name')

    if args.type == 'halo':
        m0 = (tab['fit_halo_ts'] > args.ts_threshold) & (tab['fit_halo_ts'] > tab['fit_ext_ts_ext'])
    elif args.type == 'ext':
        m0 = (tab['fit_ext_ts_ext'] > args.ts_threshold) & (tab['fit_halo_ts'] < tab['fit_ext_ts_ext'])
    else:
        m0 = (tab['fit_ext_ts_ext'] > args.ts_threshold) | (tab['fit_halo_ts'] > args.ts_threshold)
        
    m1 = np.abs(tab['glat'])> args.glat
    tab0 = tab[m0&m1]

    #tab0.reverse()

    row_str = '{name:20s} {name_3fgl:20s} {codename:25s} {assoc:25s} {glon:8.2f} {glat:8.2f} '
    row_str += '{gauss_ts_ext:8.1f} {disk_ts_ext:8.1f} {ts_ext:8.1f} {ext_r68:8.2f} {ext_r68_err:8.2f} {ts_halo:8.1f} '
    row_str += '{ts:8.1f} {ext_model:8s} {ext_index:8.2f} {idx_gauss:8d} {idx_disk:8d} {idx_halo:8d}'

    names = []
    
    print('%20s %20s %25s %25s %8s %8s %8s %8s %8s %8s %8s %8s %8s'%('name','name_3fgl','codename',
                                                                     'assoc','glon','glat',
                                                                     'ts_ext','ts_ext','ts_ext','r68','r68_err','ts_halo','ts'))

    for row in tab0:
        if row['codename'] in exclude:
            continue

        delta_ts = 2.0*(row['fit_ext_gauss_loglike']-row['fit_ext_disk_loglike'])-(row['fit_ext_gauss_ts']-row['fit_ext_disk_ts'])

        names += [row['codename']]
        
        #if row['fit_ext_gauss_ts'] > row['fit_ext_gauss_ts_ext']:
        #    continue

        #if np.abs(delta_ts) < 0.1:
        #    continue

        o = row_str.format(name=row['name'],
                           name_3fgl=row['name_3fgl'],
                           codename=row['codename'],assoc=row['assoc'],glon=row['glon'],glat=row['glat'],
                           gauss_ts_ext=row['fit_ext_gauss_ts_ext'],
                           disk_ts_ext=row['fit_ext_disk_ts_ext'],
                           ts_ext=row['fit_ext_ts_ext'],
                           ext_r68=row['fit_ext_r68'],
                           ext_r68_err=row['fit_ext_r68_err'],
                           ts_halo=row['fit_halo_ts'], ts=row['fitn_ts'][0],
                           idx_gauss=row['fit_idx_ext_gauss'],
                           idx_disk=row['fit_idx_ext_disk'],
                           idx_halo=row['fit_idx_halo'],
                           ext_model=row['fit_ext_model'],
                           ext_index=row['fit_ext_index'])
        print(o)
        continue

        print(row['fit_ext_gauss_ts_ext'],row['fit_ext_gauss_ts'], row['fit_ext_disk_ts_ext'],row['fit_ext_disk_ts'],)
        print(2.0*(row['fit_ext_gauss_loglike']-row['fit_ext_disk_loglike']),row['fit_ext_gauss_ts']-row['fit_ext_disk_ts'])
        print(row['fit_idx_ext_gauss'], row['fit_idx_ext_disk'])
        print row['fit_ext_gauss_loglike']
        print row['fit_ext_disk_loglike']

        print row['fit_dlike1_ext_gauss']
        print row['fit_dlike1_ext_disk']
        print row['fitn_dlike1']

        print('-'*80)

    if args.output is not None:
        yaml.dump(tolist(names),open(args.output,'w'))

if __name__ == "__main__":
    main()
