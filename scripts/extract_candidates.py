from astropy.table import Table
import sys
import argparse
import numpy as np
import yaml
from fermipy.utils import tolist

exclude_dict = {
    'known_sources' : [
        '3fgl_j0059.0-7242e', # SMC
        '3fgl_j1324.0-4330e', # Cen A Lobes
        '3fgl_j2051.0+3040e', # Cyg Loop
        '3fgl_j0526.6-6825e', # LMC
        '3fgl_j0524.5-6937',  # LMC
        '3fgl_j0525.2-6614',  # LMC
        '3fgl_j0456.2-6924',  # LMC
        '3fhl_j0537.9-6909',  # LMC
        'lmc_p1',
        'lmc_p2',
        'lmc_p2',
        'lmc_p4',
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
        ]
    }

def dump_source(tab, name):

    row = tab[tab['codename'] == name][0]


    print 'codename:  {:20s}'.format(row['codename'])
    print 'name:      {:20s}'.format(row['name'])
    print 'model:     {:20s}'.format(row['fit_ext_model'])
    print 'ts:        {:8.1f}'.format(row['fit_ext_ts'])
    print 'ts_ext:    {:8.1f}'.format(row['fit_ext_ts_ext'])
    print 'ext_flux:  {:8.4g}'.format(row['fit_ext_flux'])
    print 'ext_index: {:8.2f} +/- {:8.2f}'.format(np.abs(row['fit_ext_index']),np.abs(row['fit_ext_index_err']))
    print 'ext_r68:   {:8.2f} +/- {:8.2f}'.format(row['fit_ext_r68'],row['fit_ext_r68_err'])
    print 'ext_idx: {:4d} gauss_idx: {:4d} disk_idx: {:4d} halo_idx: {:4d}'.format(row['fit_idx_ext'],
                                                                                   row['fit_idx_ext_gauss'],
                                                                                   row['fit_idx_ext_disk'],
                                                                                   row['fit_idx_halo'])
    
    row_str = '{iter:>4s} {loglike:>10s} {aic_ext:>10s} {aic_ps:>10s} {glon:>8s} {glat:>8s} {ts:>8s} {ts_ext:>8s} {r68:>8s} {r68_err:>8s} {index:>8s} {daic:>8s}'
    header = row_str.format(iter='iter', loglike='loglike', glon='glon', glat='glat', ts='ts',
                            r68='r68', r68_err='r68_err', ts_ext='ts_ext', daic='daic', index='index',aic_ext='aic_ext',aic_ps='aic_ps')

    print '-'*80
    print header
    
    for i in range(5):
        print row_str.format(
            iter='{:d}'.format(i),
            loglike='{:10.1f}'.format(row['fitn_ext_gauss_loglike'][i]),
            aic_ext='{:10.1f}'.format(row['fitn_aic_ext_gauss'][i]),
            aic_ps='{:10.1f}'.format(row['fitn_aic_ps'][i]),
            glon='{:8.2f}'.format(row['fitn_ext_gauss_glon'][i]),
            glat='{:8.2f}'.format(row['fitn_ext_gauss_glat'][i]),
            ts='{:8.1f}'.format(row['fitn_ext_gauss_ts'][i]),
            r68='{:8.2f}'.format(row['fitn_ext_gauss_r68'][i]),
            r68_err='{:8.2f}'.format(row['fitn_ext_gauss_r68_err'][i]),
            index='{:8.2f}'.format(row['fitn_ext_gauss_index'][i]),
            daic='{:8.1f}'.format(row['fitn_daic_ps_ext_gauss'][i]),
            ts_ext='{:8.1f}'.format(row['fitn_ext_gauss_ts_ext'][i]))

    print '-'*80
    print header
    
    for i in range(5):
        print row_str.format(
            iter='{:d}'.format(i),
            loglike='{:10.1f}'.format(row['fitn_ext_disk_loglike'][i]),
            aic_ext='{:10.1f}'.format(row['fitn_aic_ext_disk'][i]),
            aic_ps='{:10.1f}'.format(row['fitn_aic_ps'][i]),
            glon='{:8.2f}'.format(row['fitn_ext_disk_glon'][i]),
            glat='{:8.2f}'.format(row['fitn_ext_disk_glat'][i]),
            ts='{:8.1f}'.format(row['fitn_ext_disk_ts'][i]),
            r68='{:8.2f}'.format(row['fitn_ext_disk_r68'][i]),
            r68_err='{:8.2f}'.format(row['fitn_ext_disk_r68_err'][i]),
            index='{:8.2f}'.format(row['fitn_ext_disk_index'][i]),
            daic='{:8.1f}'.format(row['fitn_daic_ps_ext_disk'][i]),
            ts_ext='{:8.1f}'.format(row['fitn_ext_disk_ts_ext'][i]))

    print '-'*80
    print header
    
    for i in range(5):
        print row_str.format(
            iter='{:d}'.format(i),
            loglike='{:10.1f}'.format(row['fitn_dlike1_halo'][i]),
            aic_ext='---',
            aic_ps='---',
            glon='---',
            glat='---',
            ts='{:8.1f}'.format(row['fitn_halo_ts'][i]),
            r68='---',
            r68_err='---',
            index='---',
            daic='{:8.1f}'.format(row['fitn_daic_ps_halo'][i]),
            ts_ext='---')
        
    

def main():

    usage = "usage: %(prog)s [config files]"
    description = "Merge FITS tables."
    parser = argparse.ArgumentParser(usage=usage,description=description)

    parser.add_argument('--output', default = None)
    parser.add_argument('--sort', default = 'name')
    parser.add_argument('--filter', default = None, type=str)
    parser.add_argument('--exclude', default = None, type=str)
    parser.add_argument('--type', default = '', type=str)
    parser.add_argument('--target', default = None, type=str)
    parser.add_argument('--ts_threshold', default = None, type=float)
    parser.add_argument('--glat', default = 0.0, type=float)   
    parser.add_argument('files', nargs='*', default = None,
                        help='One or more FITS files containing BINTABLEs.')

    args = parser.parse_args()

    tab = Table.read(sys.argv[1])

    if args.target:
        dump_source(tab, args.target)
        sys.exit(0)
    
    tab.sort(args.sort)

    if args.type == 'halo':
        m0 = (tab['fit_halo_ts'] > args.ts_threshold) & (tab['fit_halo_ts'] > tab['fit_ext_ts_ext'])
    elif args.type == 'ext':
        m0 = (tab['fit_ext_ts_ext'] > args.ts_threshold) #& (tab['fit_halo_ts'] < tab['fit_ext_ts_ext'])
    elif args.type == 'all':
        m0 = (tab['fit_ext_gauss_ts_ext'] > args.ts_threshold)
        m0 |= (tab['fit_ext_disk_ts_ext'] > args.ts_threshold)
        m0 |= (tab['fit_halo_ts'] > args.ts_threshold)
    else:
        m0 = (tab['fit_ext_ts_ext'] > args.ts_threshold) | (tab['fit_halo_ts'] > args.ts_threshold)
        
    m1 = np.abs(tab['glat'])> args.glat
    tab0 = tab[m0&m1]

    print('Found ', len(tab0), 'sources')
    
    #tab0.reverse()

    row_str = '{codename:25s} {name_3fgl:20s} {assoc:25s} {glon:>8s} {glat:>8s} '
    row_str += '{gauss_ts_ext:>8s} {disk_ts_ext:>8s} {ts_ext:>8s} {ext_r68:>8s} '
    row_str += '{ext_r68_err:>8s} {ts_halo:>8s} '
    row_str += '{ts:>8s} {ext_model:8s} {ext_index:8s} '
    row_str += '{idx_gauss:4s} {idx_disk:4s} {idx_halo:4s} {dlnl_ext:8s} {dlnl_halo:8s}'

    names = []

    print(row_str.format(codename='codename',
                         name_3fgl='name_3fgl',
                         assoc='assoc',
                         glon='glon',
                         glat='glat',
                         gauss_ts_ext='ts_ext',
                         disk_ts_ext='ts_ext',
                         ts_ext='ts_ext',
                         ext_r68='r68',
                         ext_r68_err='r68_err',
                         ts_halo='ts_halo',
                         ts='ts',
                         ext_model='model',
                         ext_index='index',
                         idx_gauss='i_g',
                         idx_disk='i_d',
                         idx_halo='i_h',
                         dlnl_ext='aic_ext',
                         dlnl_halo='aic_halo',
                         ))
    
#    print('%25s %20s %25s %8s %8s %8s %8s %8s %8s %8s %8s %8s'%('codename','name_3fgl',
#                                                                'assoc','glon','glat',
#                                                                'ts_ext','ts_ext','ts_ext','r68','r68_err','ts_halo','ts'))

    if args.exclude is None:
        exclude = []
    else:
        exclude = exclude_dict.get(args.exclude)
    
    
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
                           codename=row['codename'],assoc=row['assoc'],
                           glon='{:8.2f}'.format(row['glon']),
                           glat='{:8.2f}'.format(row['glat']),
                           gauss_ts_ext='{:8.1f}'.format(row['fit_ext_gauss_ts_ext']),
                           disk_ts_ext='{:8.1f}'.format(row['fit_ext_disk_ts_ext']),
                           ts_ext='{:8.1f}'.format(row['fit_ext_ts_ext']),
                           ext_r68='{:8.2f}'.format(row['fit_ext_r68']),
                           ext_r68_err='{:8.2f}'.format(row['fit_ext_r68_err']),
                           ts_halo='{:8.1f}'.format(row['fit_halo_ts']),
                           ts='{:8.1f}'.format(row['fitn_ts'][0]),
                           idx_gauss='{:4d}'.format(row['fit_idx_ext_gauss']),
                           idx_disk='{:4d}'.format(row['fit_idx_ext_disk']),
                           idx_halo='{:4d}'.format(row['fit_idx_halo']),
                           ext_model=row['fit_ext_model'],
                           ext_index='{:8.2f}'.format(row['fit_ext_index']),
                           dlnl_ext='{:8.2f}'.format(row['fit_daic_ps_ext']),
                           dlnl_halo='{:8.2f}'.format(row['fit_daic_ps_halo']))
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
