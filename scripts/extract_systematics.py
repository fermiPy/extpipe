

from astropy.table import Table, Column
from astropy.coordinates import SkyCoord
import numpy as np
import sys
import yaml

def main():

    ext_sources = yaml.load(open(sys.argv[1]))
    srcnames = sorted(ext_sources.keys())

    #srcnames = ['3fgl_j1209.1-5224','fhes_j2125.8+5833e']
    
    tables = [Table.read('std_psf0123_joint2a_stdmodel/table_std_psf0123_joint2a_stdmodel_cat.fits'),
              Table.read('std_psf0123_joint2a_altmodel/table_std_psf0123_joint2a_altmodel_cat.fits')]

    for i in range(8):
        tables += [Table.read('std_psf0123_joint2a_model%i/table_std_psf0123_joint2a_model%i_cat.fits'%(i,i))]


    colnames = ['fit_ext_ts_ext','fit_ext_ts','fit_ext_r68','fit_ext_flux', 'fit_ext_ra', 'fit_ext_dec',
                'fit_ext_index', 'fit_ext_glon_err', 'fit_ext_glat_err']
    ext_colnames = ['fit_idx_ext_gauss', 'fit_idx_ext_disk', 'fit_idx_ext', 'fit_ext_model']

    cols_out = {'name' : []}

    import pprint

    for codename in srcnames:

        vals = {}
        ext_vals = {}

        for c in colnames + ['sep']:
            cols_out.setdefault(c + '_val',[])
            cols_out.setdefault(c + '_min',[])
            cols_out.setdefault(c + '_max',[])
            cols_out.setdefault(c + '_med',[])
            cols_out.setdefault(c + '_err_max',[])
            cols_out.setdefault(c + '_err_std',[])

        rows = [t[t['codename']==codename][0] for t in tables if t[t['codename']==codename]]

        flux = rows[0]['fit_ext_flux']
        flux_err = rows[0]['fit_ext_flux_err']

        print '-'*80
        print codename, rows[0]['name'], flux, flux_err

        cols_out['name'] += [rows[0]['name']]

        for c in ext_colnames:
            ext_vals[c] = np.array([r[c] for r in rows])

        for c in colnames:
            vals[c] = np.array([r[c] for r in rows])
        #    if 'flux' in c:
        #        vals['rel_' + c] /= vals[c][0]

        skydir = SkyCoord(vals['fit_ext_ra'], vals['fit_ext_dec'], unit='deg')

        vals['sep'] = skydir.separation(skydir[0]).deg

        print(vals['sep'])
        print(vals['fit_ext_ts_ext'])
        #print(ext_vals['fit_idx_ext_gauss'])
        print(ext_vals['fit_idx_ext'])
        print(['G'  if m == 'gauss' else 'D' for m in ext_vals['fit_ext_model']])

        for i, r in enumerate(rows):
            continue
            print(i, r['fitn_dlike1_ext_disk'])
            print(i, r['fitn_daic_ps_ext_disk'])
            print(i, r['fitn_dlike1_ext_gauss'])
            print(i, r['fitn_daic_ps_ext_gauss'])

        for c in colnames + ['sep']:

            v = vals[c]

            #if np.sum(np.isfinite(v)):
            #    m = np.isfinite(v)
            #else:
            #    m = ~np.isfinite(v)

            #print(c, len(v), np.sum(np.isfinite(v)))
            #print(v)
            #continue
    #        print(c)
    #        print(v)
    #        print(v[m])

            wts = np.ones(len(v)-1)
            sum_wts = np.sum(wts)

            cols_out[c + '_val'] += [v[0]]
            cols_out[c + '_min'] += [np.min(v)]
            cols_out[c + '_max'] += [np.max(v)]
            cols_out[c + '_med'] += [np.median(v)]
            
            if len(v) > 1:

                err_max = 0.5*np.max(np.abs(v[1:]-v[0]))
                err_std = np.sqrt(np.sum(wts*(v[1:]-v[0])**2)/sum_wts)
                if c == 'sep':
                    err_std /= np.sqrt(2.)
                
                cols_out[c + '_err_max'] += [err_max]
                cols_out[c + '_err_std'] += [err_std]
            else:
                cols_out[c + '_err_max'] += [np.nan]
                cols_out[c + '_err_std'] += [np.nan]

            print('%30s %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f'%(c,v[0],np.min(v),np.max(v),np.median(v),
                                                              np.std(v), 0.5*np.max(np.abs(v-v[0]))))

    #    pprint.pprint(vals)

    cols = [Column(name='codename',data=srcnames)]

    for k,v in cols_out.items():
        cols += [Column(name=k,data=v)]
    
    tab = Table(cols)
    tab.write('syst_table.fits',overwrite=True)
        
if __name__ == "__main__":
    main()
