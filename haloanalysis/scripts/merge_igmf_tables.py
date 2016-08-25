import glob
from astropy.table import Table, Column
from astropy.io import fits

files = glob.glob('/nfs/farm/g/glast/u/mmeyer/projects/FermiHalo/Output/EBLm6/th_jet10.00/gam-2.00/z*/test/results_tmax1e+07.fits')




tables = []
for f in files:
    print 'loading ', f
    tables += [[Table.read(f),
                Table.read(f,'ENERGIES'),
                Table.read(f,'THETA'),
                Table.read(f,'FIXED_MODEL_PARS')]]

    
cols = []
for cname in tables[0][0].columns:

    col = tables[0][0][cname]
    shape = col.shape
    cols += [Column(name=col.name, dtype=col.dtype, shape=shape,
                    length=len(tables))]


cols += [Column(name='z', dtype='f8', length=len(tables), shape=(81,))]
    
tab = Table(cols)

for i, t in enumerate(tables):

    nebin = len(t[1]['E_cen'])

    #tab[]
    
    print i, nebin

    tab['z'][i,:] = t[3]['z']
    
    row = []
    for cname in tab.columns:

        if cname not in t[0].columns:
            continue
        
        if cname in ['inj_flux','prim_flux']:
            tab[cname][i,:,:nebin] = t[0][cname][:,:nebin]
        elif cname in ['casc_flux']:
            tab[cname][i,:,:nebin,:nebin] = t[0][cname][:,:nebin,:nebin]
        else:
            tab[cname][i,...] = t[0][cname][...]
            
        #row += [t[cname]]
        #tab[cname][i] = t[cname]
    #tab.add_row(row)

hdulist = fits.HDUList()
hdulist.append(fits.table_to_hdu(tab))
hdulist.append(fits.table_to_hdu(tables[0][1]))
hdulist.append(fits.table_to_hdu(tables[0][2]))
hdulist[1].name = 'FLUX_AND_MODEL_PARS'
hdulist[2].name = 'ENERGIES'
hdulist[3].name = 'THETA'    
hdulist.writeto('out.fits',clobber=True)

