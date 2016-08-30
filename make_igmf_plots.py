import matplotlib.pyplot as plt
from astropy.table import Table, join
import numpy as np
from haloanalysis.utils import load_source_rows
from haloanalysis.sed import SED
import matplotlib
import sys
import os


for f in sys.argv[1:]:
    tab = Table.read(f)


    prefix = os.path.splitext(f)[0]

    igmf = tab['igmf'][0]*np.ones(tab['dloglike'][0].shape)
    lcoh = tab['lcoh'][0]*np.ones(tab['dloglike'][0].shape)
    dloglike = tab['dloglike'][0]

    tab_ebounds = Table.read('tables/table_std_all_tev_sources_lnl.fits','EBOUNDS')
    tab_sed_tev = Table.read('../haloanalysis/data/CompiledTeVSources.fits')
    tab_casc = join(Table.read('tables/table_std_all_tev_sources_lnl.fits'), Table.read('tables/table_std_all_tev_sources.fits'))
    tab_casc_sed = Table.read('tables/table_std_all_tev_sources_sed.fits')

    source_full = tab['SOURCE_FULL'][0]

    rows_sed_tev = load_source_rows(tab_sed_tev, [tab['SOURCE_FULL'][0]], key='SOURCE_FULL')
    rows_sed_gev = load_source_rows(tab_casc_sed, [tab['name'][0]], key='NAME')

    sed_prim0 = SED.create_from_row(rows_sed_tev)
    sed_prim1 = SED.create_from_row2(rows_sed_gev, tab_ebounds)

    plt.figure()

    sed_prim0.plot(color='k',marker='o')
    sed_prim1.plot(color='k',marker='s')
    #plt.errorbar(sed_prim.ectr/1E6,
    #             sed_prim.ectr*sed_prim.flux,
    #             sed_prim.ectr*sed_prim.flux_err,
    #             marker='o',linestyle='None')


    lcoh_idx = 10

    m = np.isfinite(tab['prim_tev_flux'][0][lcoh_idx,lcoh_idx])




    idxs = np.arange(0,20,3)

    for i, idx in enumerate(idxs):
        plt.plot(sed_prim0.ectr/1E6, sed_prim0.ectr*tab['prim_tev_flux'][0][lcoh_idx,idx][m],
                 label='log10(B/G) = %.2f'%(igmf[lcoh_idx,idx]),color=matplotlib.cm.jet(float(i)/float(len(idxs))))

        plt.plot(sed_prim1.ectr/1E6, sed_prim1.ectr*tab['prim_flux'][0][lcoh_idx,idx],
                 #label='B = %.2f'%(igmf[lcoh_idx,idx]),
                 color=matplotlib.cm.jet(float(i)/float(len(idxs))))

        plt.plot(sed_prim1.ectr/1E6, sed_prim1.ectr*tab['casc_flux'][0][lcoh_idx,idx],
                 #label='B = %.2f'%(igmf[lcoh_idx,idx]),
                 color=matplotlib.cm.jet(float(i)/float(len(idxs))),linestyle='--')


        
    plt.gca().set_xscale('log')
    plt.gca().set_yscale('log')
    plt.gca().set_xlabel('Energy [TeV]')
    plt.gca().set_ylabel('Energy Flux [MeV cm$^{-2}$ $^{-1}$]')

    plt.legend(frameon=False,ncol=2)
    plt.gca().set_ylim(1E-9,1E-4)
    plt.gca().set_xlim(1E-3,30.)

    plt.savefig('%s_igmf_flux.png'%prefix)

    plt.figure()

    for i, idx in enumerate(idxs):
        plt.plot(sed_prim1.ectr/1E6, tab['casc_r68'][0][lcoh_idx,idx],
                 label='B = %.2f'%(igmf[lcoh_idx,idx]),
                 color=matplotlib.cm.jet(float(i)/float(len(idxs))))

    plt.plot(sed_prim1.ectr/1E6,np.sqrt((2.85*(sed_prim1.ectr/100.)**-0.8)**2+0.035**2),color='k',label='SOURCE::PSF3')

    #print np.sqrt((2.85*(sed_prim1.ectr/100.)**-0.8)**2+0.035**2)
    
    plt.legend(frameon=False,ncol=2)
    plt.gca().set_xscale('log')
    plt.gca().set_yscale('log')

    plt.gca().set_xlabel('Energy [TeV]')
    plt.gca().set_ylabel('R68 [deg]')

    plt.gca().set_ylim(1E-3,10.)
    plt.gca().set_xlim(1E-3,30.)


    
    plt.savefig('%s_igmf_r68.png'%prefix)
    #plt.figure()
    #plt.imshow(dloglike,
    #           interpolation='nearest',
    #           extent=[np.min(tab['igmf'][0]),np.max(tab['igmf'][0]),
    #                   np.min(tab['lcoh'][0]),np.max(tab['lcoh'][0])],
    #           origin='lower')



    #plt.contour(igmf,lcoh,dloglike.T,            
    #            levels=[-2.0],color='k')

    dloglike -= dloglike[-1,-1]
    
    plt.figure()
    plt.pcolormesh(lcoh,igmf,dloglike,cmap='BuGn')
    cb = plt.colorbar()
    cb.set_label('Delta-LogLikelihood')
    plt.contour(lcoh,igmf,dloglike,
                levels=[-3.0,3.0],colors=['k'])

    plt.gca().set_xlabel('log10(Lcoh/Mpc)')
    plt.gca().set_ylabel('log10(B/G)')
    plt.savefig('%s_igmf_constraints.png'%prefix)

    plt.close('all')
