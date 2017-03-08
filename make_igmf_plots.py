import matplotlib.pyplot as plt
from astropy.table import Table, join
import numpy as np
from haloanalysis.utils import load_source_rows
from haloanalysis.sed import SED
from matplotlib.offsetbox import AnchoredText
import matplotlib
import sys
import os

def preliminary(ax=None):

    if ax is None:
        ax = plt.gca()
    
    title = AnchoredText('Preliminary', loc=3,
                         prop={'size':18,'color': 'red'},
                         pad=0., borderpad=0.75,
                         frameon=False)
    ax.add_artist(title)

table_prefix = 'std_psf0123_joint2a'


for f in sys.argv[1:]:
    tab = Table.read(f)


    prefix = os.path.splitext(f)[0]

    igmf = tab['igmf'][0]*np.ones(tab['dloglike'][0].shape)
    lcoh = tab['lcoh'][0]*np.ones(tab['dloglike'][0].shape)
    dloglike = tab['dloglike'][0]

    tab_ebounds = Table.read('tables/table_%s_tev_sources_lnl.fits'%table_prefix,'EBOUNDS')
    tab_sed_tev = Table.read('../haloanalysis/data/CompiledTeVSources.fits')
    tab_casc = join(Table.read('tables/table_%s_tev_sources_lnl.fits'%table_prefix), Table.read('tables/table_%s_tev_sources.fits'%table_prefix))
    tab_casc_sed = Table.read('tables/table_%s_tev_sources_sed.fits'%table_prefix)

    source_full = tab['SOURCE_FULL'][0]

    rows_sed_tev = load_source_rows(tab_sed_tev, [tab['SOURCE_FULL'][0]], key='SOURCE_FULL')
    rows_sed_gev = load_source_rows(tab_casc_sed, [tab['name'][0]], key='NAME')

    sed_prim0 = SED.create_from_row(rows_sed_tev)
    sed_prim1 = SED.create_from_row2(rows_sed_gev, tab_ebounds)

    plt.figure()

    
    
    sed_prim0.plot(color='grey',marker='s',label='H.E.S.S. (2005-2006)')                   
    sed_prim1.plot(color='k',marker='o',label='Fermi-LAT')
    #plt.errorbar(sed_prim.ectr/1E6,
    #             sed_prim.ectr*sed_prim.flux,
    #             sed_prim.ectr*sed_prim.flux_err,
    #             marker='o',linestyle='None')


    lcoh_idx = 10

    m = np.isfinite(tab['prim_tev_flux'][0][lcoh_idx,lcoh_idx])




    idxs = np.arange(0,40,5)

    for i, idx in enumerate(idxs):

        efct0 = sed_prim0.ectr**2/sed_prim0.ewidth*1.602177e-06
        efct1 = sed_prim1.ectr**2/sed_prim1.ewidth*1.602177e-06
        
        plt.plot(sed_prim0.ectr/1E6, efct0*tab['prim_tev_flux'][0][lcoh_idx,idx][m],
                 label='B = $10^{%.1f}$ G'%(igmf[lcoh_idx,idx]),
                 color=matplotlib.cm.plasma(float(i)/float(len(idxs))))

        plt.plot(sed_prim1.ectr/1E6, efct1*tab['prim_flux'][0][lcoh_idx,idx],
                 #label='B = %.2f'%(igmf[lcoh_idx,idx]),
                 color=matplotlib.cm.plasma(float(i)/float(len(idxs))))

        plt.plot(sed_prim1.ectr/1E6, efct1*tab['casc_flux'][0][lcoh_idx,idx],
                 #label='B = %.2f'%(igmf[lcoh_idx,idx]),
                 color=matplotlib.cm.plasma(float(i)/float(len(idxs))),linestyle='--')


        
    plt.gca().set_xscale('log')
    plt.gca().set_yscale('log')
    plt.gca().set_xlabel('Energy [TeV]')
    plt.gca().set_ylabel('$\\nu\\mathrm{F}_\\nu$ [erg cm$^{-2}$ s$^{-1}$]')

    plt.legend(frameon=False,ncol=3,prop={'size' : 12})
    plt.gca().set_ylim(1E-8*1E-6,1E-4*1E-6)
    plt.gca().set_xlim(1E-3,30.)

    preliminary()
    
    plt.savefig('%s_igmf_flux.png'%prefix)

    plt.figure()

    for i, idx in enumerate(idxs[3:]):
        plt.plot(sed_prim1.ectr/1E6, tab['casc_r68'][0][lcoh_idx,idx],
                 label='B = $10^{%.1f}$ G'%(igmf[lcoh_idx,idx]),
                 color=matplotlib.cm.plasma(float(i)/float(len(idxs[3:]))))

    plt.plot(sed_prim1.ectr/1E6,np.sqrt((2.85*(sed_prim1.ectr/100.)**-0.8)**2+0.035**2),color='k',label='P8R2_SOURCE_V6::PSF3')

    #print np.sqrt((2.85*(sed_prim1.ectr/100.)**-0.8)**2+0.035**2)
    
    plt.legend(frameon=False,ncol=2,prop={'size' : 12})
    plt.gca().set_xscale('log')
    plt.gca().set_yscale('log')

    plt.gca().set_xlabel('Energy [TeV]')
    plt.gca().set_ylabel('R$_{68}$ [deg]')

    plt.gca().set_ylim(1E-2,10.)
    plt.gca().set_xlim(1E-3,1.)

    preliminary()
    
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
    plt.pcolormesh(10**lcoh,10**igmf,dloglike,cmap='BuGn_r',vmax=0.0,vmin=-64)#,shading='gouraud')
    cb = plt.colorbar()
    cb.set_label('Delta-LogLikelihood')
    plt.contour(10**lcoh,10**igmf,dloglike,
                levels=[-3.0,3.0],colors=['k'],label='test')
    plt.gca().set_xscale('log')
    plt.gca().set_yscale('log')

    plt.plot(0,0,label='95% C.L. Exclusion Region',color='k',linestyle='--')

    plt.gca().set_xlabel('L$_{\mathrm{coh}}$ [Mpc]')
    plt.gca().set_ylabel('B$_{\mathrm{IGMF}}$ [G]')

    plt.gca().legend(frameon=True,loc='lower right',prop={'size' : 12})
    
    preliminary()
    
    plt.savefig('%s_igmf_constraints.png'%prefix)

    #plt.close('all')
