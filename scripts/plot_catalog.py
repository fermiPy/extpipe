from astropy.table import Table
import matplotlib.pyplot as plt
import sys
import numpy as np
from matplotlib.offsetbox import AnchoredText
from haloanalysis.utils import load_source_rows

kGreen  = (0,1,0) 
kYellow = (1,1,0) 

def cumhist(x,bins):

    xy = np.histogram(x,bins=bins)
    counts = xy[0].astype(float)
    cumsum = np.cumsum(counts)/np.sum(counts)    
    return cumsum, xy[1][1:]
    
#def hist_to_pval(hist,x):
#    np.interp(x,hist
#              np.cumsum(fit_halo_scan_ts_max_xy[1].astype(float))/np.sum(xy[0].astype(float)),xy[1][1:])
    
def preliminary(ax=None,loc=3):

    if ax is None:
        ax = plt.gca()
    
    title = AnchoredText('Preliminary', loc=loc,
                         prop={'size':18,'color': 'red'},
                         pad=0., borderpad=0.75,
                         frameon=False)
    ax.add_artist(title)


tab = Table.read(sys.argv[1])
halo_scan_width = np.array(Table.read(sys.argv[1],hdu='SCAN_PARS')['halo_scan_width'])[0]



m = (tab['ts'] > 9.0)
m &= (tab['name'] != '3FGL J0534.5+2201')
m &= (tab['name'] != '3FGL J0534.5+2201s')
tab = tab[m]

msk_bll = (tab['class'] == 'BLL') | (tab['class'] == 'bll') #| (tab['class'] == 'fsrq') | (tab['class'] == 'FSRQ')
m0 = (tab['fit_ext_ts'] > 16.0)
m1 = (tab['fit_ext_ts'] < 16.0)

m2 = (tab['fit_halo_ts'] > 16.0)
m3 = (tab['fit_halo_ts'] < 16.0)

fit_halo_scan_ts_q84 = np.zeros(len(halo_scan_width))
fit_halo_scan_ts_q97 = np.zeros(len(halo_scan_width))
fit_halo_scan_ts_max_hist = cumhist(np.max(tab[msk_bll]['fit_halo_scan_ts'][:,:,2],axis=1),bins=np.linspace(0,25.,401))



for i in range(len(halo_scan_width)):

    hist = cumhist(tab[msk_bll]['fit_halo_scan_ts'][:,i,2],bins=np.linspace(0,25.,401))
    fit_halo_scan_ts_q84[i] = np.interp(0.84,hist[0],hist[1])
    fit_halo_scan_ts_q97[i] = np.interp(0.975,hist[0],hist[1])

    print i, fit_halo_scan_ts_q84[i], fit_halo_scan_ts_q97[i]
    
#sys.exit(0)
    

halo_width = tab['fit_halo_width']

#glat = np.abs(tab['glat'])
glat = tab['glat']
glon = tab['glon']
halo_ts = tab['fit_halo_ts']
halo_eflux = tab['fit_halo_eflux']
halo_eflux_ul95 = tab['fit_halo_eflux_ul95']
ext = tab['fit_ext_mle']
ext_err = tab['fit_ext_err']
ext_ul95 = tab['fit_ext_ul95']
eflux1000 = tab['eflux1000']
eflux10000 = tab['eflux10000']
flux1000 = tab['flux1000']
flux10000 = tab['flux10000']

glon[glon > 180.] -= 360.


# TEV Sources
names = ['3FGL J0013.9-1853',
         '3FGL J0152.6+0148',
         '3FGL J0232.8+2016',
         '3FGL J0316.6+4119',
         '3FGL J0319.8+1847',
         '3FGL J0349.2-1158',
         '3FGL J0416.8+0104',
         '3FGL J0521.7+2113',
         '3FGL J0648.8+1516',
         '3FGL J0710.3+5908',
         '3FGL J1015.0+4925',
         '3FGL J1103.5-2329',
         '3FGL J1104.4+3812',
         '3FGL J1136.6+7009',
         '3FGL J1314.7-4237',
         '3FGL J1428.5+4240',
         '3FGL J1653.9+3945',
         '3FGL J2000.0+6509',
         '3FGL J2158.8-3013',
         '3FGL J2250.1+3825',
         '3FGL J2347.0+5142',
         '3FGL J2359.3-3038']

rows = load_source_rows(tab, names, key='name')



plt.figure(figsize=(12,8))
plt.subplot(111, projection="aitoff")
plt.plot(-np.radians(glon),np.radians(glat),linestyle='None',
             marker='o',color='grey',alpha=0.3,label='__nolabel__')
plt.plot(-np.radians(rows['glon']),np.radians(rows['glat']),linestyle='None',
             marker='o',color='r',label='TeV Sources')

plt.gca().legend(frameon=False,prop={'size' : 12})
#plt.ylim(reversed(plt.ylim()))
plt.gca().grid(True)

lon_ticks = [150,120,90,60,30,0,-30,-60,-90,-120,-150]
plt.gca().set_xticklabels(['%i$^\circ$'%t for t in lon_ticks])
plt.tight_layout()
plt.savefig('catalog_tev_sources_glon_glat.png')

plt.figure()

plt.plot(rows['eflux10000'],rows['fit_halo_ts']**0.5,marker='o',linestyle='None')

plt.gca().set_xscale('log')
plt.gca().set_xlim(1E-7,1E-3)

plt.figure()

rows.sort('fit_halo_ts')
rows.reverse()

for row in rows[rows['fit_halo_ts'] < 3]:
    plt.plot(halo_scan_width,row['fit_halo_scan_ts'][:,2]**0.5,marker='None',linestyle='-',color='grey')

for row in rows[rows['fit_halo_ts'] > 3]:

    if row['assoc'] == 'TXS 0518+211':
        assoc = 'VER J0521+211'
    else:
        assoc = row['assoc']
    
    plt.plot(halo_scan_width,row['fit_halo_scan_ts'][:,2]**0.5,marker='None',linestyle='-',
             label=assoc,linewidth=1.5)


    max_ts = np.max(row['fit_halo_scan_ts'][:,2])
    pval_local =  (1.0-np.interp(max_ts, fit_halo_scan_ts_max_hist[1], fit_halo_scan_ts_max_hist[0])**len(rows))
    pval_global =  (1.0-np.interp(max_ts, fit_halo_scan_ts_max_hist[1], fit_halo_scan_ts_max_hist[0]))
    print row['assoc'], max_ts, pval_local, pval_global
    
plt.gca().fill_between(halo_scan_width,fit_halo_scan_ts_q84**0.5,
                       color=kGreen,alpha=1.0,zorder=0.9,
                       label='84% Containment')
plt.gca().fill_between(halo_scan_width,fit_halo_scan_ts_q97**0.5,
                       color=kYellow,alpha=1.0,zorder=0.8,
                       label='97.5% Containment')
    
#plt.plot(halo_scan_width,np.sum(rows['fit_halo_scan_ts'],axis=0)[:,2],marker='None',linestyle='-',
#         label='composite',color='k')
    
plt.gca().legend(frameon=False,prop={'size' : 12},ncol=2)

plt.gca().set_xscale('log')
plt.gca().set_xlim(halo_scan_width[0],halo_scan_width[-1])
plt.gca().set_ylim(0,5)

plt.gca().set_xlabel('Halo Width [deg]')
plt.gca().set_ylabel('$(\mathrm{TS}_{\mathrm{halo}})^{1/2}$')
preliminary(loc=2)
plt.savefig('catalog_tev_sources_halo_ts.png')




#-----------------------------------------------
# Halo GLON vs. GLAT


plt.figure(figsize=(12,8))
plt.subplot(111, projection="aitoff")
plt.plot(-np.radians(glon[m3]),np.radians(glat[m3]),linestyle='None',
             marker='o',color='grey',alpha=0.3,label='__nolabel__')
plt.plot(-np.radians(glon[m2]),np.radians(glat[m2]),linestyle='None',
             marker='o',color='r',label='TS$_{\mathrm{halo}}$ > 16')
plt.plot(-np.radians(glon[m0]),np.radians(glat[m0]),linestyle='None',
          marker='o',color='None',markeredgecolor='b', markeredgewidth=1.25,
          label='TS$_{\mathrm{ext}}$ > 16')


plt.gca().legend(frameon=False,prop={'size' : 12})
plt.gca().grid(True)

lon_ticks = [150,120,90,60,30,0,-30,-60,-90,-120,-150]
plt.gca().set_xticklabels(['%i$^\circ$'%t for t in lon_ticks])
preliminary()
plt.tight_layout()
plt.savefig('catalog_halo_glon_glat.png')

sys.exit(0)


#-----------------------------------------------
# Halo EFlux vs. TS

plt.figure()
plt.errorbar(halo_eflux[m3],halo_ts[m3]**0.5,linestyle='None',marker='o',color='grey',alpha=0.3,label='TS$_{\mathrm{halo}}$ < 9')
plt.errorbar(halo_eflux[m2],halo_ts[m2]**0.5,linestyle='None',marker='o',color='r',label='TS$_{\mathrm{halo}}$ > 9')

plt.gca().legend(frameon=False)
plt.gca().set_xscale('log')
plt.gca().set_xlim(1E-7,1E-3)
plt.gca().set_ylim(0.0,10.0)
plt.gca().set_xlabel('Energy Flux (E > 1 GeV) [MeV cm$^{-2}$ s$^{-1}$]')
plt.gca().set_ylabel('sqrt(TS)')

preliminary()

plt.savefig('catalog_halo_eflux1000_sqrt_ts.png')

#-----------------------------------------------
# Halo EFlux vs. GLAT

plt.figure()
plt.errorbar(glat[m3],halo_eflux_ul95[m3],yerr=0.2*halo_eflux_ul95[m3],linestyle='None',
             marker='o',color='grey',alpha=0.3,label='TS$_{\mathrm{halo}}$ < 16',
             uplims=True)
plt.errorbar(glat[m2],halo_eflux[m2],linestyle='None',
             marker='o',color='r',label='TS$_{\mathrm{halo}}$ > 16')

plt.gca().set_yscale('log')
plt.gca().set_ylim(3E-8,1E-4)
plt.gca().set_xlim(-90.0,90.0)
plt.gca().set_ylabel('Halo Energy Flux (1 - 316 GeV) [MeV cm$^{-2}$ s$^{-1}$]')
plt.gca().set_xlabel('GLAT')
plt.gca().legend(frameon=False,prop={'size' : 12})

preliminary()

plt.savefig('catalog_halo_eflux1000_abs_glat.png')

#-----------------------------------------------
# Extension vs. EFlux > 10 GeV

plt.figure()
plt.errorbar(eflux10000[m1],ext_ul95[m1],
             yerr=0.2*ext_ul95[m1],linestyle='None',marker='o',
             color='grey',uplims=True,alpha=0.3)
plt.errorbar(eflux10000[m0],ext[m0],yerr=ext_err[m0],
             linestyle='None',marker='o',
             color='r')
plt.gca().axhline(0.0276,color=plt.cm.Blues(0.4),linestyle='--',label='Syst. Uncertainty ($\Gamma = 2.4$)')
plt.gca().axhline(0.0195,color=plt.cm.Blues(1.0),linestyle='--',label='Syst. Uncertainty ($\Gamma = 1.6$)')

plt.gca().set_xscale('log')
plt.gca().set_yscale('log')
plt.gca().set_xlim(1E-10,1E-3)
plt.gca().set_ylim(0.01,2.0)
plt.gca().set_xlabel('Energy Flux (E > 10 GeV) [MeV cm$^{-2}$ s$^{-1}$]')
plt.gca().set_ylabel('Extension [deg]')
plt.gca().legend(frameon=False,prop={'size' : 12})

preliminary()

plt.savefig('catalog_ext_eflux10000_ext.png')



#-----------------------------------------------
# Extension vs. EFlux > 1 GeV

plt.figure()
plt.errorbar(eflux1000[m1],ext_ul95[m1],
             yerr=0.2*ext_ul95[m1],linestyle='None',marker='o',
             color='grey',uplims=True,alpha=0.3,label='TS$_{\mathrm{ext}}$ < 16')
plt.errorbar(eflux1000[m0],ext[m0],yerr=ext_err[m0],linestyle='None',marker='o',
             color='r',label='TS$_{\mathrm{ext}}$ > 16')
plt.gca().axhline(0.0276,color=plt.cm.Blues(0.4),linestyle='--',label='Syst. Uncertainty ($\Gamma = 2.4$)')
plt.gca().axhline(0.0247,color=plt.cm.Blues(0.7),linestyle='--',label='Syst. Uncertainty ($\Gamma = 2.0$)')
plt.gca().axhline(0.0195,color=plt.cm.Blues(1.0),linestyle='--',label='Syst. Uncertainty ($\Gamma = 1.6$)')

plt.gca().set_xscale('log')
plt.gca().set_yscale('log')
plt.gca().set_xlim(1E-7,1E-3)
plt.gca().set_ylim(0.01,2.0)
plt.gca().set_xlabel('Energy Flux (E > 1 GeV) [MeV cm$^{-2}$ s$^{-1}$]')
plt.gca().set_ylabel('Extension [deg]')
plt.gca().legend(frameon=False,prop={'size' : 12})

preliminary()

plt.savefig('catalog_ext_eflux1000_ext.png')

#-----------------------------------------------
# Extension vs. Flux > 10 GeV

plt.figure()
plt.errorbar(flux10000[m1],ext_ul95[m1],
             yerr=0.2*ext_ul95[m1],linestyle='None',marker='o',
             color='grey',uplims=True,alpha=0.3,label='TS$_{\mathrm{ext}}$ < 16')
plt.errorbar(flux10000[m0],ext[m0],yerr=ext_err[m0],
             linestyle='None',marker='o',
             color='r',label='TS$_{\mathrm{ext}}$ > 16')
plt.gca().axhline(0.0276,color=plt.cm.Blues(0.4),linestyle='--',label='Syst. Uncertainty ($\Gamma = 2.4$)')
plt.gca().axhline(0.0247,color=plt.cm.Blues(0.7),linestyle='--',label='Syst. Uncertainty ($\Gamma = 2.0$)')
plt.gca().axhline(0.0195,color=plt.cm.Blues(1.0),linestyle='--',label='Syst. Uncertainty ($\Gamma = 1.6$)')

plt.gca().set_xscale('log')
plt.gca().set_yscale('log')
plt.gca().set_xlim(1E-13,1E-7)
plt.gca().set_ylim(0.01,2.0)
plt.gca().set_xlabel('Flux (E > 10 GeV) [cm$^{-2}$ s$^{-1}$]')
plt.gca().set_ylabel('Extension [deg]')
plt.gca().legend(frameon=False,prop={'size' : 12})

preliminary()

plt.savefig('catalog_ext_flux10000_ext.png')

#-----------------------------------------------
# Extension vs. Flux > 10 GeV

plt.figure()
plt.errorbar(flux1000[m1],ext_ul95[m1],
             yerr=0.2*ext_ul95[m1],linestyle='None',marker='o',
             color='grey',uplims=True,alpha=0.3,label='TS$_{\mathrm{ext}}$ < 16')
plt.errorbar(flux1000[m0],ext[m0],yerr=ext_err[m0],
             linestyle='None',marker='o',
             color='r',label='TS$_{\mathrm{ext}}$ > 16')
plt.gca().axhline(0.0276,color=plt.cm.Blues(0.4),linestyle='--',label='Syst. Uncertainty ($\Gamma = 2.4$)')
plt.gca().axhline(0.0247,color=plt.cm.Blues(0.7),linestyle='--',label='Syst. Uncertainty ($\Gamma = 2.0$)')
plt.gca().axhline(0.0195,color=plt.cm.Blues(1.0),linestyle='--',label='Syst. Uncertainty ($\Gamma = 1.6$)')

plt.gca().set_xscale('log')
plt.gca().set_yscale('log')
plt.gca().set_xlim(1E-13,1E-7)
plt.gca().set_ylim(0.01,2.0)
plt.gca().set_xlabel('Flux (E > 10 GeV) [cm$^{-2}$ s$^{-1}$]')
plt.gca().set_ylabel('Extension [deg]')
plt.gca().legend(frameon=False,prop={'size' : 12})

preliminary()

plt.savefig('catalog_ext_flux1000_ext.png')

plt.figure()
plt.errorbar(flux1000[m1],ext_ul95[m1],
             yerr=0.2*ext_ul95[m1],linestyle='None',marker='o',
             color='grey',uplims=True,alpha=0.3,label='TS$_{\mathrm{ext}}$ < 16')
plt.errorbar(flux1000[m0],ext[m0],yerr=ext_err[m0],linestyle='None',marker='o',
             color='r',label='TS$_{\mathrm{ext}}$ > 16')

plt.gca().set_xscale('log')
plt.gca().set_yscale('log')
plt.gca().set_xlim(1E-12,1E-6)
plt.gca().set_ylim(0.01,2.0)
plt.gca().set_xlabel('Flux (E > 1 GeV) [cm$^{-2}$ s$^{-1}$]')
plt.gca().set_ylabel('Extension [deg]')
plt.gca().legend(frameon=False)

plt.savefig('catalog_ext_flux1000_ext.png')

plt.show()
