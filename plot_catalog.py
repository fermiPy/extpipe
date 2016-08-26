from astropy.table import Table
import matplotlib.pyplot as plt
import sys
import numpy as np

tab = Table.read(sys.argv[1])


m = (tab['ts'] > 9.0)
m &= (tab['name'] != '3FGL J0534.5+2201')
m &= (tab['name'] != '3FGL J0534.5+2201s')
tab = tab[m]

m0 = (tab['fit_ext_ts'] > 1.0)
m1 = (tab['fit_ext_ts'] < 1.0)

m2 = (tab['fit_halo_ts'] > 16.0)
m3 = (tab['fit_halo_ts'] < 16.0)

halo_width = tab['fit_halo_width']

glat = np.abs(tab['glat'])
halo_ts = tab['fit_halo_ts']
ext = tab['fit_ext_mle']
ext_err = tab['fit_ext_err']
ext_ul95 = tab['fit_ext_ul95']
eflux1000 = tab['eflux1000']
eflux10000 = tab['eflux10000']


plt.figure()
plt.errorbar(eflux1000[m3],halo_ts[m3]**0.5,linestyle='None',marker='o',color='grey',alpha=0.3)
plt.errorbar(eflux1000[m2],halo_ts[m2]**0.5,linestyle='None',marker='o',color='r')

plt.gca().set_xscale('log')
plt.gca().set_xlim(1E-7,1E-3)
plt.gca().set_ylim(0.0,10.0)
plt.gca().set_xlabel('Energy Flux (E > 1 GeV) [MeV cm$^{-2}$ s$^{-1}$]')
plt.gca().set_ylabel('sqrt(TS)')

plt.savefig('halo_eflux1000_sqrt_ts.png')

plt.figure()
plt.errorbar(eflux1000[m3],glat[m3],linestyle='None',marker='o',color='grey',alpha=0.3)
plt.errorbar(eflux1000[m2],glat[m2],linestyle='None',marker='o',color='r')

plt.gca().set_xscale('log')
plt.gca().set_xlim(1E-7,1E-3)
plt.gca().set_ylim(0.0,90.0)
plt.gca().set_xlabel('Energy Flux (E > 1 GeV) [MeV cm$^{-2}$ s$^{-1}$]')
plt.gca().set_ylabel('abs(GLAT)')

plt.savefig('halo_eflux1000_abs_glat.png')

plt.figure()
plt.errorbar(eflux10000[m1],ext_ul95[m1],
             yerr=0.2*ext_ul95[m1],linestyle='None',marker='o',
             color='grey',uplims=True,alpha=0.3)
plt.errorbar(eflux10000[m0],ext[m0],yerr=ext_err[m0],
             linestyle='None',marker='o',
             color='r')

plt.gca().set_xscale('log')
plt.gca().set_yscale('log')
plt.gca().set_xlim(1E-10,1E-3)
plt.gca().set_ylim(0.01,2.0)
plt.gca().set_xlabel('Energy Flux (E > 10 GeV) [MeV cm$^{-2}$ s$^{-1}$]')
plt.gca().set_ylabel('Extension [deg]')

plt.savefig('ext_eflux10000_ext.png')

plt.figure()
plt.errorbar(eflux1000[m1],ext_ul95[m1],
             yerr=0.2*ext_ul95[m1],linestyle='None',marker='o',
             color='grey',uplims=True,alpha=0.3)
plt.errorbar(eflux1000[m0],ext[m0],yerr=ext_err[m0],linestyle='None',marker='o',
             color='r')

plt.gca().set_xscale('log')
plt.gca().set_yscale('log')
plt.gca().set_xlim(1E-7,1E-3)
plt.gca().set_ylim(0.01,2.0)
plt.gca().set_xlabel('Energy Flux (E > 1 GeV) [MeV cm$^{-2}$ s$^{-1}$]')
plt.gca().set_ylabel('Extension [deg]')

plt.savefig('ext_eflux1000_ext.png')

plt.show()
