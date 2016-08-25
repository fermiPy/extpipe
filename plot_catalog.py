from astropy.table import Table
import matplotlib.pyplot as plt
import sys

tab = Table.read(sys.argv[1])


m = (tab['ts'] > 25.0)
m &= (tab['name'] != '3FGL J0534.5+2201')
m &= (tab['name'] != '3FGL J0534.5+2201s')
#m = tab['fit_ext_ts'] > 4.0
tab = tab[m]

m0 = (tab['fit_ext_ts'] > 9.0)
m1 = (tab['fit_ext_ts'] < 9.0)

ext = tab['fit_ext_mle']
ext_err = tab['fit_ext_err']
ext_ul95 = tab['fit_ext_ul95']
eflux1000 = tab['eflux1000']
eflux10000 = tab['eflux10000']


plt.figure()
plt.errorbar(eflux10000[m1],ext_ul95[m1],
             yerr=0.2*ext_ul95[m1],linestyle='None',marker='o',
             color='grey',uplims=True,alpha=0.3)
plt.errorbar(eflux10000[m0],ext[m0],yerr=ext_err[m0],
             linestyle='None',marker='o',
             color='r')

plt.gca().set_xscale('log')
plt.gca().set_yscale('log')
plt.gca().set_xlim(1E-9,1E-3)
plt.gca().set_ylim(0.01,2.0)

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

plt.show()
