import sys
import numpy as np
from astropy.table import Table, Column
from gammatools.core.histogram import *
from scipy.stats import chi2



def make_dist_hist(tabs,axis,masks,cols,labels,xlabel='',ylabel='',log=False):
        
    plt.figure()
    for t,m,l,c in zip(tabs,masks,labels,cols):
        
        h = Histogram(axis)
        h.fill(t[m][c])
        h = h.normalize()
        h.plot(hist_style='stepfilled',alpha=0.3,linewidth=1,label=l)

        print h.quantile(0.5)

    plt.gca().legend(frameon=False)
    plt.gca().set_xlabel(xlabel)
    plt.gca().set_ylabel(ylabel)

    if log:
        plt.gca().set_xscale('log')
    

def make_ts_hist(tab,masks,cols,labels,cumulative=True):

    if not isinstance(cols,list):
        cols = [cols]*len(masks)
    
    plt.figure()
    for m,l,c in zip(masks,labels,cols):
        h = Histogram(Axis.create(0.0,20,100))
        h.fill(tab[m][c])
        h = h.normalize()
        if cumulative:
            h = h.cumulative(lhs=False)
        h.plot(hist_style='step',alpha=0.3,linewidth=2,label=l)
        
    dof  = 2
    label = r"$\chi^2_{1} / 2$"
    kwargs = dict( label = label, lw=1.5, c='k',dashes=(5,2))

    
    if cumulative:
        plt.gca().plot(h.axis(0).center,
                       0.5*(1-chi2.cdf(h.axis(0).edges[:-1],1)),**kwargs)
    else:
        plt.gca().plot(h.axis(0).center,
                       h.axis(0).width*chi2.pdf(h.axis(0).center,1),**kwargs)
        
    label = r"$\chi^2_{2} / 2$"
    kwargs = dict( label = label, lw=1.5, c='r',dashes=(5,2))
#    plt.gca().plot(h.axis(0).center,
#                   0.5*(1-chi2.cdf(h.axis(0).center,2)),**kwargs)

#    if cumulative:
#        plt.gca().plot(h.axis(0).center,
#                       0.5*(1-chi2.cdf(h.axis(0).edges[:-1],2)),**kwargs)
#    else:
#        plt.gca().plot(h.axis(0).center,
#                       h.axis(0).width*chi2.pdf(h.axis(0).center,2),**kwargs)
    
    plt.gca().set_yscale('log')
    plt.gca().set_ylim(1E-4,1)
    plt.gca().legend(frameon=False)
    plt.gca().set_xlabel('TS$_\mathrm{ext}$')
    plt.gca().set_ylabel('Cumulative Fraction')


tab0 = Table.read(sys.argv[1])
tab1 = Table.read(sys.argv[2])

m0 = (tab0['npred'] > 3.) & (tab0['ts'] > 9.) & (tab0['ts'] < 10000.) & (np.abs(tab0['glat']) > 5.) 
m1 = (tab1['npred'] > 3.) & (tab1['ts'] > 9.) & (tab1['ts'] < 10000.) & (np.abs(tab1['glat']) > 5.) 
#m0 &= ((tab0['class'] == 'bll') | (tab0['class'] == 'fsrq'))
#m0 &= (tab0['dfde1000_index'] < 1.8)

m0_hlat = (np.abs(tab0['glat']) > 10.)
m0_bll = ((tab0['class'] == 'bll') | (tab0['class'] == 'BLL'))
m0_fsrq = ((tab0['class'] == 'fsrq') | (tab0['class'] == 'FSRQ'))
m0_psr = ((tab0['class'] == 'psr') | (tab0['class'] == 'PSR'))
m0_bcu = ((tab0['class'] == 'bcu') | (tab0['class'] == 'BCU'))
m0_unkn = ((tab0['class'] == 'unkn'))

m1_hlat = (np.abs(tab1['glat']) > 10.)
m1_bll = ((tab1['class'] == 'bll') | (tab1['class'] == 'BLL'))
m1_fsrq = ((tab1['class'] == 'fsrq') | (tab1['class'] == 'FSRQ'))
m1_psr = ((tab1['class'] == 'psr') | (tab1['class'] == 'PSR'))
m1_bcu = ((tab1['class'] == 'bcu') | (tab1['class'] == 'BCU'))
m1_unkn = ((tab1['class'] == 'unkn'))
        

m0_nsrc2 = (tab0['fit_nsrc'] == 2)
m0_nsrc3 = (tab0['fit_nsrc'] == 3)
m0_nsrc4 = (tab0['fit_nsrc'] == 4)

plt.figure()
plt.scatter(2*tab0['fit_dlike_ext'][m0_nsrc2],tab0['fitn_ext_ts'][m0_nsrc2][:,1],color='b',edgecolor='k',label='2 PS')
plt.scatter(2*tab0['fit_dlike_ext'][m0_nsrc3],tab0['fitn_ext_ts'][m0_nsrc3][:,2],color='r',edgecolor='k',label='3 PS')
plt.scatter(2*tab0['fit_dlike_ext'][m0_nsrc4],tab0['fitn_ext_ts'][m0_nsrc4][:,3],color='g',edgecolor='k',label='4 PS')
plt.axvline(0.0,color='k',linestyle='--',zorder=-10)
plt.axvline(-6.0,color='r',linestyle='--',zorder=-10)
plt.axhline(0.0,color='k',zorder=-10)
plt.gca().set_xlim(-30,30)
plt.gca().set_ylim(-5,50)
plt.gca().set_xlabel('2$\\times$ Delta LogLikelihood (N PS vs. N-1 PS + Extension)')
plt.gca().set_ylabel('TS$_{\mathrm{ext}}$')

plt.gca().legend(frameon=False)

plt.savefig('fit_dlike_ext.png')




plt.figure()
plt.scatter(2*tab0['fit_dlike_halo'][m0_nsrc2],tab0['fitn_halo_ts'][m0_nsrc2][:,1],color='b',edgecolor='k',label='2 PS')
plt.scatter(2*tab0['fit_dlike_halo'][m0_nsrc3],tab0['fitn_halo_ts'][m0_nsrc3][:,2],color='r',edgecolor='k',label='3 PS')
plt.scatter(2*tab0['fit_dlike_halo'][m0_nsrc4],tab0['fitn_halo_ts'][m0_nsrc4][:,3],color='g',edgecolor='k',label='4 PS')
plt.axvline(0.0,color='k',linestyle='--',zorder=-10)
plt.axvline(-4.0,color='r',linestyle='--',zorder=-10)
plt.axhline(0.0,color='k',zorder=-10)
plt.gca().set_xlim(-30,30)
plt.gca().set_ylim(-5,50)
plt.gca().set_xlabel('2$\\times$ Delta LogLikelihood (N PS vs. N-1 PS + Halo)')
plt.gca().set_ylabel('TS$_{\mathrm{halo}}$')

plt.gca().legend(frameon=False)

plt.savefig('fit_dlike_halo.png')

plt.show()

#h0.plot(hist_style='stepfilled',alpha=0.3,label='nominal position',color='k')
#h1.plot(hist_style='stepfilled',alpha=0.3,label='with localization')

#h1.plot(hist_style='step',alpha=0.3,label='relocalized',linewidth=2)
#h2.plot(hist_style='step',alpha=0.3,label='multiple',linewidth=2)


eflux_kwargs = dict(xlabel='Halo Flux 95% CL Upper Limit [cm$^{-2}$ s$^{-1}$]',log=True)
ext_kwargs = dict(xlabel='Extension 95% CL Upper Limit [deg]',log=True)


if 0:
    make_dist_hist([tab0,tab1],Axis(10**np.linspace(-11,-8,100)),
                   [m0&m1&m0_hlat,m0&m1&m1_hlat],
                   ['fit_halo_2.00_1.000_eflux_ul95','fit_halo_2.00_1.000_eflux_ul95'],['E > 1 GeV','E > 10 GeV'],
                   **eflux_kwargs)

    plt.savefig('halo_flux_ul_1.000.png')

    make_dist_hist([tab0,tab1],Axis(10**np.linspace(-11,-8,100)),
                   [m0&m1&m0_hlat,m0&m1&m1_hlat],
                   ['fit_halo_2.00_0.316_eflux_ul95','fit_halo_2.00_0.316_eflux_ul95'],['E > 1 GeV','E > 10 GeV'],
                   **eflux_kwargs)

    plt.savefig('halo_flux_ul_0.316.png')

    make_dist_hist([tab0,tab1],Axis(10**np.linspace(-11,-8,100)),
                   [m0&m1&m0_hlat,m0&m1&m1_hlat],
                   ['fit_halo_2.00_0.100_eflux_ul95','fit_halo_2.00_0.100_eflux_ul95'],['E > 1 GeV','E > 10 GeV'],
                   **eflux_kwargs)

    plt.savefig('halo_flux_ul_0.100.png')


    make_dist_hist([tab0,tab1],Axis(10**np.linspace(-2,0,100)),
                   [m0&m1&m0_hlat,m0&m1&m1_hlat],
                   ['ext_ul95','ext_ul95'],['E > 1 GeV','E > 10 GeV'],
                   **ext_kwargs)

    plt.savefig('ext_ul_0.100.png')


    
    
# Halo Size
if 0:
    make_ts_hist(tab0,
                 [m0&m0_hlat,m0&m0_hlat,m0&m0_hlat],
                 ['fit_halo_2.00_0.100_ts','fit_halo_2.00_0.316_ts','fit_halo_2.00_1.000_ts'],['0.1 deg','0.316 deg','1.0 deg'])

    plt.savefig('halo_ts_size.png')
    
    make_ts_hist(tab0,
                 [m0&m0_hlat,m0&m0_hlat,m0&m0_hlat],
                 ['fit_halo_1.50_0.316_ts','fit_halo_2.00_0.316_ts','fit_halo_2.50_0.316_ts'],['1.5','2.0','2.5','3.0'])

    plt.savefig('halo_ts_index.png')

    
# BLL vs. FSRQ
if 0:
    make_ts_hist(tab0,
                 [m0&m0_fsrq&m0_hlat,m0&m0_bll&m0_hlat,m0&m0_psr&m0_hlat,m0&m0_bcu&m0_hlat,m0&m0_unkn&m0_hlat],
                 ['fit_ext_ts','fit_ext_ts','fit_ext_ts','fit_ext_ts','fit_ext_ts'],['FSRQ','BLLac','PSR','BCU','Unknown'])


    plt.savefig('fit_ext_ts_class_1gev.png')
    
    make_ts_hist(tab1,
                 [m1&m1_fsrq&m1_hlat,m1&m1_bll&m1_hlat,m1&m1_psr&m1_hlat,m1&m1_bcu&m1_hlat,m1&m1_unkn&m1_hlat],
                 ['fit_ext_ts','fit_ext_ts','fit_ext_ts','fit_ext_ts','fit_ext_ts'],['FSRQ','BLLac','PSR','BCU','Unknown'])

    plt.savefig('fit_ext_ts_class_10gev.png')

    make_ts_hist(tab0,
                 [m0&(tab0['ts'] < 25.),m0&(tab0['ts'] < 100.)&(tab0['ts'] > 25.),
                  m0&(tab0['ts'] < 300.)&(tab0['ts'] > 100.),m0&(tab0['ts'] < 1000.)&(tab0['ts'] > 300.),m0&(tab0['ts'] > 1000.)],
                 'fit_ext_ts',['9 < ts < 25','25 < ts < 100','100 < ts < 300','300 < ts < 1000','ts > 1000'])

    plt.savefig('fit_ext_ts_source_ts_1gev.png')

# Localize vs. Non Localized
if 0:
    make_ts_hist(tab0,
                 [m0&(tab0['ts'] > 25.)&m0_hlat,m0&(tab0['ts'] > 25.)&m0_hlat,m0&(tab0['ts'] > 25.)&m0_hlat],
                 ['ext0_ts','ext1_ts','fit_ext_ts'],['Nominal','Localized','Multi-Source'])

    plt.savefig('fit_ext_ts_local_1gev.png')
    
    # Localize vs. Non Localized
    make_ts_hist(tab1,
                 [m1&(tab1['ts'] > 25.)&m1_hlat,m1&(tab1['ts'] > 25.)&m1_hlat,m1&(tab1['ts'] > 25.)&m1_hlat],
                 ['ext0_ts','ext1_ts','fit_ext_ts'],['Nominal','Localized','Multi-Source'])

    plt.savefig('fit_ext_ts_local_10gev.png')


if 0:
    make_ts_hist(tab0,
                 [m0&(tab0['glat'] < 30.)&(tab0['ts'] > 25.),m0&(tab0['glat'] > 30.0)&(tab0['ts'] > 25.)],
                 'fit_ext_ts',['|b| < 30.','|b| > 30.'])

    plt.savefig('fit_ext_ts_latitude_1gev.png')
    
#make_ts_hist(tab0,
#            [m0&(tab0['dfde1000_index'] < 2.0),m0&(tab0['dfde1000_index'] > 2.0)],
#            ['fit_ext_ts','fit_ext_ts'],['< 2.0','> 2.0'])

#make_ts_hist(tab0,
#            [m0&(tab0['ts'] > 9.),m0&(tab0['ts'] > 25.),m0&(tab0['ts'] > 100.),m0&(tab0['ts'] > 300.),m0&(tab0['ts'] > 1000.)],
#            'fit_ext_ts',['ts > 9','ts > 25','ts > 100','ts > 300','ts > 1000'])




#make_ts_hist(tab0,
#            [m0&(tab0['ts'] > 9.),m0&(tab0['ts'] > 25.),m0&(tab0['ts'] > 100.),m0&(tab0['ts'] > 300.),m0&(tab0['ts'] > 1000.)],
#            'ext_ts',['ts > 9','ts > 25','ts > 100','ts > 300','ts > 1000'],cumulative=False)

#make_ts_hist(tab0,
#            [m0&(tab0['ts'] > 25.),m0&(tab0['ts'] > 100.),m0&(tab0['ts'] > 300.),m0&(tab0['ts'] > 1000.)],
#            'fit2_halo_2.00_1.000_ts',['ts > 25','ts > 100','ts > 300','ts > 1000'])


#make_ts_hist(tab0,
#            [m0&(tab0['ts'] > 25.),m0&(tab0['ts'] > 100.),m0&(tab0['ts'] > 300.),m0&(tab0['ts'] > 1000.)],
#            'fit2_halo_2.00_0.316_ts',['ts > 25','ts > 100','ts > 300','ts > 1000'])

#make_ts_hist(tab0,
#            [m0&(tab0['ts'] > 25.),m0&(tab0['ts'] > 100.),m0&(tab0['ts'] > 300.),m0&(tab0['ts'] > 1000.)],
#            'fit2_halo_2.00_0.100_ts',['ts > 25','ts > 100','ts > 300','ts > 1000'])


sys.exit(0)

h2 = Histogram(Axis.create(0,0.5,100))
h3 = Histogram(Axis.create(0,0.5,100))

plt.figure()
h2.fill(tab0['ext0_ul95'])
h3.fill(tab0['ext1_ul95'])


h2.plot(hist_style='stepfilled',alpha=0.3,label='nominal position',color='k')
h3.plot(hist_style='stepfilled',alpha=0.3,label='relocalized')
plt.gca().legend(frameon=False)
plt.gca().set_xlabel('Extension UL [deg]')

plt.savefig('extul.png')

h2 = Histogram(Axis.create(0,0.5,100))
h3 = Histogram(Axis.create(0,0.5,100))

plt.figure()
h2.fill(tab0['ext1_ul95'])
h3.fill(tab1['ext1_ul95'])

h3 *= h2.sum()[0]/h3.sum()[0]


h2.plot(hist_style='stepfilled',alpha=0.3,label='combined',color='k')
h3.plot(hist_style='stepfilled',alpha=0.3,label='joint')
plt.gca().legend(frameon=False)

plt.gca().set_xlabel('Extension UL [deg]')
