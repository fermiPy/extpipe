import sys
import os
import argparse
import math

import numpy as np
import matplotlib.pyplot as plt
from collections import OrderedDict

from astropy.table import Table, Column
from astropy.io import fits
from astropy import units as u
from astropy.coordinates import SkyCoord

from scipy.special import gamma

parser = argparse.ArgumentParser()

parser.add_argument('--output', default = 'test')
parser.add_argument('--input', default = 'test.fits',
                    help='input fits file from aggregate output')
parser.add_argument('--show_plots', dest='show_plots', action='store_true',
                    help='include this argument to display plots')
parser.add_argument('--save_plots', dest='save_plots', action='store_true',
                    help='include this argument to save plots')
parser.set_defaults(show_plots=False)
parser.set_defaults(save_plots=False)

args = parser.parse_args()

#print name of input file
print "Analyzing file:", args.input

hdulist=fits.open(args.input)
data=hdulist[1].data

#list of object classes, full list here: 
#http://heasarc.gsfc.nasa.gov/W3Browse/fermi/fermilpsc.html
objcls=['agn','bcu','BCU','bin','bll','BLL','css','fsrq','FSRQ','gal',
        'glc','hmb','nlsy1','NLSY1','nov','psr','PSR','pwn','rdg','RDG',
        'sbg','sey','sfr','snr','spp','ssrq','unkn']

def column(matrix, i):
    return[row[i] for row in matrix]

def save(path, ext='png', close=True, verbose=True):
    """Save a figure from pyplot.
    Parameters
    ----------
    path : string
        The path (and filename, without the extension) to save the
        figure to.
    ext : string (default='png')
        The file extension. This must be supported by the active
        matplotlib backend (see matplotlib.backends module).  Most
        backends support 'png', 'pdf', 'ps', 'eps', and 'svg'.
    close : boolean (default=True)
        Whether to close the figure after saving.  If you want to save
        the figure multiple times (e.g., to multiple formats), you
        should NOT close it in between saves or you will have to
        re-plot it.
    verbose : boolean (default=True)
        Whether to print information about when and where the image
        has been saved.
    """
        
    # Extract the directory and filename from the given path
    directory = os.path.split(path)[0]
    filename = "%s.%s" % (os.path.split(path)[1], ext)
    if directory == '':
        directory = '.'

    # If the directory does not exist, create it
    if not os.path.exists(directory):
        os.makedirs(directory)

    # The final path to save to
    savepath = os.path.join(directory, filename)

    if verbose:
        print("Saving figure to '%s'..." % savepath),

    # Actually save the figure
    plt.savefig(savepath)
    
    # Close it
    if close:
        plt.close()

    if verbose:
        print("Done")

def chi2func(x, df, norm):
    y=np.where(x>0.,norm*(((x)**(df/2-1) * np.exp(-x/2))) / ((2.0**(df/2.0))*gamma(df/2)),0.)
    return y

plt.figure(figsize=(10,20),facecolor='w',edgecolor='w')
plt.figure(1)

#plot ra, dec of objects
plt.subplot(321,projection="aitoff")
plt.grid(True)
c = SkyCoord(column(data,'ra')*u.degree, column(data,'dec')*u.degree, frame='icrs')
plt.plot(c.ra.wrap_at(180*u.deg).radian,c.dec.wrap_at(180*u.deg).radian,'ro')
plt.ylabel('dec (deg)')
#plt.xlabel('ra (deg)')

#plot glat, glon of objects
plt.subplot(322,projection="aitoff")
plt.grid(True)
gc = SkyCoord(column(data,'glon')*u.degree, column(data,'glat')*u.degree, frame='galactic')
plt.plot(gc.l.wrap_at(180*u.deg).radian,gc.b.wrap_at(180*u.deg).radian,'go')
plt.ylabel('glat (deg)')
#plt.xlabel('glon (deg)')

#plot number of objects in each class
objnum=[]
#print column(data,'class')
for i in range(0,len(objcls)):
    holder=[s for s in column(data,'class') if objcls[i] in s]
    print objcls[i], len(holder)
    objnum.append(len(holder))

class_info = OrderedDict(zip(objcls, objnum))
plt.subplot(312)
plt.bar(range(len(class_info)), class_info.values(),align='center',log=True)
plt.xticks(range(len(class_info)),class_info.keys())
plt.xlabel('class')
locs, labels = plt.xticks()
plt.setp(labels, rotation=90)
plt.xlim([-1.,len(class_info)])
plt.ylim([0.1,2000.])

#plot eflux at 1 GeV vs. dfde index at 1 GeV
plt.subplot(313)
plt.plot(column(data,'eflux1000'),column(data,'dfde1000_index'),'bo')
plt.ylabel('dfde 1GeV index')
plt.xlabel('eflux 1GeV')
plt.xscale('log')

if args.save_plots:
    save('figures/%s' % args.output,'png',False,True)

#ts distribution vs. chi2(1) distribution
ext_list=['ext0_ts','ext1_ts','ext2_ts']
num_ext=len(ext_list)
width=num_ext*5
plt.figure(figsize=(width,5),facecolor='w',edgecolor='w')
plt.figure(2)

for i in range(0,num_ext):
    ax=plt.subplot(1,num_ext,i+1)
    ext_data=column(data,'%s' % ext_list[i])
    clean_data=[x for x in ext_data if not math.isnan(x)]#remove nan from data
    #print clean_data
    n, bins, patches = plt.hist(clean_data, 
                                int(math.ceil( max(column(data,'%s' % ext_list[i])) )),
                                normed=1, facecolor='green')
    bincenters = 0.5*(bins[1:]+bins[:-1])
    chi2_1vals=chi2func(bincenters, df=1.0, norm=1.0)
    chi2_2vals=chi2func(bincenters, df=2.0, norm=1.0)
    plt.plot(bincenters,chi2_1vals,'r--',linewidth=2.0,label='$\chi^2_1$/2')
    plt.plot(bincenters,chi2_2vals,'b--',linewidth=2.0,label='$\chi^2_2$/2')
    legend = ax.legend(loc='upper right',frameon=False)
    plt.ylabel('PDF')
    plt.xlabel('TS$_{%s}$' % ext_list[i][0:4])
    plt.yscale('log')
    plt.ylim([0.00001,1.])

if args.save_plots:
    save('figures/%s_ts' % args.output,'png',False,True)

#loop over and make figures for each specific class of objects
fignum=3
#for i in range(0):
for i in range(0,len(objcls)):
    class_data=[s for s in data if objcls[i] in s]
    if len(class_data)==0:
        continue

    plt.figure(figsize=(width,10),facecolor='w',edgecolor='w')
    plt.figure(fignum)
    #plot glat, glon of objects
    plt.subplot(2,num_ext,1,projection="aitoff")
    plt.grid(True)
    gc = SkyCoord(column(class_data,'glon')*u.degree, column(class_data,'glat')*u.degree, frame='galactic')
    plt.plot(gc.l.wrap_at(180*u.deg).radian,gc.b.wrap_at(180*u.deg).radian,'go')
    plt.ylabel('glat (deg)')
    #plt.xlabel('glon (deg)')
    
    #plot eflux vs dfde index
    plt.subplot(2,num_ext,2)
    plt.plot(column(class_data,'eflux1000'),column(class_data,'dfde1000_index'),'bo')
    plt.ylabel('dfde 1GeV index')
    plt.xlabel('eflux 1GeV')
    plt.xscale('log')
    plt.title('%s' % objcls[i])

    plt.subplot(2,num_ext,3)
    plt.plot(column(class_data,'ext0_ts'),column(class_data,'ext1_ts'),'ro')
    plt.ylabel('ext1 ts')
    plt.xlabel('ext0 ts')
    plt.title('%s' % objcls[i])

    fignum+=1

    for j in range(0,num_ext):
        if len(class_data)==1:continue
        ax=plt.subplot(2,num_ext,j+num_ext+1)
        ext_data=column(class_data,'%s' % ext_list[j])
        clean_data=[x for x in ext_data if not math.isnan(x)]#remove nan from data
        maxbin=int(math.ceil( max(column(class_data,'%s' % ext_list[j])) ))
        if maxbin==0: maxbin=1
        n, bins, patches = plt.hist(clean_data, maxbin,
                                    normed=1, facecolor='green')
        bincenters = 0.5*(bins[1:]+bins[:-1])
        chi2_1vals=chi2func(bincenters, df=1.0, norm=1.0)
        chi2_2vals=chi2func(bincenters, df=2.0, norm=1.0)
        plt.plot(bincenters,chi2_1vals,'r--',linewidth=2.0,label='$\chi^2_1$/2')
        plt.plot(bincenters,chi2_2vals,'b--',linewidth=2.0,label='$\chi^2_2$/2')
        legend = ax.legend(loc='upper right',frameon=False)
        plt.ylabel('PDF')
        plt.xlabel('TS$_{%s}$' % ext_list[j][0:4])
        plt.yscale('log')
        plt.xlim([0,maxbin+1])
        plt.ylim([0.00001,1.])
    
    if args.save_plots:
        save('figures/%s_%s' % (args.output, objcls[i]),'png',False,True)


print "You observed a total of", len(data), "objects"

if args.show_plots:
    print "Displaying plots, please close them to exit"
    plt.show()


