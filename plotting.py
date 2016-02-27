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
from scipy.special import gammainc
from scipy.stats import chi2


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

def plot_xy(data,x='eflux1000',y='dfde1000_index',rangex=None, rangey=None, logx=False, logy=False, **kwargs):
    plt.plot(column(data,x),column(data,y),'bo')
    if rangex is not None: 
        plt.xlim([rangex[0],rangex[1]])
    if rangey is not None:
        plt.ylim([rangey[0],rangey[1]])
    plt.xlabel(x)
    plt.ylabel(y)
    if logx: plt.xscale('log')
    if logy: plt.yscale('log')


def plot_aitoff(data, coords=['glon','glat'], frame='galactic', subplot=[1,1,1], **kwargs):
    #plot coords of objects
    plt.subplot(subplot[0],subplot[1],subplot[2],projection="aitoff")
    plt.grid(True)
    c = SkyCoord(column(data,coords[0])*u.degree, column(data,coords[1])*u.degree, frame=frame)

    if frame=='galactic':
        plt.plot(c.l.wrap_at(180*u.deg).radian,c.b.wrap_at(180*u.deg).radian,'ro')
        plt.ylabel('glat (deg)')
    elif frame=='icrs':
        plt.plot(c.ra.wrap_at(180*u.deg).radian,c.dec.wrap_at(180*u.deg).radian,'ro')
        plt.ylabel('dec (deg)')
    else:
        print "currently only galactic and icrs frames are supported"


def plot_class_stats(data,objcls,out_stats=True,**kwargs):
    #plot number of objects in each object class (objcls)
    objnum=[]
    #print column(data,'class')
    for i in range(0,len(objcls)):
        holder=[s for s in column(data,'class') if objcls[i] in s]
        if out_stats: print objcls[i], len(holder)
        objnum.append(len(holder))
        
    class_info = OrderedDict(zip(objcls, objnum))
    plt.bar(range(len(class_info)), class_info.values(),align='center',log=True)
    plt.xticks(range(len(class_info)),class_info.keys())
    plt.xlabel('class')
    locs, labels = plt.xticks()
    plt.setp(labels, rotation=90)
    plt.xlim([-1.,len(class_info)])
    plt.ylim([0.1,2000.])


def plot_ts_vs_chi2(data,ext_list='ext1_ts', ndf_chi2=[1], subplot=[1,2,1], **kwargs):
    ax=plt.subplot(subplot[0],subplot[1],subplot[2])
    ext_data=column(data,'%s' % ext_list)
    clean_data=[x for x in ext_data if not math.isnan(x)]#remove nan from data
    n, bins, patches = plt.hist(clean_data, 
                                int(math.ceil( max(column(data,'%s' % ext_list)) )),
                                normed=1, facecolor='green')
    bincenters = 0.5*(bins[1:]+bins[:-1])
    chi2_vals=[]
    colors=['r','b','g']
    for j in range(0,len(ndf_chi2)):
        chi2_vals.append(chi2.pdf(bincenters, ndf_chi2[j]))
        plt.plot(bincenters,chi2_vals[j],'%s--' % colors[j],linewidth=2.0,label='$\chi^2_%i$/2' % ndf_chi2[j])
    legend = ax.legend(loc='upper right',frameon=False)
    plt.ylabel('PDF')
    plt.xlabel('TS$_{%s}$' % ext_list[0:4])
    plt.yscale('log')
    plt.ylim([0.00001,2.])
    
    ax=plt.subplot(subplot[0],subplot[1],subplot[2]+1)
    n, bins, patches = plt.hist(clean_data, 
                                int(math.ceil( max(column(data,'%s' % ext_list)) )),
                                normed=1, facecolor='green',cumulative=-1)
    chi2_sfvals=[]
    for j in range(0,len(ndf_chi2)):
        chi2_sfvals.append(chi2.sf(bincenters, ndf_chi2[j]))
        plt.plot(bincenters,chi2_sfvals[j],'%s--' % colors[j],linewidth=2.0,label='$\chi^2_%i$/2' % ndf_chi2[j])
    legend = ax.legend(loc='upper right',frameon=False)
    plt.ylabel('1-CDF')
    plt.xlabel('TS$_{%s}$' % ext_list[0:4])
    plt.yscale('log')
    plt.ylim([0.00001,2.])
    
