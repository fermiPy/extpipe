import sys
import glob
import os
import argparse
import yaml
import pprint

import numpy as np
import matplotlib.pyplot as plt

from astropy.table import Table, Column
from astropy.io import fits

hdulist=fits.open('test.fits')
data=hdulist[1].data

#list of object classes
objcls=['fsrq','bll','bcu','PSR','unkn']

#Ssave/Show plots at the end?
show_plots = True
save_plots = True

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


#plot ra, dec of objects
plt.figure(facecolor='w',edgecolor='w')
plt.figure(1)
plt.subplot(221)
plt.plot(column(data,'ra'),column(data,'dec'),'ro')
plt.ylabel('dec (deg)')
plt.xlabel('ra (deg)')

#plot glat, glon of objects
plt.subplot(222)
plt.plot(column(data,'glon'),column(data,'glat'),'go')
plt.ylabel('glat (deg)')
plt.xlabel('glon (deg)')

#plot number of objects in each class
objnum=[]
for i in range(0,len(objcls)):
    holder=[s for s in column(data,'class') if objcls[i] in s]
    #print holder, objcls[i]
    objnum.append(len(holder))

class_info = dict(zip(objcls, objnum))
plt.subplot(223)
plt.bar(range(len(class_info)), class_info.values(),align='center')
plt.xticks(range(len(class_info)),class_info.keys())
plt.xlabel('class')

#plot eflux at 1 GeV vs. dfde index at 1 GeV
plt.subplot(224)
plt.plot(column(data,'eflux1000'),column(data,'dfde1000_index'),'bo')
plt.ylabel('dfde 1GeV index')
plt.xlabel('eflux 1GeV')
plt.xscale('log')

if save_plots:
    save('figures/test','png',False,True)

#loop over and make figures for each specific class of objects
for i in range(0,len(objcls)):
    class_data=[s for s in data if objcls[i] in s]

    plt.figure(figsize=(10,5),facecolor='w',edgecolor='w')
    plt.figure(i+2)
    #plot glat, glon of objects
    plt.subplot(121)
    plt.plot(column(class_data,'glon'),column(class_data,'glat'),'go')
    plt.ylabel('glat (deg)')
    plt.xlabel('glon (deg)')

    #plot eflux vs dfde index
    plt.subplot(122)
    plt.plot(column(class_data,'eflux1000'),column(class_data,'dfde1000_index'),'bo')
    plt.ylabel('dfde 1GeV index for %s' % objcls[i])
    plt.xlabel('eflux 1GeV for %s' % objcls[i])
    plt.xscale('log')

    if save_plots:
        save('figures/test_%s' % objcls[i],'png',False,True)


print "You observed a total of", len(data), "objects"

if show_plots:
    plt.show()


