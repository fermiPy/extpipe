import sys
import os
import argparse
import math

import numpy as np
import matplotlib.pyplot as plt

from astropy.io import fits

import plotting

parser = argparse.ArgumentParser()

parser.add_argument('--output', default = 'test')
parser.add_argument('--input', default = 'test.fits',
                    help='input fits file from aggregate output')
parser.add_argument('--show_plots', dest='show_plots', action='store_true',
                    help='include this argument to display plots')
parser.add_argument('--save_plots', dest='save_plots', action='store_true',
                    help='include this argument to save plots')
parser.add_argument('--plot_cdf', dest='plot_cdf', action='store_true',
                    help='in the TS distribution, plot CDF instead of PDF')
parser.set_defaults(show_plots=False)
parser.set_defaults(save_plots=False)
parser.set_defaults(plot_cdf=False)

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

ext_list='ext1_ts' #['ext0_ts','ext1_ts','ext2_ts']


plt.figure(figsize=(8,8),facecolor='w',edgecolor='w')
plt.figure(1)
plotting.plot_aitoff(data,subplot=[2,2,1])

plt.subplot(2,2,2)
plotting.plot_xy(data,x='ts',y='ext2_ts',rangex=[0,1000.],rangey=[0,1000.])

plt.subplot(2,1,2)
plotting.plot_class_stats(data,objcls,out_stats=True)

if args.save_plots:
    plotting.save('figures/%s' % args.output,'png',False,True)

plt.figure(figsize=(8,4),facecolor='w',edgecolor='w')
plt.figure(2)
#this has to be its own figure... tbd
plotting.plot_ts_vs_chi2(data,ext_list='ext2_ts', ndf_chi2=[1.0,2.0])

if args.save_plots:
    plotting.save('figures/%s_ts' % args.output,'png',False,True)

plt.figure(figsize=(8,6),facecolor='w',edgecolor='w')
plotting.plot_xy(data, x='ext2_ts', y=plotting.delta_dlike(data,x='fit1_dlike',y='fit2_dlike'), mathy='yes')
    
if args.save_plots:
    plotting.save('figures/%s_ts_vs_2dll' % args.output,'png',False,True)


#fignum=3
for i in range(0,len(objcls)):
    class_data=[s for s in data if objcls[i] in s]
    if len(class_data)<2:continue
    #print objcls[i]
    plt.figure(figsize=(8,8),facecolor='w',edgecolor='w')
    #plt.figure(fignum)
    plotting.plot_aitoff(class_data,subplot=[2,2,1])
    plt.subplot(2,2,2)
    plotting.plot_xy(class_data,x='ts',y='ext2_ts',rangey=[0,10.])
    plotting.plot_ts_vs_chi2(class_data,ext_list='ext0_ts', ndf_chi2=[1.0,2.0], subplot=[2,2,3])

    #fignum+=1
    plt.suptitle('%s' % objcls[i])

    if args.save_plots:
        plotting.save('figures/%s_%s' % (args.output, objcls[i]),'png',False,True)

    plt.figure(figsize=(8,6),facecolor='w',edgecolor='w')
    plotting.plot_xy(class_data, x='ext2_ts', y=plotting.delta_dlike(class_data,x='fit1_dlike',y='fit2_dlike'), mathy='yes')
    
    if args.save_plots:
        plotting.save('figures/%s_%s_ts_vs_2dll' % (args.output, objcls[i]),'png',False,True)



print "You observed a total of", len(data), "objects"

if args.show_plots:
    print "Displaying plots, please close them to exit"
    plt.show()


