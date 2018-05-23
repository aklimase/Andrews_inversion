#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  9 16:42:35 2018

@author: temp


"""

import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.gridspec import GridSpec
import numpy as np
import glob
import dread
from obspy import read
import os.path as path
import seaborn as sns
#import mpl_defaults
#from mpl_defaults import plot_defaults
plt.style.use("classic")
mpl.rcParams['font.size'] = 30

mpl.rcParams['figure.subplot.hspace'] = 0.3
mpl.rcParams['figure.subplot.left'] = 0.12
mpl.rcParams['figure.subplot.right'] = 0.98
mpl.rcParams['figure.subplot.top'] = 0.98
mpl.rcParams['figure.subplot.bottom'] = 0.07


sns.set_context("poster")
sns.set(rc={'axes.labelsize':28, 'xtick.labelsize':28, 'ytick.labelsize':28, "font.size":30})
sns.set(font_scale=2.5)
sns.set_style('whitegrid')

top_dir = '/Volumes/USGS_Data/project'
box = 'all_paths'
boxpath = top_dir + '/boxes/' + box
#outfile_path = boxpath + '/secondo_rebin3'
f = '/Volumes/USGS_Data/project/catalogs/all_paths_M2.5_USGS_Catalog.txt'
data = np.genfromtxt(f, dtype = None, comments = '#', delimiter = '|', names = True)

                     
depth = data['Depthkm']
mag = data['Magnitude']

gs = GridSpec(3,1)

fig = plt.figure(figsize = (14,16))
#fig.text(0.01, 0.5,'counts', va='center', rotation='vertical')

plt.subplot(gs.new_subplotspec((0, 0), colspan=1))
plt.hist(depth, bins = np.arange(0,28, 2), histtype = 'bar')
plt.xlabel('Event Depth (km)')
plt.ylabel('Counts')
plt.xlim(0,28)
plt.xticks(np.arange(0,29, 4))
plt.yticks(np.arange(0,1000, 200))

plt.subplot(gs.new_subplotspec((1, 0), colspan=1))
plt.hist(mag, bins = np.arange(2.5,5.75, 0.25))
plt.xlabel('Event Magnitude')
plt.ylabel('Counts')
plt.xlim(2.5,5.75)
#plt.ylim(0,1700)
plt.xticks(np.arange(2.5,5.7, 0.5))
plt.yticks(np.arange(0,1900, 500))

f = '/Volumes/USGS_Data/project/catalogs/rrup_Catalog.txt'
data = np.genfromtxt(f, dtype = None, comments = '#', delimiter = '\t', names = True)

Rrup = data['Rrupkm']
m = data['mag']
depth = data['depthkm']

plt.subplot(gs.new_subplotspec((2, 0), colspan=1))
plt.hist(Rrup, bins = np.arange(0,250, 12.5))
plt.xlabel(r'Record $R_{rup}$ (km)')
plt.ylabel('Counts')
#plt.xlim(2.5,5.75)
plt.xticks(np.arange(0,250, 25))
plt.yticks(np.arange(0,8000, 2000))
plt.show()
plt.savefig('/Volumes/USGS_Data/notes/kappa_paper/histo.png')

##############################
#gs = GridSpec(1,2,width_ratios=[1, 0.05])

#ax = plt.subplot(gs.new_subplotspec((0, 0), colspan=1))
#fig = plt.subplot(figsize = (14,14))

#hexplot = sns.jointplot(x, y, kind="hex")
#sns.plt.subplots_adjust(left=0.2, right=0.8, top=0.8, bottom=0.2)  # shrink fig so cbar is visible
#cax = hexplot.fig.add_axes([.85, .25, .05, .4])  # x, y, width, height
#sns.plt.colorbar(cax=cax)

g=sns.jointplot(x = Rrup, y= m, stat_func=None, size=16)

#Clear the axes containing the scatter plot
g.ax_joint.cla()

#Generate some colors and markers
norm=plt.Normalize(0,np.max(depth))
#colors = plt.cm.plasma(norm(depth.fillna(0).values))
colors = plt.cm.plasma_r(norm(depth))

#Plot each individual point separately
#for i in range(len(Rrup)):
for i in range(500):
    g.ax_joint.scatter(Rrup[i], m[i], color=colors[i])

g.set_axis_labels(r'Record $R_{rup}$ (km)', 'Magnitude', fontsize=30)
g.ax_joint.set_xlim(0,250)

#
#cmap = plt.cm.plasma_r
#s_m = mpl.cm.ScalarMappable(cmap = cmap, norm=norm)
#s_m.set_array([])

#ax = plt.subplot(gs.new_subplotspec((0, 1), colspan=1))
# Make space for the colorbar
plt.subplots_adjust(left=0.07, right=0.98, top=0.98, bottom=0.3)  # shrink fig so cbar is visible
#cax = g.fig.add_axes([.8, .1, .03, .6])  # x, y, width, height
#plt.colorbar(s_m, cax=cax)

axHist = g.fig.add_axes([.07, .06, 0.78, 0.18])
# get hist data
N, bins, patches = axHist.hist(depth, bins=np.arange(0,24,1)) #, orientation=u'horizontal')
#plt.xlim(0, max(depth))
plt.xlabel(r'Depth (km)', fontsize = 30)
plt.yticks(np.arange(0,6000, 2000))

# set a color for every bar (patch) according 
# to bin value from normalized min-max interval
for bin, patch in zip(bins, patches):
    color = plt.cm.plasma_r(norm(bin))
    patch.set_facecolor(color)

plt.show()
plt.savefig('/Volumes/USGS_Data/notes/kappa_paper/scatter_r.png')


#cbar = plt.colorbar(s_m, cax=cbar_ax)
##########################################################

#norm=plt.Normalize(0,np.max(depth))
##colors = plt.cm.plasma(norm(depth.fillna(0).values))
#colors = plt.cm.plasma(norm(depth))
#
#h = sns.jointplot(x = Rrup, y= m, c = colors, space = 0, size = 14)
#plt.show()
#plt.savefig('/Volumes/USGS_Data/notes/kappa_paper/histo.png')

###############################################
#Rrup = []
#ev = []
#magl = []
#dep = []
#outfilename = '/Volumes/USGS_Data/project/catalogs/rrup_Catalog.txt'
#top_dir = '/Volumes/USGS_Data/project'
#box = 'all_paths'
#boxpath = top_dir + '/boxes/' + box
#record_path = glob.glob(boxpath + '/record_spectra_rebin3/Event_*/*.out')
#
#for i in range(len(record_path)):
#    #read in uncorrected data for header info
#    record = (record_path[i].split('/')[-1])
#    ev.append(record.split('.')[0])
#    base = path.basename(record)
#    network, station, channel, loc = base.split('_')[0:4]
#    yyyy, month, day, hh, mm, ss = base.split('_')[4:]
#    ss = ss.split('.')[0]
#    eventid = yyyy + '_' + month + '_' + day + '_' + hh + '_' + mm + '_' + ss
#
#    #read in uncorrected data for header info
#    raw_file = boxpath + '/uncorrected/Event_'+ eventid + '/' + network + '_' + station + '_HHN_' + loc + '_' + eventid + '.SAC'
#    stream = read(raw_file)
#    tr = stream[0]
#
#    evlon =  tr.stats.sac.evlo #deg
#    evlat =  tr.stats.sac.evla #deg
#    evdepth = tr.stats.sac.evdp #km
#    stlon = tr.stats.sac.stlo #deg
#    stlat = tr.stats.sac.stla #deg
#    stdepth = tr.stats.sac.stdp #km
#    
#    mag = tr.stats.sac.mag
#
#    #find distance between event and station
#    dist =  dread.compute_rrup(evlon, evlat, evdepth, stlon, stlat, stdepth) #in km
#    Rrup.append(dist)
#    magl.append(mag)
#    dep.append(evdepth)
#
#outfile = open(outfilename, 'w')
#out = (np.array([ev, Rrup, magl, dep])).T
#outfile.write('#Record \t Rrup(km) \t mag \t depth(km) \n')
#np.savetxt(outfile, out, fmt = '%s', delimiter='\t')
#outfile.close()
    
##    

