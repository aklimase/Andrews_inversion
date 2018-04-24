#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  9 19:13:51 2017

@author: alexisklimasewski
"""

import numpy as np
import glob
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.patches as mpatches
import os.path as path
from obspy import read  
import dread
from pyproj import Geod

g = Geod(ellps='clrk66')

mpl.rcParams.update({'font.size': 22})
plt.rcParams["font.family"] = "serif"

#indices of frequency limits for inversion

boxnum = 4

freq_list=[(38,46),(43,49),(39,48),(38,46), (40,46)]
f1 = freq_list[boxnum][0]
f2 = freq_list[boxnum][1]

top_dir = '/Volumes/USGS_Data/project'

boxpath = top_dir + '/boxes/'
box_list = [('Imperial_Valley', 'PFO_TPFO_PMD'), ('Imperial_Valley', 'SWS_ERR'), ('Riverside', 'FRD_RDM'), ('Salton_Trough', 'SWS_ERR'), ('all_paths', '')]
box = box_list[boxnum][0] + '' + box_list[boxnum][1]
residual_dir =  top_dir + '/boxes/all_paths/residual_spectra/' + box + '/'

print box
#list of files in the directory
record_list = glob.glob(residual_dir + '*.out')

data = np.genfromtxt(record_list[0], comments = '#', dtype = float)
freq = data.T[0][f1:f2]
  
I = 1
J = 50

residual = np.zeros((len(record_list), J))
decay_fit = np.zeros((I, f2-f1))

G = np.zeros((J, 2))
d = np.zeros((J, 1))

fig = plt.figure(figsize = (10,10))
plt.title(r'$\mathrm{log(record)-log(station*event)}$')
plt.xscale('log')
#plt.yscale('log')
plt.xlim(0.1, 50)
plt.ylim(-4,4)
plt.ylabel(r'log residual')
plt.xlabel(r'$\mathrm{frequency(Hz)}$')
plt.tick_params(axis='both', which='major', labelsize=15)
plt.tick_params(axis='both', which='both', length = 8, width = 2)

C = np.zeros(len(record_list))
tstar = []
L2norm = []

for i in range(len(record_list)):#in lin space
    print record_list[i]
    data = np.genfromtxt(record_list[i], comments = '#', dtype = float)
    freq = data.T[0][f1:f2]
    Amplitude = data.T[1][f1:f2] 
    residual[i:,] = data.T[1]
    plt.hold(True)
    
    A0 = np.max(Amplitude)#lin
    C[i] = np.max(Amplitude)#lin
    plt.plot(data.T[0], np.log(data.T[1]), color = 'green', alpha= 0.2)#y = linear A
    plt.hold(True)
    
#    for j in range(len(freq)):
#        d[j][0] = np.log(Amplitude[j]/A0)#lin to log
#        G[j][0] = -1*np.pi*(freq[j])
#
#    G_inv = np.linalg.pinv(G, rcond=1e-10)
#    m = float(np.dot(G_inv,d))
#    tstar.append(m)
#
#    t = m
#    decay = C[i]*np.exp(-np.pi*t*(freq))
#    decay_fit[i:,] = decay
#    plt.plot(freq, np.log(decay), color = 'blue', alpha =0.2)
#    plt.hold(True)
    
    for j in range(len(freq)):
        d[j][0] = np.log(Amplitude[j])#lin to log
        G[j][0] = -1*np.pi*(freq[j])
        G[j][1] = 1

    G_inv = np.linalg.pinv(G, rcond=1e-10)
    m = np.dot(G_inv,d)

    t = m[0]
    A0 = m[1]
    tstar.append(float(t))
    decay = np.exp(A0)*np.exp(-np.pi*t*(freq))
    decay_fit[i:,] = decay
    plt.plot(freq, np.log(decay), color = 'blue', alpha =0.2)
    plt.hold(True)
    
    
    
    L2 = np.mean(np.sqrt((decay-Amplitude)**2)) #avg for all frequencies
    L2norm.append(L2)

print'frequency range: ', (freq[0], freq[-1])
print'L2 norm average: ', np.mean(L2norm)

##############plot residual mean on top
## convert back to log space for the mean
mean = np.mean(np.log(residual), axis = 0)
std = np.std(np.log(residual), axis = 0)

mean_lin = np.exp(mean)[f1:f2]

for j in range(len(freq)):
    d[j][0] = np.log(mean_lin[j])#lin to log
    G[j][0] = -1*np.pi*(freq[j])
    G[j][1] = 1

G_inv = np.linalg.pinv(G, rcond=1e-10)
m = np.dot(G_inv,d)
tstar_avg = m[0]
A0 = m[1]

decay_avg = np.exp(A0)*np.exp(-np.pi*tstar_avg*(freq))
print('t* of average: ' , tstar_avg)
plt.plot(freq, np.log(decay_avg), color = 'red')
plt.errorbar(data.T[0], mean, yerr = std, zorder = 2000, color = 'black', elinewidth=1, capsize = 5, markeredgewidth=1, fmt='o')
plt.axhline(y=0.0, color='black', linestyle='-')

green_patch = mpatches.Patch(color='green', label='residual')
blue_patch = mpatches.Patch(color='blue', label='t* decay')
black_patch = mpatches.Patch(color='black', label='residual average')
red_patch = mpatches.Patch(color='red', label='t* decay of average')
plt.legend(handles=[green_patch, blue_patch, black_patch, red_patch], prop={'size': 16})

plt.show()
plt.savefig(top_dir + '/boxes/all_paths/t*plots/t*_' + box + '.png')

################ tstar histogram
avg = np.mean(tstar)
sig = np.std(tstar)
med = np.median(tstar)

print 'mean :', np.round(avg*1000,5)
print 'median: ', np.round(med*1000,5)
print 'std: ', np.round(sig*1000,5)


fig = plt.figure(figsize = (20,15))
# the histogram of the data
plt.hist(tstar, 12, facecolor='green', alpha=0.75)

plt.xlabel(r'$\mathrm{t^*\ value}$')
plt.ylabel(r'$\mathrm{Number\ of\ Paths}$')
plt.title('Histogram of t* ' + r'$\mathrm{\mu}$' + ' =' + str(round(avg*1000, 5)) + ' ms')
plt.grid(True)

plt.show()
plt.savefig(top_dir + '/boxes/all_paths/t*plots/t*_histo_' + box + '.png')

#############################################################
#write output file

eventid_list = []
network_station_list = []
d_list = []
az_list = []
mag_list = []
dep_list = []
base_list = []
for i in range(len(record_list)):
    record = (record_list[i].split('/')[-1])
    base = path.basename(record)
    network, station, channel, loc = base.split('_')[0:4]
    yyyy, month, day, hh, mm, ss = base.split('_')[4:]
    ss = ss.split('.')[0]
    eventid = yyyy + '_' + month + '_' + day + '_' + hh + '_' + mm + '_' + ss
    #read in uncorrected data for header info
    raw_file = boxpath + 'all_paths/uncorrected/Event_'+ eventid + '/' + network + '_' + station + '_HHN_' + loc + '_' + eventid + '.SAC'
    stream = read(raw_file)
    tr = stream[0]
    
    base_list.append(str(base.split('.')[0]))
    
    evlon =  tr.stats.sac.evlo #deg
    evlat =  tr.stats.sac.evla #deg
    evdepth = tr.stats.sac.evdp #km
    stlon = tr.stats.sac.stlo #deg
    stlat = tr.stats.sac.stla #deg
    stdepth = tr.stats.sac.stdp #km
    mag = tr.stats.sac.mag
    
    eventid_list.append(str(eventid))
    network_station_list.append(str(network + '_' + station))
    
    az12,az21,dist = g.inv(evlon,evlat,stlon,stlat)
    mag_list.append(float(mag))
    dep_list.append(float(evdepth))
    az_list.append(float(az21))

    
    #find distance between event and station
    dist =  dread.compute_rrup(evlon, evlat, evdepth, stlon, stlat, stdepth) #in km
    d_list.append(dist)
    #km to cm



outfile = open(top_dir + '/boxes/all_paths/t*plots/' + box + '_tstar.out', 'w')
out = (np.array([base_list, eventid_list, network_station_list, tstar, d_list, az_list, mag_list, dep_list]).T)
outfile.write('#record_name \t eventid \t network_station \t tstar  \t distance_km \t azimuth \t magnitude \t depth_km \n')
np.savetxt(outfile, out, fmt='%s', delimiter='\t')
outfile.close()


#fmt=['%s','%s', '%s', '%.3f', '%.3f', '%.3f', '%.3f', '%.3f']

