#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 25 09:56:29 2017

@author: escuser

plot output from running inversion with all the box data separately for each box
"""

import matplotlib.pyplot as plt
import numpy as np
import os.path as path
import glob
import dread
from obspy import read

boxpath = '/Users/escuser/project/boxes/'
box_list = ['Imperial_Valley_PFO_TPFO_PMD', 'Imperial_Valley_SWS_ERR', 'Riverside_FRD_RDM', 'Salton_Trough_SWS_ERR']
#box_list = ['Imperial_Valley_SWS_ERR', 'Imperial_Valley_PFO_TPFO_PMD']
residual_list = []

#read in corrected event directories to find event matches
for i in range(len(box_list)):
    record_paths = glob.glob(boxpath + box_list[i] + '/record_spectra/Event_*/*.out')#full path for only specified channel
    print(box_list[i], len(record_paths))
    out_dir = boxpath + 'secondo_all/'
    
    residual = np.zeros((len(record_paths), 50))
##################################
#    residual = np.zeros((10, 50))


    #for each event/station, get the information
    for j in range(len(record_paths)):
#    for j in range(15):
        base = path.basename(record_paths[j])
        network, station, channel, loc = base.split('_')[0:4]
        yyyy, month, day, hh, mm, ss = base.split('_')[4:]
        ss = ss.split('.')[0]
        eventid = yyyy + '_' + month + '_' + day + '_' + hh + '_' + mm + '_' + ss    
        
        #read in raw data for record info
        #correct for distance
        raw_file = boxpath + box_list[i] +  '/uncorrected/Event_'+ eventid + '/' + network + '_' + station + '_HHN_' + loc + '_' + eventid + '.SAC'
        stream = read(raw_file)
        tr = stream[0]
        
        evlon =  tr.stats.sac.evlo #deg
        evlat =  tr.stats.sac.evla #deg
        evdepth = tr.stats.sac.evdp #km
        stlon = tr.stats.sac.stlo #deg
        stlat = tr.stats.sac.stla #deg
        stdepth = tr.stats.sac.stdp #km
        
        #find distance between event and station
        dist =  dread.compute_rrup(evlon, evlat, evdepth, stlon, stlat, stdepth) #in km
        #km to cm
        dist = dist*100000
    
        #record spectra
        record_data = np.genfromtxt(record_paths[j], dtype = float, comments = '#', delimiter = None, usecols = (0,1))#only read in first two cols
        f_bins = record_data[:,0]
        
        record_spec = record_data[:,1]*dist
        #event spectra
#        print(out_dir + eventid + '.out')
        event_data = np.genfromtxt(out_dir + eventid + '.out', dtype = float, comments = '#', delimiter = None, usecols = (0,1))
        event_spec = event_data[:,1]
        #station spectra
#        print(out_dir + station + '.out')
        station_data = np.genfromtxt(out_dir + station + '.out', dtype = float, comments = '#', delimiter = None, usecols = (0,1))
        station_spec = station_data[:,1]

    
        residual[j:,] = np.log(record_spec) - np.log(station_spec*event_spec)
        residual_list.append(np.log(record_spec) - np.log(station_spec*event_spec))
    
    mean = np.mean(residual, axis = 0)
    print('residual: ', len(residual))
    std = np.std(residual, axis = 0)
    
#    residual_list.append(residual)
    
    fig = plt.figure(figsize = (20,15))
    title = 'log(Record spectra) - log(event*site) ' + box_list[i]
    plt.title(title, fontsize = 20)
    plt.xscale('log')
    plt.xlim(0.1, 50)
#    plt.ylim(-5,5)
    plt.ylabel('residual', fontsize = 15)
    plt.xlabel('frequency (Hz)', fontsize = 15)
    
    for k in range(len(residual[:,0])):
        plt.plot(f_bins, residual[k], alpha = 0.7)
        plt.hold(True)

    plt.errorbar(f_bins, mean, yerr = std, zorder = 2000, color = 'black', elinewidth=2, capsize = 10, markeredgewidth=2, fmt='o')
    plt.axhline(y=0.0, color='black', linestyle='-')
#    plt.savefig(boxpath + 'secondo_all/IV_only_residuals_' + box_list[i] + '.png')
    plt.show()

fig = plt.figure(figsize = (20,15))
title = 'log(Record spectra) - log(event*site) entire inversion'
plt.title(title, fontsize = 20)
plt.xscale('log')
plt.xlim(0.1, 50)
#plt.ylim(-5,5)
plt.ylabel('residual', fontsize = 15)
plt.xlabel('frequency (Hz)', fontsize = 15)

for m in range(len(residual_list)):
    plt.plot(f_bins, residual_list[m], alpha = 0.7)
    plt.hold(True)
    
mean = np.mean(residual_list, axis = 0)
std = np.std(residual_list, axis = 0)    
plt.errorbar(f_bins, mean, yerr = std, zorder = 4000, color = 'black', elinewidth=2, capsize = 10, markeredgewidth=2, fmt='o')
plt.axhline(y=0.0, color='black', linestyle='-')
#plt.savefig(boxpath + 'secondo_all/IV_only_entire_inversion.png')
plt.show()