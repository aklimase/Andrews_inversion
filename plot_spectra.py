#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jul  6 15:01:34 2017

@author: escuser

Read in the .out spectra files and plot along with info for each event, and plot the record spectra for each combo
"""

import matplotlib.pyplot as plt
import numpy as np
import os.path as path
import glob

#read in anza_event_catalog.txt
catalog = '/Users/escuser/Documents/Alexis_Data/anza_event_catalog.txt'

eventid = []
mag = []

cat = np.genfromtxt(catalog, dtype = str, comments = '#', delimiter = '|')
eventid.append(cat[:,0])
eventid_string = eventid[0]
mag.append(cat[:,10])
mag_string = mag[0]

print(eventid_string)

## Convert eventid_string to a float:
#eventid_cat = np.zeros(len(eventid_string))
eventid_cat = []
#mag_cat = np.zeros(len(eventid_string))

for eventi in range(len(eventid_string)):
    eventid_i = eventid_string[eventi].split('ci')[1]
    eventid_cat.append(eventid_i)
#    mag_cat[eventi] = float(mag_string[eventi])
print(eventid_cat)

#remove the ci from the name

record_paths = glob.glob('/Users/escuser/Documents/Alexis_Data/cut_sac_files/record_spectra/*.out')#full path for only specified channel

event_dir = '/Users/escuser/Documents/Alexis_Data/cut_sac_files/event_site_spectra/'

station_dir = '/Users/escuser/Documents/Alexis_Data/cut_sac_files/event_site_spectra/'
station_files = glob.glob(station_dir + '*.out')

for i in range(15, 25):  ##for every record
    #North component
    base = path.basename(record_paths[i])
    eventid = base.split('.')[0]
    network = base.split('.')[1]
    station = base.split('.')[2]
    
    #look up event info and print
    for j in range(len(eventid_cat)):
        if eventid_cat[j] == eventid:
            print('event: '+ eventid + '   magnitude: ' + mag_string[j])

    #record spectra
    record_data = np.genfromtxt(record_paths[i], dtype = float, comments = '#', delimiter = None, usecols = (0,1))#only read in first two cols
    f_bins = record_data[:,0]
    #####################################
    #power spectra so take square root
    record_spec = np.sqrt(record_data[:,1])
    #event spectra
    event_data = np.genfromtxt(event_dir + eventid + '.out', dtype = float, comments = '#', delimiter = None, usecols = (0,1))
    event_spec = event_data[:,1]
    
    station_data = np.genfromtxt(station_dir + station + '.out', dtype = float, comments = '#', delimiter = None, usecols = (0,1))
    station_spec = station_data[:,1]
    
    
    #plot the station spectra
    plt.figure(figsize = (12,16))
    plt.subplot(311)
    plt.title('record ' + base, fontsize = 15)
    plt.loglog(f_bins, record_spec)
    plt.ylim(10**(-10), 10**(0))
    plt.subplot(312)
    plt.title('event ' +  eventid, fontsize = 15)
    plt.loglog(f_bins, event_spec, color = 'green')
    plt.ylim(10**(-10), 10**(0))
    plt.subplot(313)
    plt.title('station ' + station, fontsize = 15)
    plt.loglog(f_bins, station_spec, color='r')
    plt.ylim(10**(-10), 10**(0))
    plt.show()
