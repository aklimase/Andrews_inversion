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

for i in range(len(record_paths)):  ##for every record
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
    record_spec = np.sqrt(record_data[:,1])
    #event spectra
    event_data = np.genfromtxt(event_dir + eventid + '.out', dtype = float, comments = '#', delimiter = None, usecols = (0,1))
    event_spec = event_data[:,1]
    
    
#    #plot the station spectra
#    f, (ax1, ax2, ax3) = plt.subplots(3, sharex=True, sharey=True)
#    ax1.plot(x, y)
#    ax1.set_title('Sharing both axes')
#    ax2.scatter(x, y)
#    ax3.scatter(x, 2 * y ** 2 - 1, color='r')
## Fine-tune figure; make subplots close to each other and hide x ticks for
## all but bottom plot.
#f.subplots_adjust(hspace=0)
#plt.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False)






#plt.figure(figsize = (10,8))
#plt.loglog(freq_list, amp, color='cornflowerblue')
#plt.title('station: ' + stationlist[i])
#plt.xlabel('frequency (Hz)')
#plt.ylabel('velocity spectrum')
#plt.show()