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

#read in anza_event_catalog.txt
catalog = '/Users/escuser/Documents/Alexis_Data/anza_event_catalog.txt'
#make a dictionary for each catalog
data = np.genfromtxt(catalog, dtype = str, comments = '#', delimiter = '|')
print(data.shape)

d = zip(data)
print(d)


record_paths = '/Users/escuser/Documents/Alexis_Data/cut_sac_files/record_spectra'
event_paths = '/Users/escuser/Documents/Alexis_Data/cut_sac_files/event_site_spectra'
station_paths = '/Users/escuser/Documents/Alexis_Data/cut_sac_files/event_site_spectra'

#for i in range(len(record_paths)):  ##for every record
#    #North component
#    base = path.basename(record_paths[i])
#    eventid = base.split('.')[0]
#    network = base.split('.')[1]
#    station = base.split('.')[2]
#    
#    #look up event info and print
#    
#    
#    #plot the record spectra
#    
#    #plot the event spectra
#    
#    #plot the station spectra
#
#
#
#plt.figure(figsize = (10,8))
#plt.loglog(freq_list, amp, color='cornflowerblue')
#plt.title('station: ' + stationlist[i])
#plt.xlabel('frequency (Hz)')
#plt.ylabel('velocity spectrum')
#plt.show()