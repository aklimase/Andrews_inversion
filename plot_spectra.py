#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jul  6 15:01:34 2017

@author: escuser

inputs: reads in the record spectra, and secondo calculated event and site spectra

method: reads in record, event, and site spectra and calculates a list of the L1norm for each record 
and frequency bin, calculates the mean and standard deviation of L1norm for each frequency bin, 
corrects the record for distance, 
calculates model predicted record spectra by multiplying event and site spectra, 
plots the record spectra and the event*site spectra, separate event and site spectra, 
and the event*site spectra along with 1 sigma L1norms

outputs: write plots as .png s to spectra_plots dir
"""

import matplotlib.pyplot as plt
import numpy as np
import os.path as path
import glob
import dread
from obspy import read  
from spec_func import L1norm

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

boxpath = '/Users/escuser/Documents/Alexis_Data/cut_sac_files'

record_paths = glob.glob('/Users/escuser/Documents/Alexis_Data/cut_sac_files/record_spectra/*.*.out')#full path for only specified channel

#record_paths = glob.glob('/Users/escuser/Documents/Alexis_Data/cut_sac_files/record_spectra/10891517.*.out')#full path for only specified channel
#record_paths.extend(glob.glob('/Users/escuser/Documents/Alexis_Data/cut_sac_files/record_spectra/14598996.*.out'))#full path for only specified channel


event_dir = '/Users/escuser/Documents/Alexis_Data/cut_sac_files/event_site_spectra/'

station_dir = '/Users/escuser/Documents/Alexis_Data/cut_sac_files/event_site_spectra/'
station_files = glob.glob(station_dir + '*.out')

L1_norm = L1norm(record_paths)
L1_norm_mean = np.mean(L1_norm, axis=0)
L1_norm_std = np.std(L1_norm, axis=0)

for i in range(len(record_paths)):  ##for every record
#for i in range(5):  ##for every record
    #North component
    base = path.basename(record_paths[i])
#    eventid = base.split('.')[0]
#    network = base.split('.')[1]
#    station = base.split('.')[2]
    
    eventid, network, station, channel, extn = base.split('.')
    #read in raw data
    #correct for distance
    raw_file = boxpath + '/rawdata/' + eventid + '.' + network + '.' +  station + '.HHN.sac'
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
    
    #look up event info and print
    for j in range(len(eventid_cat)):
        if eventid_cat[j] == eventid:
            print('event: ' + eventid + '   magnitude: ' + mag_string[j])


    #record spectra
    record_data = np.genfromtxt(record_paths[i], dtype = float, comments = '#', delimiter = None, usecols = (0,1))#only read in first two cols
    f_bins = record_data[:,0]
    
    record_spec = record_data[:,1]*dist
    #event spectra
    event_data = np.genfromtxt(event_dir + eventid + '.out', dtype = float, comments = '#', delimiter = None, usecols = (0,1))
    event_spec = event_data[:,1]
    #station spectra
    station_data = np.genfromtxt(station_dir + station + '.out', dtype = float, comments = '#', delimiter = None, usecols = (0,1))
    station_spec = station_data[:,1]
    
    
    #plot the station spectra
    fig = plt.figure(figsize = (12,25))
    fig.text(0.04, 0.5, 'Velocity amplitude (cm/s)', va='center', rotation='vertical', fontsize = 15)
    plt.subplot(411)
    plt.title('record ' + base, fontsize = 15)
    plt.grid()
    plt.loglog(f_bins, record_spec, color = 'b', label = 'record')
    plt.loglog(f_bins, station_spec*event_spec, color='black', label = 'event*site')
    plt.legend(loc=1, borderaxespad=0.)
    plt.subplot(412)
    plt.title('event ' +  eventid, fontsize = 15)
    plt.grid()
    plt.loglog(f_bins, event_spec, color = 'green')
    plt.subplot(413)
    plt.title('station ' + station, fontsize = 15)
    plt.grid()
    plt.loglog(f_bins, station_spec, color='r')
    plt.subplot(414)
    plt.title('Residuals', fontsize = 15)
    plt.grid()
    plt.loglog(f_bins, station_spec*event_spec, color='black', label = 'event*site')
    plt.loglog(f_bins, station_spec*event_spec + L1_norm_std, color='black', ls = '--', label = '1 sigma L1')
    plt.loglog(f_bins, station_spec*event_spec - L1_norm_std, color='black', ls = '--')
    plt.legend(loc=1, borderaxespad=0.)
    plt.xlabel('Frequency (Hz)', fontsize = 15)
    plt.savefig(boxpath + '/spectra_plots/' + eventid + '.' + network + '.' +  station + '.png')
#    plt.show()
