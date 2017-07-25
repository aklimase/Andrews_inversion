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


box = 'Imperial_Valley_PFO_TPFO_PMD'
#box = 'Imperial_Valley_SWS_ERR'
#box = 'Riverside_FRD_RDM'
#box = 'Salton_Trough_SWS_ERR'
print box

boxpath = '/Users/escuser/project/boxes/' + box

record_paths = glob.glob(boxpath + '/record_spectra/Event_*/*.out')#full path for only specified channel

event_dir = boxpath + '/secondo/'
event_files = glob.glob(event_dir + '*.out')

station_dir = boxpath + '/secondo/'
station_files = glob.glob(station_dir + '*.out')

#print record_paths

#L1_norm = L1norm(record_paths)
#L1_norm_mean = np.mean(L1_norm, axis=0)
#L1_norm_std = np.std(L1_norm, axis=0)

#make a 2d array with each row a record and each col a freq
residual = np.zeros((len(record_paths), 50))

for i in range(len(record_paths)):  ##for every record
#for i in range(5):  ##for every record
    #North component
    base = path.basename(record_paths[i])
#    eventid = base.split('.')[0]
#    network = base.split('.')[1]
#    station = base.split('.')[2]
    
#    eventid, network, station, channel, extn = base.split('.')
    
#    record = record_paths[i]
    network, station, channel, loc = base.split('_')[0:4]
    yyyy, month, day, hh, mm, ss = base.split('_')[4:]
    ss = ss.split('.')[0]
    eventid = yyyy + '_' + month + '_' + day + '_' + hh + '_' + mm + '_' + ss    
    
    #read in raw data for record info
    #correct for distance
    raw_file = boxpath + '/uncorrected/Event_'+ eventid + '/' + network + '_' + station + '_HHN_' + loc + '_' + eventid + '.SAC'
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
#    fig = plt.figure(figsize = (12,25))
#    fig.text(0.04, 0.5, 'Velocity amplitude (cm/s)', va='center', rotation='vertical', fontsize = 15)
#    plt.subplot(411)
#    plt.title('record ' + base, fontsize = 15)
#    plt.grid()
#    plt.loglog(f_bins, record_spec, color = 'b', label = 'record')
#    plt.loglog(f_bins, station_spec*event_spec, color='black', label = 'event*site')
#    plt.legend(loc=1, borderaxespad=0.)
#    plt.subplot(412)
#    plt.title('event ' +  eventid, fontsize = 15)
#    plt.grid()
#    plt.loglog(f_bins, event_spec, color = 'green')
#    plt.subplot(413)
#    plt.title('station ' + station, fontsize = 15)
#    plt.grid()
#    plt.loglog(f_bins, station_spec, color='r')
#    plt.subplot(414)
#    plt.title('Residuals', fontsize = 15)

###############################################################################
    residual[i:,] = np.log(record_spec) - np.log(station_spec*event_spec)
    
#    plt.grid()
#    plt.plot(f_bins, station_spec*event_spec - record_spec, color='black', label = 'event*site')
#    plt.xscale('log')
##    plt.xlim(0.1, 50)
#    plt.ylim(-500,500)
##    plt.loglog(f_bins, station_spec*event_spec + L1_norm_std, color='black', ls = '--', label = '1 sigma L1')
##    plt.loglog(f_bins, station_spec*event_spec - L1_norm_std, color='black', ls = '--')
#    plt.legend(loc=1, borderaxespad=0.)
#    plt.xlabel('Frequency (Hz)', fontsize = 15)
#    plt.savefig(boxpath + '/spectra_plots/' + eventid + '.' + network + '.' +  station + '.png')
#    plt.close()
##    plt.show()

print(len(residual[0]))
#print(residual[0:,])#first row?
print(len(residual[:,0]))
mean = np.mean(residual, axis = 0)
std = np.std(residual, axis = 0)

fig = plt.figure(figsize = (20,15))
title = 'Event*site - Record spectra ' + box
plt.title(title, fontsize = 20)
plt.xscale('log')
plt.ylabel('residual', fontsize = 15)
plt.xlabel('frequency (Hz)', fontsize = 15)

for i in range(len(residual[:,0])):
    plt.plot(f_bins, residual[i])
    plt.ylim(-50,50)
    plt.hold(True)
    
plt.plot(f_bins, mean, color = 'black')
plt.errorbar(f_bins, mean, yerr = std, zorder = 2000, color = 'black')
#plt.plot(f_bins, mean-std, color = 'black')
plt.savefig(boxpath + '/residuals.png')
plt.show()
    
