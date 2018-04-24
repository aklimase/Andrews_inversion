#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 16 10:06:20 2017

@author: escuser

inputs: reads in paths to the instrument corrected HHN and HHE channel sac files

method: reads in trace data from sac files, converts m/s to cm/s to agree with secondo,
calls bin_spec function to return evenly spaced log bins and binned data,
takes the average of the N and E components

outputs: writes bins and binned spectra into the record_spectra directory
"""

import matplotlib.pyplot as plt
plt.style.use("ggplot")
from obspy import read
from mtspec import mtspec
import os
import os.path as path
import glob
import numpy as np
from spec_func import bin_spec
#import time
#import dread

top_dir = '/Volumes/USGS_Data/project'
box = 'all_paths'

boxpath = top_dir + '/boxes/' + box
event_dirs = glob.glob(boxpath + '/corrected/Event_*')

outpath = 'record_spectra_rebin3'

events = []
for i in range(len(event_dirs)):
    events.append(path.basename(event_dirs[i]))
##make event directories within corrected local data
##make a directory for each event
for i in range(len(events)):
    if not path.exists(boxpath + '/' + outpath + '/' + events[i]):
        os.makedirs(boxpath + '/' + outpath + '/'  + events[i])


for i in range(3000, len(event_dirs)):
#for i in range(1,2):
#    t1 = time.time()
    event = events[i][6:]
#    print 'binning and fft of event: '+ event
    recordpaths = glob.glob(boxpath + '/corrected/Event_' + event +'/*_*_HHN*.SAC')#full path for only specified channel
#    print recordpaths
    stns = [(x.split('/')[-1]).split('_')[1] for x in recordpaths]
    for j in range(len(stns)):
        recordpath_E = glob.glob(boxpath + '/corrected/Event_' + event +'/*_' + stns[j] + '_HHE*.SAC')
        recordpath_N = glob.glob(boxpath + '/corrected/Event_' + event +'/*_' + stns[j] + '_HHN*.SAC')
        if(len(recordpath_E) == 1 and len(recordpath_N) == 1):
            #North component
            base_N = path.basename(recordpath_E[0])
            base_E = path.basename(recordpath_N[0])
    
            network = base_N.split('_')[0]
            station = base_N.split('_')[1]
            full_channel_N = base_N.split('_')[2]
            full_channel_E = base_E.split('_')[2]
            #mtspec returns power spectra (square of spectra)
            stream = read(recordpath_N[0])
            tr = stream[0]
            data = tr.data

            ##########################################################################
            ## m/s
            data = data
            
            spec_amp, freq = mtspec(data, delta = 0.01, time_bandwidth = 4, number_of_tapers=7, quadratic = True)
            #power spectra
            spec_array_N = np.array(spec_amp)
            freq_array_N = np.array(freq)
        
            stream = read(recordpath_E[0])
            tr = stream[0]
            data = tr.data
            #m 
            data = data
            
            spec_amp, freq = mtspec(data, delta = 0.01, time_bandwidth = 4, number_of_tapers=7, quadratic = True)
            #power spectra
            spec_array_E = np.array(spec_amp)
            freq_array_E = np.array(freq)

            
            if(len(spec_array_E)==len(spec_array_N)) and len(spec_array_E)>1300:
                ####here we bin into evenly spaced bins with frequency
                #spectra is power spectra so add the two components
                data_NE_2 = spec_array_E + spec_array_N
                #now data is NE power spectra
                #take the square root for normal velocity spectra
                data_NE = np.sqrt(data_NE_2)
                
                #0.1-end
                bins, binned_data = bin_spec(data_NE[6:-1], freq[6:-1], num_bins = 75)
                
#                
#                #if less than 0.5, put into 5 bins
#                bins1, binned_data1 = bin_spec(data_NE[0:32], freq[0:32], num_bins = 5)
#                #if between 0.5 and 15.0 put into 25 bins
#                bins2, binned_data2 = bin_spec(data_NE[32:945], freq[32:945], num_bins = 25)
#                #if between 15.0 and 50 put into 30 bins
#                bins3, binned_data3  = bin_spec(data_NE[945:,], freq[945:,], num_bins = 30)
#                
#
#                bins = np.concatenate((bins1, bins2, bins3))
#                binned_data = np.concatenate((binned_data1, binned_data2, binned_data3))
                
                
                #make sure that all spec is a number
                if (np.isnan(binned_data).any() == False):
                ##write to file
                    outfile = open(boxpath + '/' + outpath + '/Event_'+ event + '/'+ network + '_' + station + '_' + 'HHNE' + '__' + event + '.out', 'w')
                    data = np.array([bins, binned_data])
                
                    data = data.T
                    outfile.write('#bins \t \t vel_spec_NE_m \n')
                    np.savetxt(outfile, data, fmt=['%E', '%E'], delimiter='\t')
                    outfile.close()
#    t2 = time.time()
#    print 'time for event: (s)', (t2-t1)
        



