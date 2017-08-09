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

#read in all instrument corrected data and compute spectra
#boxpath = '/net/anzanas.wr.usgs.gov/data/users/alexis/ANZA_boxes/Riverside_FRD_RDM'
#boxpath = '/Users/escuser/Documents/Alexis_Data/cut_sac_files'
#boxpath = '/Users/escuser/project/boxes/Riverside_FRD_RDM'
#boxpath = '/Users/escuser/project/boxes/Imperial_Valley_SWS_ERR'
#boxpath = '/Users/escuser/project/boxes/Salton_Trough_SWS_ERR'

#boxpath = '/Users/alexisklimasewski/Documents/USGS/Riverside_FRD_RDM'
boxpath = '/Users/escuser/project/boxes/all_paths'


event_dirs = glob.glob(boxpath + '/corrected/Event_*')

events = []
for i in range(len(event_dirs)):
    events.append(path.basename(event_dirs[i]))
#make event directories within corrected local data
#make a directory for each event
for i in range(len(events)):
    if not path.exists(boxpath + '/record_spectra/' + events[i]):
        os.makedirs(boxpath + '/record_spectra/' + events[i])


#read in only N, E, and return geometrical average
#loop for every event and station

#hard_start = 0
hard_start = events.index('Event_2010_05_31_05_19_24')
print(hard_start)

for i in range(hard_start, len(events)):
#    event = recordpaths_N.split('/')[-2]#event dir second to last item in path
    event = events[i][6:]
    print 'binning and fft of event: '+ event
    recordpaths = glob.glob(boxpath + '/corrected/Event_' + event +'/*_*_HHN*.SAC')#full path for only specified channel
#    recordpaths_E = glob.glob(boxpath + '/corrected/Event_' + event + '/*_*_HHE*.SAC')#full path for only specified channel
    stns = [(x.split('/')[-1]).split('_')[1] for x in recordpaths]
#    network = [(x.split('/')[-1]).split('_')[0] for x in recordpaths]
#    files = glob.glob(event_dirs[i] + '/*.SAC')
    for j in range(len(stns)):
        recordpath_E = glob.glob(boxpath + '/corrected/Event_' + event +'/*_' + stns[j] + '_HHE*.SAC')
        recordpath_N = glob.glob(boxpath + '/corrected/Event_' + event +'/*_' + stns[j] + '_HHN*.SAC')
        if(len(recordpath_E) == 1 and len(recordpath_N) == 1):
            #North component
            base_N = path.basename(recordpath_E[0])
            base_E = path.basename(recordpath_N[0])
    
        #    eventid = base_N.split('.')[0]
            network = base_N.split('_')[0]
            station = base_N.split('_')[1]
            full_channel_N = base_N.split('_')[2]
            full_channel_E = base_E.split('_')[2]
            #mtspec returns power spectra (square of spectra)
            stream = read(recordpath_N[0])
            tr = stream[0]
            data = tr.data
            ##########################################################################
            ## convert to cm/s from m/s before mtspec
            data = data*100
            
            spec_amp, freq = mtspec(data, delta = 0.01, time_bandwidth = 4, number_of_tapers=7, quadratic = True)
            spec_array_N = np.array(spec_amp)
            freq_array_N = np.array(freq)
        
            stream = read(recordpath_E[0])
            tr = stream[0]
            data = tr.data
            spec_amp, freq = mtspec(data, delta = 0.01, time_bandwidth = 4, number_of_tapers=7, quadratic = True)
            spec_array_E = np.array(spec_amp)
            freq_array_E = np.array(freq)
            
            if(len(spec_array_E)==len(spec_array_N)):
                ####here we bin into evenly spaced bins with frequency
                #spectra is power spectra so add the two components
                data_NE_2 = spec_array_E + spec_array_N
                #now data is NE power spectra
                #take the square root for normal velocity spectra
                data_NE = np.sqrt(data_NE_2)
                
                bins, binned_data = bin_spec(data_NE, freq, num_bins = 50)
                
                #make sure that all spec is a number
                if (np.isnan(binned_data).any() == False):
                ##write to file
                    outfile = open(boxpath + '/record_spectra/'+ events[i] + '/'+ network + '_' + station + '_' + 'HHNE' + '__' + event + '.out', 'w')
                    data = np.array([bins, binned_data])
                
                    data = data.T
                    outfile.write('#bins \t \t vel_spec_NE_cm \n')
                    np.savetxt(outfile, data, fmt=['%E', '%E'], delimiter='\t')
                    outfile.close()
        











