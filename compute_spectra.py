#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 16 10:06:20 2017

@author: escuser

Should read in the instrument corrected data and write spectra files with log space binned frequencies
square the velocity spectrum
note: not the log!
"""

import matplotlib.pyplot as plt
plt.style.use("ggplot")
from obspy import read
from mtspec import mtspec
import os.path as path
import glob
import numpy as np
from library import bin

#read in all instrument corrected data and compute spectra
#boxpath = '/net/anzanas.wr.usgs.gov/data/users/alexis/ANZA_boxes/Riverside_FRD_RDM'
boxpath = '/Users/escuser/Documents/Alexis_Data/cut_sac_files'
recordpaths_N = glob.glob(boxpath + '/icorrdata/*.AZ.*.HHN.sac')#full path for only specified channel
recordpaths_E = glob.glob(boxpath + '/icorrdata/*.AZ.*.HHE.sac')#full path for only specified channel
print(len(recordpaths_N), len(recordpaths_E))

#read in only N, E, and return geometrical average
#loop for every event and station

for i in range(len(recordpaths_N)):  ##for every record
    #North component
    base_N = path.basename(recordpaths_N[i])
    base_E = path.basename(recordpaths_E[i])
    eventid = base_N.split('.')[0]
    network = base_N.split('.')[1]
    station = base_N.split('.')[2]
    full_channel_N = base_N.split('.')[3]
    full_channel_E = base_E.split('.')[3]
    
    print(eventid, network, station)
    
    #mtspec returns power spectra (square of spectra)
    stream = read(recordpaths_N[i])
    tr = stream[0]
    data = tr.data
    spec_amp, freq = mtspec(data, delta = 0.01, time_bandwidth = 4, number_of_tapers=7, quadratic = True)#nfft for number of zeros to pad with
    spec_array_N = np.array(spec_amp)
    freq_array_N = np.array(freq)

    stream = read(recordpaths_E[i])
    tr = stream[0]
    data = tr.data
    spec_amp, freq = mtspec(data, delta = 0.01, time_bandwidth = 4, number_of_tapers=7, quadratic = True)#nfft for number of zeros to pad with
    spec_array_E = np.array(spec_amp)
    freq_array_E = np.array(freq)
    
    ####here we bin into evenly spaced bins with frequency
    #spectra is power spectra so add the two components
    bins, binned_data = bin((spec_array_E + spec_array_N), freq, num_bins = 50)
    
    ##write to file
    outfile = open(boxpath + '/record_spectra/'+ eventid + '.' + network + '.' + station + '.' + 'HHNE' + '.out', 'w')
#    data = np.array([freq_array_N, np.sqrt(spec_array_E*spec_array_N), spec_array_N, spec_array_E])
    data = np.array([bins, binned_data])

    data = data.T
#    outfile.write('#freq \t \t spec_amp_NE \t spec_amp_N \t spec_amp_E  \n')
#    np.savetxt(outfile, data, fmt=['%E', '%E', '%E', '%E'], delimiter='\t')
    outfile.write('#bins \t \t binned_NE_power\n')
    np.savetxt(outfile, data, fmt=['%E', '%E'], delimiter='\t')
    outfile.close()












