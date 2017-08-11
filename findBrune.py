#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  7 09:53:15 2017

@author: escuser

compute “B” for a M3 earthquake, and then for every ~M3 earthquake
in your dataset (call it “A”), find A/B for every frequency range you inverted for.  
Then sum up A/B over all frequency ranges and find which earthquake has the minimum value for that
because then that would suggest that earthquake is closest to a brune spectrum
"""

import glob
import numpy as np
import obspy
import matplotlib.pyplot as plt


mag_ub = 2.77
mag_lb = 2.75
beta = 3500. #m/s
stressdrop = 5e6 #pascals
U = 0.63
rho = 2750. #kg/m^3



event_spectra_dir = '/Users/escuser/project/boxes/all_paths/secondo/'
event_spectra = glob.glob(event_spectra_dir + '*.out')


#find events in catalog that are in mag range
catalog = '/Users/escuser/project/all_paths_M2.5_USGS_Catalog.txt'
cat = np.genfromtxt(catalog, comments = '#', delimiter = '|', dtype = None, usecols = [1,10])

event = []
magl = []
for i in range(len(cat)):
    m = cat[i][1]
    if m > mag_lb and m < mag_ub:
        magl.append(cat[i][1])
        time = obspy.core.utcdatetime.UTCDateTime(cat[i][0].split('.')[0])
        ev = str(time.year).zfill(4) + '_' + str(time.month).zfill(2) + '_' + str(time.day).zfill(2) + '_' + str(time.hour).zfill(2) + '_' + str(time.minute).zfill(2) + '_' + str(time.second).zfill(2)
        event.append(ev)
       
#compute B for all of the events in out directory
cf_list = []
Brune_list = []
spec_list = []
ev_list = []
for i in range(len(event)):
    f = event_spectra_dir + event[i] + '.out'
    if f in event_spectra:
#        print(event[i])
        ev_list.append(event[i])
        data = np.genfromtxt(f, dtype = float, comments = '#', delimiter = None, usecols = (0,1))#only read in first two cols
        freq = data[:,0]
        spec = data[:,1]/100. #m/s
        ml = magl[i]
        M = 0.884 + 0.667*ml
#        M = magl[i]
        M0 = 10.**((3./2.)*M + 9.1)
        fc = beta*((stressdrop/(8.47*M0))**(1./3))
#        fc = 30.
        omega0 = ((M0*U)/(4*rho*np.pi*beta**(3.0)))*0.1
        Brune = (2*np.pi*(freq)*omega0)/(1.+((1./fc)*freq)**2.)
        cf_list.append(np.log(spec)-np.log(Brune))
        Brune_list.append(Brune)
        spec_list.append(spec)

for i in range(0,20):
    fig = plt.figure(figsize = (10,7))
    plt.ylabel('Velocity amplitude (m/s)', fontsize = 16)
    plt.xlim(0.5,70)
#    plt.loglog(freq, cf_list[i], label = 'correction factor')
    plt.loglog(freq , spec_list[i], label = 'event spectra')
    plt.grid()
    plt.loglog(freq, Brune_list[i], label = 'Brune spectra')
    plt.legend(loc = 'lower right', fontsize = 16)
    plt.xlabel('Frequency (Hz)', fontsize = 16)
    plt.title(ev_list[i])
    plt.tick_params(axis='both', which='major', labelsize=15)
    plt.tick_params(axis='both', which='both', length = 5, width = 1)
    plt.show()


#for each event, find A/B for all other events
#sum up all A/B
cfarray = np.array(cf_list)
#ind = cfarray.index(0.5)
sum_list = map(sum,cfarray[:,np.arange(21,50)]**2.0)
ind = sum_list.index(min(sum_list))
print(ev_list[ind])
print(cf_list[ind])


outfile = open('/Users/escuser/project/boxes/all_paths/constraint_' + ev_list[ind] + '.out', 'w')
out = (np.array([freq, cf_list[ind]])).T
outfile.write('#freq_bins \t log(cf) \n')
np.savetxt(outfile, out, fmt=['%E', '%E'], delimiter='\t')
outfile.close()
#pick the event with the smallest sum

