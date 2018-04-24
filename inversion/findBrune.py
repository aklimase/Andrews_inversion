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
import mpl_defaults
from mpl_defaults import plot_defaults
plot_defaults()

#change these magnitude upper and lower bounds
#these are for Joe's ~3 eventds
mag_ub = 3.06#2.77
mag_lb = 3.06#2.75
beta = 3500. #3500m/s
stressdrop = 5e6 #pascals
U = 0.63#0.63
rho = 2750. #kg/m^3

top_dir = '/Volumes/USGS_Data/project'
box = 'all_paths'
boxpath = top_dir + '/boxes/' + box
##############################################
secondo_dir = 'secondo_rebin3'
writefile = 'no' #or no
#############################################

event_spectra_dir = boxpath + '/' + secondo_dir+ '/'
event_spectra = glob.glob(event_spectra_dir + '[2]*.out')

#find events in catalog that are in mag range
catalog = top_dir + '/catalogs/all_paths_M2.5_USGS_Catalog.txt'
cat = np.genfromtxt(catalog, comments = '#', delimiter = '|', dtype = None, usecols = [1,10])
event = []
magl = []
for i in range(len(cat)):
    m = cat[i][1]
    if m >= mag_lb and m <= mag_ub:
        magl.append(cat[i][1])
        time = obspy.core.utcdatetime.UTCDateTime(cat[i][0].split('.')[0])
        ev = str(time.year).zfill(4) + '_' + str(time.month).zfill(2) + '_' + str(time.day).zfill(2) + '_' + str(time.hour).zfill(2) + '_' + str(time.minute).zfill(2) + '_' + str(time.second).zfill(2)
        event.append(ev)
       
#compute Brune spectra for all of the events in directory
cf_list = []
Brune_list = []
spec_list = []
ev_list = []

spec_demean_list = []
Brune_demean_list = []

spec_demean_list_log = []
Brune_demean_list_log = []
 
for i in range(len(event)):
    f = event_spectra_dir + event[i] + '.out'
    if f in event_spectra:
        ev_list.append(event[i])
        data = np.genfromtxt(f, dtype = float, comments = '#', delimiter = None, usecols = (0,1))#only read in first two cols
        freq = data[:,0]
        #these record spectra are in cm
        spec = (data[:,1])
        ml = magl[i]
        #if less than 3, convert local magnitude to moment magnitude
        if ml < 3.0:
            M = 0.884 + 0.667*ml#0.884 + 0.667*ml
        else:
            M = magl[i]
        #compute Brune in SI units
        #Moment from the moment magnitude
        M0 = 10.**((3./2.)*M + 9.1)
        #corner frequency
        fc = beta*(stressdrop/(8.47*M0))**(1./3.)
        print fc
        omega0 = (M0*U)/(4.*rho*np.pi*(beta**(3.0)))
        #brune spectra over all frequencies
        Brune = (2.*np.pi*(freq)*omega0)/(1.+((1./fc)*freq)**2.)
        #stay in meters

        cf_list.append(np.log10(spec)-np.log10(Brune))

        Brune_list.append(Brune)
        spec_list.append(spec)

        

for i in range(len(ev_list)):
    fig = plt.figure(figsize = (10,7))
    plt.ylabel('Velocity amplitude (m^2)', fontsize = 16)
    plt.xlim(0.5,70)
#    plt.ylim(-1,1)
    plt.loglog(freq , spec_list[i], color = 'green', label = 'event spectra')
#    plt.loglog(freq , spec_demean_list[i], color = 'green', ls = '--', label = 'demeaned event spectra')
    plt.grid()
    plt.loglog(freq, Brune_list[i], color = 'blue', label = 'Brune spectra')
#    plt.loglog(freq , Brune_demean_list[i],color = 'blue', ls = '--', label = 'demeaned Brune spectra')
    plt.legend(loc = 'lower right', fontsize = 16)
    plt.xlabel('Frequency (Hz)', fontsize = 16)
    plt.title(ev_list[i])
    plt.tick_params(axis='both', which='major', labelsize=15)
    plt.tick_params(axis='both', which='both', length = 5, width = 1)
    plt.text(0.7, .1, 'Avg log(diff) 0.1-36Hz: ' + str(round(np.mean(cf_list[i][20:70]),3)), fontsize = 16)
#    plt.savefig(boxpath + '/secondo_Joe/Demeaned_' + ev_list[i] + '.png')
    plt.show()


#for each event, find A/B for all other events
#sum up all A/B over freqencies we are fitting
cfarray = np.array(cf_list)
#ind = cfarray.index(0.5)
sum_list = map(sum,cfarray[:,np.arange(20,70)]**2.0) ###found the best fit from 0.5-36Hz
print sum_list
##find the minimum in log space
ind = sum_list.index(min(sum_list))
print(ev_list[ind])

#write the constraint file in linear space to agree with the event and station spectra
if writefile == 'yes':
    outfile = open(boxpath + '/constraint_rebin3_' + ev_list[ind] + '.out', 'w')
    out = (np.array([freq, (10.**(cf_list[ind]))]).T)
    outfile.write('#freq_bins \t cf_m \n')
    np.savetxt(outfile, out, fmt=['%E', '%E'], delimiter='\t')
    outfile.close()

