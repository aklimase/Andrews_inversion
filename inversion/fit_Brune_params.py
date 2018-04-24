#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 16 13:53:33 2018

@author: temp
"""

import glob
import numpy as np
import obspy
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy.optimize import curve_fit

#plt.rcParams["font.family"] = "serif"
plt.rcParams['axes.axisbelow'] = True
mpl.rcParams.update({'font.size': 18})

#change these magnitude upper and lower bounds
#these are for Joe's ~3 eventds

beta = 3500. #3500m/s
U = 0.63#0.63
rho = 2750. #2750 kg/m^3

#read in all event spectra files
top_dir = '/Volumes/USGS_Data/project'
box = 'all_paths'
boxpath = top_dir + '/boxes/' + box
event_spectra_dir = boxpath + '/secondo_constrained_2010_05_25_19_49_51/'
#event_spectra_dir = boxpath + '/secondo_meters/'

event_spectra = glob.glob(event_spectra_dir + '[2]*.out')


#find events in catalog that are in mag range
catalog = top_dir + '/catalogs/all_paths_M2.5_USGS_Catalog.txt'
cat = np.genfromtxt(catalog, comments = '#', delimiter = '|', dtype = None, usecols = [1,10])
event = []
magl = []
for i in range(len(cat)):
    m = cat[i][1]
    magl.append(cat[i][1])
    time = obspy.core.utcdatetime.UTCDateTime(cat[i][0].split('.')[0])
    ev = str(time.year).zfill(4) + '_' + str(time.month).zfill(2) + '_' + str(time.day).zfill(2) + '_' + str(time.hour).zfill(2) + '_' + str(time.minute).zfill(2) + '_' + str(time.second).zfill(2)
    event.append(ev)


#define Brune velocity spectra
def Brune(f, fc, omega0):
     return (2.*np.pi*f*omega0)/(1+(f/fc)**2.)

#for each event, use the curve fit to find the moment and stress drop from corner frequency and omega0
#for i in range(len(event)):
for i in range(3050, 3075):
    f = event_spectra_dir + event[i] + '.out'
#    f = event_spectra[i]
    data = np.genfromtxt(f, dtype = float, comments = '#', delimiter = None, usecols = (0,1))#only read in first two cols
    freq = data[:,0][20:45]#try 30 = 2.2Hz
    spec = data[:,1][20:45]
    
    popt, pcov = curve_fit(Brune, freq, spec, bounds=(0, [30., 1000.]))
    fc_fit, omega0_fit = popt
    
    M0 = (omega0_fit*4.*np.pi*rho*beta**3.)/U
    stressdrop = ((fc_fit/beta)**3.)*(8.47*M0)
    
    M = (2./3.)*(np.log10(M0) - 9.05)
    if M > 3:
        ml = M
    else:
        ml = (M - 0.884)/.667
    
    
    fig = plt.figure(figsize = (10,7))
    plt.ylabel('Velocity amplitude (m)', fontsize = 16)
    plt.xlim(0.01,70)
    plt.ylim(.001,10)
    plt.loglog(data[:,0] , data[:,1], color = 'green', label = 'event spectra')
    plt.grid()
    plt.loglog(freq, Brune(freq, fc_fit, omega0_fit), color = 'blue', label = 'Brune spectra')
    plt.legend(loc = 'upper left', fontsize = 16)
    plt.xlabel('Frequency (Hz)', fontsize = 16)
    plt.title(event_spectra[i].split('/')[7].split('.')[0])
    plt.tick_params(axis='both', which='major', labelsize=15)
    plt.tick_params(axis='both', which='both', length = 5, width = 1)
    plt.text(0.1, .005, 'ml cat: ' + "{:.3f}".format(magl[i]) + ' ml calc: ' + "{:.3f}".format(ml), fontsize = 16)
    plt.text(0.1, .003, 'M0: ' + "{:.3e}".format(M0) + ' stress drop: ' + "{:.3f}".format(stressdrop/(10.**6.)) + ' MPa', fontsize = 16)
    plt.text(0.1, .002, 'fc: ' + "{:.3f}".format(fc_fit) + ' omega0: ' + "{:.3f}".format(omega0_fit), fontsize = 16)

#    plt.savefig(boxpath + '/secondo_Joe_test/Powerspec_distancecorrect_spectra_' + ev_list[i] + '.png')
    plt.show()

    print 'ml cat: ' + "{:.3f}".format(magl[i]) + ' ml calc: ' + "{:.3f}".format(ml)
    print 'fc fit: ' + "{:.3f}".format(fc_fit)

     
#calculate ml and compare to actual ml