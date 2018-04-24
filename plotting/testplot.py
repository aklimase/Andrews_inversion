#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 12 17:54:14 2018

@author: temp
"""

import matplotlib.pyplot as plt
import numpy as np
import os.path as path
import glob
from obspy import read  
from matplotlib.gridspec import GridSpec
import dread
import matplotlib as mpl

plt.rcParams["font.family"] = "serif"
plt.rcParams['axes.axisbelow'] = True
mpl.rcParams.update({'font.size': 18})

#all in SI units
#mag_ub = 3.9#2.77
#mag_lb = 2.8#2.75
beta = 3500. #3500m/s
stressdrop = 5e6 #pascals
U = 0.63#0.63
rho = 2750. #kg/m^3



top_dir = '/Volumes/USGS_Data/project'

box = 'all_paths'

boxpath = top_dir + '/boxes/' + box

#Joe_ev = top_dir + '/catalogs/Joe_evlist.txt'
#events = np.genfromtxt(Joe_ev, dtype = None, comments = '#', delimiter = None, names = True)
#ev_list = events['evid']


#velocity time series N,E
seismo_files = glob.glob(boxpath + '/corrected/Event_*/*HHN*.SAC')


#record_files = []
##record_files_dist = []
#for ev in ev_list:
#    record_files.extend(glob.glob(boxpath + '/record_spectra_meters/Event_' + ev + '/*.out'))


record_files = glob.glob(boxpath + '/record_spectra_meters/Event_*/*.out')

event_dir = boxpath + '/secondo_meters/'
#event_files = glob.glob(boxpath + '/secondo/*.out')

station_dir = boxpath + '/secondo_meters/'
#station_files = glob.glob(boxpath + '/secondo/*.out')

#print record_paths
print(len(record_files))

hardstart = boxpath + '/record_spectra_meters/Event_2010_02_06_09_12_27/AZ_BZN_HHNE__2010_02_06_09_12_27.out'
ind = record_files.index(hardstart)
#
#ind = 0
for i in range(ind, ind+20):  ##for every record
#for i in range(ind, ind+20): 
#    #North component
    base = path.basename(record_files[i])
    box = record_files[i].split('/')[-4]
    network, station, channel, loc = base.split('_')[0:4]
    yyyy, month, day, hh, mm, ss = base.split('_')[4:]
    ss = ss.split('.')[0]
    eventid = yyyy + '_' + month + '_' + day + '_' + hh + '_' + mm + '_' + ss    
#    
#    #read in raw data for record info
#    #correct for distance
    raw_file_N = boxpath + '/corrected/Event_'+ eventid + '/' + network + '_' + station + '_HHN_' + loc + '_' + eventid + '.SAC'
    stream = read(raw_file_N)
    tr = stream[0]
    dataN = tr.data

    raw_file_E = boxpath + '/corrected/Event_'+ eventid + '/' + network + '_' + station + '_HHE_' + loc + '_' + eventid + '.SAC'
    stream = read(raw_file_E)
    tr = stream[0]
    dataE = tr.data
#    print(raw_file)
#    
    evlon =  tr.stats.sac.evlo #deg
    evlat =  tr.stats.sac.evla #deg
    evdepth = tr.stats.sac.evdp #km
    stlon = tr.stats.sac.stlo #deg
    stlat = tr.stats.sac.stla #deg
    stdepth = tr.stats.sac.stdp #km

    ml = tr.stats.sac.mag
#    
#    #find distance between event and station
    dist_km =  dread.compute_rrup(evlon, evlat, evdepth, stlon, stlat, stdepth) #in km
#    #km to m
    dist_m = dist_km*1000.

    record_data = np.genfromtxt(record_files[i], dtype = float, comments = '#', delimiter = None, usecols = (0,1))#only read in first two cols
    f_bins = record_data[:,0]

    #record spectra
    record_spec = record_data[:,1]

    #record spectra distance corrected
#    record_data_dist = np.genfromtxt(record_files_dist[i], dtype = float, comments = '#', delimiter = None, usecols = (0,1))#only read in first two cols

#    record_spec_dist = record_data_dist[:,1]
    #event spectra
    event_data = np.genfromtxt(event_dir + eventid + '.out', dtype = float, comments = '#', delimiter = None, usecols = (0,1))
    event_spec = event_data[:,1]
    #station spectra
    station_data = np.genfromtxt(station_dir + station + '.out', dtype = float, comments = '#', delimiter = None, usecols = (0,1))
    station_spec = station_data[:,1]

#    Bias_data = np.genfromtxt( boxpath + '/secondo_Bias/Bias.out', dtype = float, comments = '#', delimiter = None, usecols = (0,1))
#    Bias = Bias_data[:,1]

    ##calculate the Brune spectra
    #if less than 3, convert local magnitude to moment magnitude
    if ml < 3.0:
        M = 0.884 + 0.667*ml#0.884 + 0.667*ml
    else:
        M = ml
    #compute Brune in SI units
    #Moment from the moment magnitude
    M0 = 10.**((3./2.)*M + 9.1)
    #corner frequency
    fc = beta*(stressdrop/(8.47*M0))**(1./3.)
    omega0 = (M0*U)/(4.*rho*np.pi*(beta**(3.0)))
#   omega0 = M0
#    print 'local magnitude', ml, 'Moment magnitude', M, 'Moment ', M0, 'corner frequency (Hz)', fc, 'Omega0 ', omega0
    #brune spectra over all frequencies
    Brune = (2.*np.pi*(f_bins)*omega0)/(1.+((1./fc)*f_bins)**2.)
    #m^2 
#    Brune = Brune*(10.**4.)


#    plot the station spectra
    fig = plt.figure(figsize = (22,16))
    plt.tight_layout(pad=0.01, w_pad=0.01, h_pad=0.01)
    plt.subplots_adjust(left=0.06, right=0.99, top=0.95, bottom=0.12)




    gs = GridSpec(3,2)

    plt.subplot(gs.new_subplotspec((0, 0), colspan=2))
    plt.title('event ' +  eventid, fontsize = 15)
    plt.plot(dataE, color='black', label = 'EW')
    plt.plot(dataN, color = 'red', label = 'NS')
    plt.legend(loc=1, borderaxespad=0.)
#    fig.yaxis.get_major_formatter().set_powerlimits((0, 1))
    plt.gca().yaxis.get_major_formatter().set_powerlimits((0, 1))
#    fig.text(0.04, 0.3, 'Velocity amplitude (cm)', va='center', rotation='vertical', fontsize = 20)
#    fig.text(0.5, 0.05, 'distance: (km) ' + str(round(dist_km,1)) + ' distance: (cm) ' + str(round(dist_cm,2), fontsize = 20)

    plt.ylabel('Velocity (m/s)')

    plt.subplot(gs.new_subplotspec((1, 0), colspan=1))
    plt.xlim(0.1,70)
    plt.grid()
#    plt.loglog(f_bins, record_spec, color = 'b', label = 'record')
    plt.loglog(f_bins, record_spec, label = 'record')
    plt.tick_params(axis='both', which='major', labelsize=15)
    plt.tick_params(axis='both', which='both', length = 5, width = 1)
    plt.legend()
    fig.text(0.01, 0.4, 'Velocity amplitude (m)', va='center', rotation='vertical', fontsize = 12)


    fig.text(0.4, 0.01, 'distance: (km) ' + str(round(dist_km,1)) + ' distance: (m) ' + '{:.2e}'.format(dist_m), fontsize = 16)

    plt.subplot(gs.new_subplotspec((1, 1), colspan=1))
    plt.xlim(0.1,70)
    plt.grid()
#    plt.title('distance: (km) ' + str(round(dist_km,1)) + ' distance: (cm) ' + str(round(dist_cm,2)))

    plt.loglog(f_bins, record_spec*dist_m, color = 'b', label = 'record dist corr')
    plt.loglog(f_bins, station_spec*event_spec, color='black', label = 'event*site')
    plt.legend(loc=1, borderaxespad=0.)
#    plt.xlim(0.1, 50)
    plt.tick_params(axis='both', which='major', labelsize=15)
    plt.tick_params(axis='both', which='both', length = 5, width = 1)

    plt.subplot(gs.new_subplotspec((2, 0), colspan=1))
    plt.xlim(0.1,70)
    plt.xlabel('event ' + eventid, fontsize = 16)
    plt.grid()
    plt.loglog(f_bins, event_spec, color = 'green')
    plt.loglog(f_bins, Brune, ls = '--', color = 'green')
#    plt.xlim(0.1, 50)
    plt.tick_params(axis='both', which='major', labelsize=15)
    plt.tick_params(axis='both', which='both', length = 5, width = 1)


    plt.subplot(gs.new_subplotspec((2, 1), colspan=1))
    plt.xlim(0.1,70)
    plt.xlabel('station ' + station, fontsize = 16)
    plt.grid()
    plt.loglog(f_bins, station_spec, color='r')
#    plt.xlim(0.1, 50)
    plt.tick_params(axis='both', which='major', labelsize=15)
    plt.tick_params(axis='both', which='both', length = 5, width = 1)
#    plt.title('Frequency (Hz)', fontsize = 15)
    plt.savefig(boxpath + '/plots_unconstrained/' + eventid + '.' + network + '.' +  station + '.png')
    plt.show()
#    plt.close()
#plt.close('all')


