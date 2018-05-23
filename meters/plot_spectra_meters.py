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
import obspy
import glob
from obspy import read  
from matplotlib.gridspec import GridSpec
import dread
import matplotlib.ticker as mtick
import waveforms as wf

import matplotlib as mpl

plt.style.use("classic")
plt.rcParams['axes.axisbelow'] = True
mpl.rcParams.update({'font.size': 20})

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

#velocity time series N,E
seismo_files = glob.glob(boxpath + '/corrected/Event_*/*HHN*.SAC')

record_files = glob.glob(boxpath + '/record_spectra_rebin3/Event_*/*.out')
        
event_dir = boxpath + '/secondo_rebin3/'

station_dir = boxpath + '/secondo_rebin3/'

print(len(record_files))

hardstart = boxpath + '/record_spectra_rebin3/Event_2013_11_30_11_36_35/AZ_BZN_HHNE__2013_11_30_11_36_35.out'
ind = record_files.index(hardstart)
#ind = 0
#
for i in range(ind, ind+1):  ##for every record
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
    raw_file_N = boxpath + '/cutdata_s/Event_'+ eventid + '/' + network + '_' + station + '_HHN_' + loc + '_' + eventid + '.SAC'
    
    prefilt = (0,0.4,35,50)
    respfile =  top_dir + '/boxes/' + box + '/respfiles/' + network + '.' + station + '.HH*.resp'
    tsunit = 'VEL'
    icorr_sacfile = top_dir + '/boxes/' + box + '/icorrN.SAC'
    
    stream = read(raw_file_N)
    tr = stream[0]
#    tr.detrend(type = 'demean')#removes mean of data
#    tr.detrend(type = 'simple')#rtrend linear from first and last samples
    #rewrite to a sac file
    tr.write('temp.sac', format = 'sac')
    sacfile = 'temp.sac'
    
    wf.remove_response(sacfile, respfile, icorr_sacfile,prefilt,tsunit)
    
    stream = read(icorr_sacfile)
    tr = stream[0]
    dataN = tr.data

    
    raw_file_E = boxpath + '/cutdata_s/Event_'+ eventid + '/' + network + '_' + station + '_HHE_' + loc + '_' + eventid + '.SAC'
    
    respfile =  top_dir + '/boxes/' + box + '/respfiles/AZ.BZN.HH*.resp'
    tsunit = 'VEL'
    icorr_sacfile = top_dir + '/boxes/' + box + '/icorrE.SAC'
    stream = read(raw_file_E)
    tr = stream[0]
#    tr.detrend(type = 'demean')#removes mean of data
#    tr.detrend(type = 'simple')#rtrend linear from first and last samples
    #rewrite to a sac file
    tr.write('temp.sac', format = 'sac')
    sacfile = 'temp.sac'
    
    wf.remove_response(raw_file_E, respfile,icorr_sacfile,prefilt,tsunit)
    
    stream = read(icorr_sacfile)
    tr = stream[0]
    dataE = tr.data
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

    #event spectra
    event_data = np.genfromtxt(event_dir + eventid + '.out', dtype = float, comments = '#', delimiter = None, usecols = (0,1))
    event_spec = event_data[:,1]
    #station spectra
    station_data = np.genfromtxt(station_dir + station + '.out', dtype = float, comments = '#', delimiter = None, usecols = (0,1))
    station_spec = station_data[:,1]
    
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

    Brune = (2.*np.pi*(f_bins)*omega0)/(1.+((1./fc)*f_bins)**2.)
    shift = np.mean(np.log10(event_spec[27:70])-np.log10(Brune[27:70]))
#    Brune = 10.**(np.log10(Brune)+shift)

#    plot the station spectra
    fig = plt.figure(figsize = (22,14))
    plt.suptitle('event ' +  eventid + ', magnitude ' + str(ml) + ', recorded on station ' + station + ', distance ' + str(round(dist_km,2)) + ' km ')

    plt.tight_layout(pad=0.01, w_pad=0.01, h_pad=0.01)
    plt.subplots_adjust(left=0.06, right=0.99, top=0.92, bottom=0.05)
    

    gs = GridSpec(3,2)
    
    plt.subplot(gs.new_subplotspec((0, 0), colspan=1))
#    plt.title('event ' +  eventid + ', magnitude ' + str(ml) + ', recorded on station ' + station + ', distance ' + str(round(dist_km,2)) + ' km ')
    x = np.arange(0, len(dataE), 1)/100.
    plt.plot(x, dataE, color='black', label = 'EW')
#    plt.plot(x, dataN, color = 'red', label = 'NS')
    plt.xlim(-2, 65)
    
    plt.legend(loc=1, borderaxespad=0.)
    plt.gca().yaxis.get_major_formatter().set_powerlimits((0, 1))
    plt.ylabel('Velocity (m/s)')
    plt.ylim(min(dataN), max(dataN))
    plt.annotate('seconds', va='center',  xy = (0.48,-0.12), xycoords = 'axes fraction')
    plt.annotate('(a)', xy = (0.01,0.90), xycoords = 'axes fraction', fontsize = 24, weight = 'bold')
    
    plt.subplot(gs.new_subplotspec((0, 1), colspan=1))
#    plt.title('event ' +  eventid + ', magnitude ' + str(ml) + ', recorded on station ' + station + ', distance ' + str(round(dist_km,2)) + ' km ')
    x = np.arange(0, len(dataE), 1)/100.
#    plt.plot(x, dataE, color='black', label = 'EW')
    plt.plot(x, dataN, color = 'black', label = 'NS')
    plt.xlim(-2, 65)
    
    plt.legend(loc=1, borderaxespad=0.)
    plt.gca().yaxis.get_major_formatter().set_powerlimits((0, 1))
    plt.ylim(min(dataN), max(dataN))
#    plt.gca().yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2e'))
#    plt.ylabel('Velocity (m/s)')
#    plt.xlabel('counts')
    plt.annotate('seconds', va='center',  xy = (0.48,-0.12), xycoords = 'axes fraction')
    plt.annotate('(b)', xy = (0.01,0.9), xycoords = 'axes fraction', fontsize = 24, weight = 'bold')

    plt.subplot(gs.new_subplotspec((1, 0), colspan=1))
    plt.xlim(0.5,70)
    plt.grid()
    plt.loglog(f_bins, record_spec, label = 'record', c= 'cornflowerblue', lw = 2)
    plt.tick_params(axis='both', which='major')
    plt.tick_params(axis='both', which='both', length = 5, width = 1)
    plt.legend(loc=3, borderaxespad=0.)
    plt.annotate('(c)',  xy = (0.01,0.9), xycoords = 'axes fraction', fontsize = 24, weight = 'bold')
    fig.text(0.01, 0.4, 'Velocity amplitude (m)', va='center', rotation='vertical')


    plt.subplot(gs.new_subplotspec((1, 1), colspan=1))
    plt.xlim(0.5,70)
    plt.ylim(0.0001, 0.3)
    plt.grid()
    plt.loglog(f_bins, record_spec*dist_m, label = 'record dist corr', c= 'cornflowerblue', lw =2)
    plt.loglog(f_bins, station_spec*event_spec, color='black', label = 'event*site', lw = 2)
    plt.legend(loc=3, borderaxespad=0.)
    plt.tick_params(axis='both', which='major')
    plt.tick_params(axis='both', which='both', length = 5, width = 1)
    plt.annotate('(d)',  xy = (0.01,0.9), xycoords = 'axes fraction', fontsize = 24, weight = 'bold')


    plt.subplot(gs.new_subplotspec((2, 0), colspan=1))
    plt.xlim(0.5,70)
    plt.ylim(0.1, 5)
    plt.xlabel('Frequency (Hz)')
    plt.grid()
    plt.loglog(f_bins, event_spec, color = 'green', label = 'event', lw = 2)
#    shift = np.mean((event_spec[27:70]))-np.mean((Brune[27:70]))
    plt.loglog(f_bins, Brune, ls = '--', color = 'green', label = 'Brune source shift', lw = 2)
    plt.legend(loc=3, borderaxespad=0.)
    plt.tick_params(axis='both', which='major')
    plt.tick_params(axis='both', which='both', length = 5, width = 1)
    plt.annotate('(e)',  xy = (0.01,0.9), xycoords = 'axes fraction', fontsize = 24, weight = 'bold')


    
    plt.subplot(gs.new_subplotspec((2, 1), colspan=1))
    plt.xlim(0.5,70)
    plt.ylim(0.0001, 0.1)
    plt.grid()
    plt.loglog(f_bins, station_spec, color='r', label = 'station ' + station, lw = 2)
    plt.legend(loc=3, borderaxespad=0.)
    plt.tick_params(axis='both', which='major')
    plt.tick_params(axis='both', which='both', length = 5, width = 1)
    plt.xlabel('Frequency (Hz)')
    plt.annotate('(f)', xy = (0.01,0.9), xycoords = 'axes fraction', fontsize = 24, weight = 'bold')
    plt.savefig(boxpath + '/secondo_rebin3/plot_' + eventid + '.' + network + '.' +  station + '.png')
    plt.savefig('/Volumes/USGS_Data/notes/kappa_paper/fig2_2013_noshift.png')
    plt.show()
#    plt.close()
#plt.close('all')