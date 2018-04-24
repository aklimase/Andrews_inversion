#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 15 11:13:48 2017
 for each record, take record - event*station and write a file
@author: alexisklimasewski
"""


import numpy as np
import os.path as path
import glob
import dread
from obspy import read
import obspy
from pyproj import Geod

g = Geod('+ellps=clrk66')

top_dir = '/Volumes/USGS_Data/project'

boxpath = top_dir + '/boxes/'
box_list = [('Imperial_Valley', 'PFO_TPFO_PMD'), ('Imperial_Valley', 'SWS_ERR'), ('Riverside', 'FRD_RDM'), ('Salton_Trough', 'SWS_ERR'), ('all_paths', '')]
residual_list = []

evpaths = glob.glob(boxpath + 'all_paths/record_spectra/Event_*')
ev = [x.split('/')[-1] for x in evpaths]

print(len(evpaths))

#read in corrected event directories to find event matches
for i in range(4,5):
    print box
    box = box_list[i][0]
    stations = box_list[i][1]
    st_list = stations.split('_')
    #read in catalog file
    catalog = top_dir + '/' + box + '_M2.5_USGS_Catalog.txt'

    time_cat = np.genfromtxt(catalog, comments = '#', delimiter = '|', dtype = None, usecols = [1])
    #need to truncate milliseconds
    f = np.genfromtxt(catalog, comments = '#', delimiter = '|', dtype = None, usecols = [2,3,4,10])
#    lat, lon, depth, mag = f.T
    #times of all events in the box catalog
    time_cat = [obspy.core.utcdatetime.UTCDateTime(x.split('.')[0]) for x in time_cat]
    
    #check all events in record_spectra to see if time_cat event is there
    #if yes add to record_paths
    box_events = []
    #make a record path list from box events in time_cat
    for k in range(len(time_cat)):
        yyyy = str(time_cat[k].year).zfill(4)
        month = str(time_cat[k].month).zfill(2)
        day = str(time_cat[k].day).zfill(2)
        hh = str(time_cat[k].hour).zfill(2)
        mm = str(time_cat[k].minute).zfill(2)
        ss = str(time_cat[k].second).zfill(2)
        box_events.append('Event_' + yyyy + '_' + month + '_' + day + '_' + hh + '_' + mm + '_' + ss)
    
    print(len(box_events))
    record_paths = []
    for p in range(len(box_events)):
        if box_events[p] in ev:
            for st in st_list:
                f = '*' + st + '*.out'
#                print glob.glob((boxpath + 'all_paths_subset/record_spectra/' + box_events[p] + '/' + f))
                record_paths.extend(glob.glob(boxpath + 'all_paths/record_spectra/' + box_events[p] + '/' + f))
    print(box, len(record_paths))
    
    out_dir = boxpath + 'all_paths/secondo/'
    
    residual = np.zeros((len(record_paths), 50))

    mag_list = []
    dist_list = []
    az_list = []
    dep_list = []
    
    #for each event/station, get the information
    for j in range(len(record_paths)):
        base = path.basename(record_paths[j])
        network, station, channel, loc = base.split('_')[0:4]
        yyyy, month, day, hh, mm, ss = base.split('_')[4:]
        ss = ss.split('.')[0]
        eventid = yyyy + '_' + month + '_' + day + '_' + hh + '_' + mm + '_' + ss    
        
        #read in raw data for record info
        #correct for distance
        raw_file = boxpath + 'all_paths/uncorrected/Event_' + eventid + '/' + network + '_' + station + '_HHN_' + loc + '_' + eventid + '.SAC'
        stream = read(raw_file)
        tr = stream[0]
        
        evlon =  tr.stats.sac.evlo #deg
        evlat =  tr.stats.sac.evla #deg
        evdepth = tr.stats.sac.evdp #km
        stlon = tr.stats.sac.stlo #deg
        stlat = tr.stats.sac.stla #deg
        stdepth = tr.stats.sac.stdp #km
        mag = tr.stats.sac.mag
        az12,az21,dist = g.inv(evlon,evlat,stlon,stlat)
        mag_list.append(float(mag))
        dep_list.append(float(evdepth))
        az_list.append(az21)

        #find distance between event and station
        dist =  dread.compute_rrup(evlon, evlat, evdepth, stlon, stlat, stdepth) #in km
        dist_list.append(dist) # in km
        
        #km to cm
        dist = dist*100000
    
        #record spectra
        record_data = np.genfromtxt(record_paths[j], dtype = float, comments = '#', delimiter = None, usecols = (0,1))#only read in first two cols
        f_bins = record_data[:,0]
        
        record_spec = record_data[:,1]*dist
        #event spectra
#        print(out_dir + eventid + '.out')
        event_data = np.genfromtxt(out_dir + eventid + '.out', dtype = float, comments = '#', delimiter = None, usecols = (0,1))
        event_spec = event_data[:,1]
        #station spectra
#        print(out_dir + station + '.out')
        station_data = np.genfromtxt(out_dir + station + '.out', dtype = float, comments = '#', delimiter = None, usecols = (0,1))
        station_spec = station_data[:,1]

        residual[j:,] = np.log(record_spec) - np.log(station_spec*event_spec)
    
        outfile = open(top_dir + '/boxes/all_paths/residual_spectra/' + box + '/' + base, 'w')
        out = (np.array([f_bins, np.exp(np.log(record_spec) - np.log(station_spec*event_spec))])).T
        outfile.write('#freq_bins \t residual_spec_cm \n')
        np.savetxt(outfile, out, fmt=['%E', '%E'], delimiter='\t')
        outfile.close()
    