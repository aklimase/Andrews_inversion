#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 19 10:44:03 2017

@author: escuser

read in corrected sac files
writes event location to sac header using the catalog text file
write corrected sac files to local data directory
"""
import glob
import os
import os.path as path
import numpy as np
import obspy
from obspy import read


top_dir = '/Volumes/USGS_Data/project'
boxpath = top_dir + '/boxes/all_paths'

#list of paths to the SAC files without header
event_dirs = glob.glob(top_dir + '/boxes/all_paths/uncorrected_local/Event_*')

#read in catalog file with all of the header information
catalog = top_dir + '/catalogs/all_paths_M2.5_USGS_Catalog.txt'

#get the event location information for the header from the catalog file
data = np.genfromtxt(catalog, comments = '#', delimiter = '|', dtype = None, usecols = [2,3,4,10], names = True, encoding = None)
time_cat = data['Time']
#trim the trailing 0s
time_cat = [obspy.core.utcdatetime.UTCDateTime(x.split('.')[0]) for x in time_cat]
lat = data['Latitude']
lon = data['Longitude']
depth = data['Depthkm']
mag = data['Magnitude']

#list of all of my event directories
events = []
for i in range(len(event_dirs)):
    events.append(path.basename(event_dirs[i]))
    
#make event directories within corrected local data
#make a directory for each event so each record has a directory
#leave this out if you don't need files to go in event specific directories
for i in range(len(events)):
    if not path.exists(boxpath + '/uncorrected/' + events[i]):
        os.makedirs(boxpath + '/uncorrected/' + events[i])

#for each event directory, find match in catalog file
for i in range(len(events)):
    base = path.basename(events[i])
    #getting the time of my events from the event id string
    (yyyy, month, day, hh, mm, ss) = [int(s) for s in base.split('_')[1:]]
    #utc time object
    time = obspy.core.utcdatetime.UTCDateTime(yyyy, month, day, hh, mm, ss)
    #check and if the time is in the catalog, if so get the index
    if time in time_cat:
        ind = time_cat.index(time)
    ################################################################
    #go into each sac file in directory and write location to header
    #for that event get the SAC files for that event at each of the stations
    files = glob.glob(event_dirs[i] + '/*.SAC')
    print 'writing: ' + events[i]
    for j in range(len(files)):
        sac_base = path.basename(files[j])
        #use obspy to read in files
        stream = read(files[j])
        tr = stream[0]
        #event location
        tr.stats.sac.evlo = lon[ind]
        tr.stats.sac.evla = lat[ind]
        tr.stats.sac.evdp = depth[ind]
        tr.stats.sac.mag = mag[ind]
        #origin time
        tr.stats.sac.nzhour = hh
        tr.stats.sac.nzmin = mm
        tr.stats.sac.nzsec = ss
        tr.stats.sac.nzmsec = 0
        #station info included in SAC files
        stream.write(boxpath + '/uncorrected/' + events[i] + '/' + sac_base, format='SAC')
