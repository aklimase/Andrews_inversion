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
#from obspy import core.utcdatetime

#box = 'Riverside_FRD_RDM'
box = 'Salton_Trough_SWS_ERR'
#box = 'Imperial_Valley_PFO_TPFO_PMD'
#box = 'Imperial_Valley_SWS_ERR'


boxpath = '/net/anzanas.wr.usgs.gov/data/users/alexis/ANZA_boxes/' + box
event_dirs = glob.glob(boxpath + '/uncorrected/Event_*')

localdir = '/Users/escuser/project'

#read in catalog file
catalog = localdir + '/Imperial_Valley_M2.5_USGS_Catalog.txt'

time_cat = np.genfromtxt(catalog, comments = '#', delimiter = '|', dtype = None, usecols = [1])
#need to truncate milliseconds
f = np.genfromtxt(catalog, comments = '#', delimiter = '|', dtype = None, usecols = [2,3,4,10])
lat, lon, depth, mag = f.T
time_cat = [obspy.core.utcdatetime.UTCDateTime(x.split('.')[0]) for x in time_cat]

events = []
for i in range(len(event_dirs)):
    events.append(path.basename(event_dirs[i]))
#make event directories within corrected local data
#make a directory for each event
for i in range(len(events)):
    if not path.exists(localdir + '/boxes/' + box + '/uncorrected/' + events[i]):
        os.makedirs(localdir + '/boxes/' + box + '/uncorrected/' + events[i])



#for each event directory, find match in catalog file
for i in range(len(events)):
    base = path.basename(events[i])
    (yyyy, month, day, hh, mm, ss) = [int(s) for s in base.split('_')[1:]]
    time = obspy.core.utcdatetime.UTCDateTime(yyyy, month, day, hh, mm, ss)
#    print(time)
#    print(cat_dict[str(time)])
#    #check and see if ev
    if time in time_cat:
        ind = time_cat.index(time)
    #go into each sac file in directory and write location to header
    files = glob.glob(event_dirs[i] + '/*.SAC')
    print 'writing: ' + events[i]
    for j in range(len(files)):
        sac_base = path.basename(files[j])
        stream = read(files[j])
        tr = stream[0]
        tr.stats.sac.evlo = lon[ind]
        tr.stats.sac.evla = lat[ind]
        tr.stats.sac.evdp = depth[ind]
        tr.stats.sac.mag = mag[ind]
        #origin time
        tr.stats.sac.nzhour = hh
        tr.stats.sac.nzmin = mm
        tr.stats.sac.nzsec = ss
        tr.stats.sac.nzmsec = 0
        stream.write(localdir + '/boxes/' + box + '/uncorrected/' + events[i] + '/' + sac_base, format='SAC')
