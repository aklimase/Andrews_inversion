#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 19 13:35:41 2017

@author: escuser
based off of icorr.py

inputs: path to cut s sac files (cutdata\_s) and path to response files (response\_files)

used the instrument corrections from the iris client, run through obspy, filers from 35-50 Hz, no lower filer, keep long period noise, when uncommented, writes the .resp instrument response files for network and stations \\

outputs: writes instrument corrected sac files to icorrdata

"""

#### 
# Test waveforms.py to download resp file and use it to remove response
# VJS 6/2017
###

import waveforms as wf
import os.path as path
import os
from obspy import read
import glob

############################################################################
#makes response with wf.download_response
#demeans, detrends and prefilters data
#removes instrument response

# Set parameters for resp file:
location = '*'

# Unit for time series?
tsunit = 'VEL'

channel = 'HH*'

nyquistf = 50

# Start and end date/time for resp file:
starttime = '1998-01-01T00:00.000'
endtime = '2599-12-31T23:59.000'


# Pre-filter for instrument response correction:
#   The last number should be the nyquist frequency; second to last reduces
#   the filter towards the nyquist; first two are for the "water level"
# for BH Nyquist frequency is 20, for HH Nyquist frequency is 50
prefilt = (0.0,0.001,0.7*nyquistf,nyquistf)


#boxpath = '/net/anzanas.wr.usgs.gov/data/users/alexis/ANZA_boxes/Riverside_FRD_RDM'
#boxpath = '/Users/escuser/Documents/Alexis_Data/cut_sac_files'
#eventpaths = glob.glob(boxpath + '/cutdata_s/*.AZ.*.' + channel + '.sac')#full path for only specified channel



#box = 'Imperial_Valley_PFO_TPFO_PMD'
#box = 'Imperial_Valley_SWS_ERR'
#box = 'Riverside_FRD_RDM'
#box = 'Salton_Trough_SWS_ERR'
box = 'all_paths_subset'

boxpath = '/Users/escuser/project/boxes/' + box
event_dirs = glob.glob(boxpath + '/cutdata_s/Event_*')

eventpaths = glob.glob(boxpath + '/cutdata_s/Event_*/*.SAC')#full path
print 'Number of files: ', len(eventpaths)

cut_dir = boxpath + '/cutdata_s/'

events = [os.path.basename(x) for x in event_dirs]

#make a directory for each event
for i in range(len(events)):
    if not path.exists(cut_dir + '/' + events[i]):
        os.makedirs(cut_dir + '/' + events[i])


l = []   
for i in range(len(eventpaths)):
    #get the path of every file in the event directory and split into network.station, and channel
    base = path.basename(eventpaths[i])
    network = base.split('_')[0]
    station = base.split('_')[1]
    pair = network + '.' + station
    #if network and station not part of the list yet add
#    if network not in networklist:
#        networklist.append(network)
    if pair not in l:
        l.append(network + '.' + station)
#print(networklist)
print(l)

#make response files in response directory for each combo of networks and stations
#response_path = '/net/anzanas.wr.usgs.gov/data/users/alexis/ANZA_boxes/response_files'
response_path = boxpath + '/respfiles'

####################################################
#comment out if already have resp files
##makes a response file for each station and channel
##for network in networklist:
for i in range(len(l)):
    respfile = response_path + '/' + l[i] + '.' + channel + '.resp'
    network, station = l[i].split('.')
    wf.download_response(network,station,location,channel,starttime,endtime,respfile)

        
#read in all uncorrected sac files and loop through and for each look in response directory, then remove response and save to directory of corrected files
icorr_path = boxpath + '/corrected'
#make a directory for each event
for i in range(len(events)):
    if not path.exists(icorr_path + '/' + events[i]):
        os.makedirs(icorr_path + '/' + events[i])
        
#read in cut data, rmean and rtrend, find .resp file, correct and add to icorr dir
for i in range(len(eventpaths)):#in this case event paths are all sac files
    base = path.basename(eventpaths[i])
    print 'correcting file: ' + base
    folder = eventpaths[i].split('/')[-2]
    network = base.split('_')[0]
    station = base.split('_')[1]
    full_channel = base.split('_')[2]
    #find response file
    respfile = response_path + '/' + network + '.' + station + '.' + channel + '.resp'
    #first rmean and rtrend
    stream = read(eventpaths[i])
    tr = stream[0]
    #check and make sure that the trace isn't empty
    if(tr.stats.npts > 0):
        tr.detrend(type = 'demean')#removes mean of data
        tr.detrend(type = 'simple')#rtrend linear from first and last samples
        #rewrite to a sac file
        tr.write('temp.sac', format = 'sac')
        sacfile = 'temp.sac'
        icorr_sacfile = icorr_path + '/' + folder + '/'+ base
        #uncorrected_sac_file,resp_file,corrected_sac_file,water_level_bds,resp_unit
        wf.remove_response(sacfile,respfile,icorr_sacfile,prefilt,tsunit)#prefilt values for HH
    
    