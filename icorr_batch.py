#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 19 13:35:41 2017

@author: escuser
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
import matplotlib.pyplot as plt

############################################################################
#makes response with wf.download_response
#demeans, detrends and prefilters data
#removes instrument response

# Set parameters for resp file:
location = '*'

# Unit for time series?
#tsunit = 'VEL'
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
# for magnitude 4-5 the corner frequency is below 1 Hz so play around with LB
prefilt = (0.0,0.001,0.7*nyquistf,nyquistf)

### Could batch this - so could do it in two steps:
# 1.   Read in all files you have, and for every station, download a resp file
#       for all HH channels at once (so set location = '*', channel = 'HH*').
#      Set it to save all these to a resp file directory.
#
# 2.   Read in all uncorrected sac files in a directory, loop through and for each
#       one, look in the resp file directory and find the one that corresponds to 
#       your station/instrument in the loop (can use "split" to get this info, I 
#       can help if you want).  Then, use wf.remove_response inside the loop for 
#       each file, save to a directory of corrected sac files. 


#boxpath = '/net/anzanas.wr.usgs.gov/data/users/alexis/ANZA_boxes/Riverside_FRD_RDM'
boxpath = '/Users/escuser/Documents/Alexis_Data/cut_sac_files'
#eventpaths = glob.glob(boxpath + '/rawdata/*.AZ.*.' + channel + '.sac')#full path for only specified channel
eventpaths = glob.glob(boxpath + '/cutdata_s/*.TA.*.' + channel + '.sac')#full path for only specified channel

events = []
for i in range(len(eventpaths)):
    #events.append(eventpaths[i].split('/')[9])#take last part of path aka each folder name
    eventid = (eventpaths[i].split('/')[7]).split('.')[0]
    if eventid not in events:
        events.append(eventid)
print(len(events))

##read in all files to find networks and stations
networklist = []
stationlist = []

     
for i in range(len(eventpaths)):
    #get the path of every file in the event directory and split into network.staation, and channel
    base = path.basename(eventpaths[i])
    network = base.split('.')[1]
    station = base.split('.')[2]
    #if network and station not part of the list yet add
    if network not in networklist:
        networklist.append(network)
    if station not in stationlist:
        stationlist.append(station)
print(networklist)
print(len(stationlist))

#make response files in response directory for each combo of networks and stations
#response_path = '/net/anzanas.wr.usgs.gov/data/users/alexis/ANZA_boxes/response_files'
response_path = '/Users/escuser/Documents/Alexis_Data/cut_sac_files/response_files'

for network in networklist:
    for station in stationlist:
#        respfile = response_path + '/' + network + '_' + station + '_' + 'HH*.resp'
#        channel = 'HH*'
#        wf.download_response(network,station,location,channel,starttime,endtime,respfile)
#        respfile = response_path + '/' + network + '_' + station + '_' + 'BH*.resp'
#        channel = 'BH*'
#        wf.download_response(network,station,location,channel,starttime,endtime,respfile)
        respfile = response_path + '/' + network + '.' + station + '.' + channel + '.resp'
        wf.download_response(network,station,location,channel,starttime,endtime,respfile)

        
##read in all uncorrected sac files and loop through and for each look in response directory, then remove response and save to directory of corrected files
#icorr_path = boxpath + '/icorrdata'
##make a directory for each event
#for i in range(len(events)):
#    if not os.path.exists(icorr_path + '/' + events[i]):
#        os.makedirs(icorr_path + '/' + events[i])
##read in raw data, rmean and rtrend, find .resp file, correct and add to icorr dir
##for i in range(len(eventpaths)):#in this case event paths are all sac files
#for i in range(len(eventpaths)):#in this case event paths are all sac files
#    base = path.basename(eventpaths[i])
#    eventid = base.split('.')[0]
#    network = base.split('.')[1]
#    station = base.split('.')[2]
#    full_channel = base.split('.')[3]
#    #find response file
#    respfile = response_path + '/' + network + '.' + station + '.' + channel + '.resp'
#    #first rmean and rtrend
#    stream = read(eventpaths[i])
#    tr = stream[0]
#    tr.detrend(type = 'demean')#removes mean of data
#    tr.detrend(type = 'simple')#rtrend linear from first and last samples
#    #rewrite to a sac file
#    tr.write('temp.sac', format = 'sac')
#    sacfile = 'temp.sac'
#    icorr_sacfile = icorr_path + '/' + eventid + '.' + network + '.' + station + '.' + full_channel + '.sac'
#    print(icorr_sacfile)
#    #uncorrected_sac_file,resp_file,corrected_sac_file,water_level_bds,resp_unit
##    icorr_sacfile = sacfile # no corrections yet
#    wf.remove_response(sacfile,respfile,icorr_sacfile,prefilt,tsunit)#prefilt values for HH
    
    