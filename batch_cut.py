#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 27 16:20:24 2017

@author: escuser
"""

from library import cut_swave
import glob
import os.path as path
from obspy import read
import matplotlib.pyplot as plt


#for every event and station, cut the icorr file and plot the before and after
location = '*'
tsunit = 'VEL'
channel = 'HH*'

boxpath = '/Users/escuser/Documents/Alexis_Data/cut_sac_files'
#eventpaths = glob.glob(boxpath + '/icorrdata/*.AZ.*.' + channel + '.sac')#full path for only specified channel
eventpaths = glob.glob(boxpath + '/rawdata/*.*.*.' + channel + '.sac')#full path for only specified channel
print(len(eventpaths))

cutfilepath = boxpath + '/cutdata_s'


for i in range(len(eventpaths)):#in this case event paths are all sac files
#for i in range(0, 100, 2):#in this case event paths are all sac files
#loop through in groups of 3
    base = path.basename(eventpaths[i])
    print(eventpaths[i])
    eventid = base.split('.')[0]
    network = base.split('.')[1]
    station = base.split('.')[2]
    full_channel = base.split('.')[3]
    cutfile = cutfilepath + '/' + eventid + '.' + network + '.' + station + '.' + full_channel + '.sac'
    
    #calling cut function
    cuttime = cut_swave(eventpaths[i], cutfile)

    #read in traces from eventpaths[i] and cutfile
    stream = read(eventpaths[i])
    tr = stream[0]
    trace_starttime = tr.stats.starttime
    cut_xval = (cuttime - trace_starttime)/0.01
    data = tr.data
#    fig = plt.figure(figsize = (25,20))
#    plt.subplot(4,1,1)
#    plt.plot(data, color='black')
#    plt.title(base)
#    plt.xlim(0, len(data))
#    plt.axvline(x=cut_xval, c = 'r')#first cut
#    plt.axvline(x=cut_xval + 120/0.01, c = 'r')#second cut
    
    stream_cut = read(cutfile)
    tr_cut = stream_cut[0]
    data_cut = tr_cut.data
#    plt.subplot(4,1,2)
#    plt.plot(data_cut, color='black')
#    plt.xlim(0, len(data_cut))
    
    
#    
#    base = path.basename(eventpaths[i+1])
#    eventid = base.split('.')[0]
#    network = base.split('.')[1]
#    station = base.split('.')[2]
#    full_channel = base.split('.')[3]
#    cutfile = cutfilepath + '/' + eventid + '.' + network + '.' + station + '.' + full_channel + '.sac'
#    
#    #calling cut function
#    cuttime = cut_swave(eventpaths[i+1], cutfile)
#
#    stream = read(eventpaths[i+1])
#    tr = stream[0]
#    trace_starttime = tr.stats.starttime
#    cut_xval = (cuttime - trace_starttime)/0.01
#    data = tr.data
##    plt.subplot(4,1,3)
##    plt.plot(data, color='black')
##    plt.title(base)
##    plt.xlim(0, len(data))
##    plt.axvline(x=cut_xval, c = 'r')#first cut
##    plt.axvline(x=cut_xval + 120/0.01, c = 'r')#second cut
#    
#    stream_cut = read(cutfile)
#    tr_cut = stream_cut[0]
#    data_cut = tr_cut.data
##    plt.subplot(4,1,4)
##    plt.plot(data_cut, color='black')
##    plt.xlim(0, len(data_cut))
##    
##    plt.show()



