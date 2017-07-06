#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 28 14:46:27 2017

@author: escuser
"""

#plot raw data with lines, cut data, and corrected data
import matplotlib.pyplot as plt
import glob
import os.path as path
from obspy import read

channel = 'HH*'

boxpath = '/Users/escuser/Documents/Alexis_Data/cut_sac_files'
#eventpaths = glob.glob(boxpath + '/icorrdata/*.AZ.*.' + channel + '.sac')#full path for only specified channel
rawpaths = glob.glob(boxpath + '/rawdata/*.AZ.*.' + channel + '.sac')#full path for only specified channel
cutfilepaths = glob.glob(boxpath + '/cutdata_s/*.AZ.*.' + channel + '.sac')#full path for only specified channel
icorrfilepaths = glob.glob(boxpath + '/icorrdata/*.AZ.*.' + channel + '.sac')#full path for only specified channel

cutfilepath = boxpath + '/cutdata_s'#full path for only specified channel
icorrfilepath = boxpath + '/icorrdata'#full path for only specified channel


for i in range(20,40):#in this case event paths are all sac files
#for i in range(0, 100, 2):#in this case event paths are all sac files
#loop through in groups of 3
    base = path.basename(rawpaths[i])
    eventid = base.split('.')[0]
    network = base.split('.')[1]
    station = base.split('.')[2]
    full_channel = base.split('.')[3]
    
    rawfile = rawpaths[i]
    cutfile = cutfilepath + '/' + eventid + '.' + network + '.' + station + '.' + full_channel + '.sac'
    icorrfile = icorrfilepath + '/' + eventid + '.' + network + '.' + station + '.' + full_channel + '.sac'

    fig = plt.figure(figsize = (25,20))
    fig.suptitle(base, fontsize = 20)
    #read in traces from eventpaths[i] and cutfile
    stream = read(rawfile)
    tr = stream[0]
    data = tr.data
    
    plt.subplot(3,1,1)
    plt.plot(data, color='black')
    plt.title('raw')
    plt.xlim(0, len(data))
    
    stream = read(cutfile)
    tr = stream[0]
    data = tr.data
    
    plt.subplot(3,1,2)
    plt.plot(data, color='black')
    plt.title('cut')
    plt.xlim(0, len(data))
    
    stream = read(icorrfile)
    tr = stream[0]
    data = tr.data
    
    plt.subplot(3,1,3)
    plt.plot(data, color='black')
    plt.title('instrument corrected')
    plt.xlim(0, len(data))
    

    plt.show()
