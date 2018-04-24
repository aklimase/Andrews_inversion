#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 28 14:46:27 2017

@author: escuser

inputs: path to raw data, cut data, and corrected data

reads in all three sac files for a record and makes a plot showing all three 

outputs: generates plot to the screen, can change to save images as .png s
""" 
import matplotlib.pyplot as plt
import glob
import os
import os.path as path
from obspy import read

top_dir = '/Volumes/USGS_Data/project'

#box = 'Imperial_Valley_PFO_TPFO_PMD'
#box = 'Imperial_Valley_SWS_ERR'
#box = 'Riverside_FRD_RDM'
#box = 'Salton_Trough_SWS_ERR'
box = 'all_paths'

boxpath = top_dir + '/boxes/' + box


#event_dirs = glob.glob(boxpath + '/corrected/Event_*')

#cutfilepaths = glob.glob(boxpath + '/cutdata_s/Event_*/*.SAC')

uncorr = boxpath + '/uncorrected/'
cut = boxpath + '/cutdata_s/'
#corr = boxpath + '/corrected/'

#cutfilepath = boxpath + '/cutdata_s'#full path for only specified channel
correctedfilepath = glob.glob(boxpath + '/corrected/Event_2011_04*/*.SAC')#full path for only specified channel

event_dirs = glob.glob(boxpath + '/corrected/Event_*')
events = []
for i in range(len(event_dirs)):
    events.append(path.basename(event_dirs[i]))
#make event directories within corrected local data
#make a directory for each event
for i in range(len(events)):
    if not path.exists(boxpath + '/plots/' + events[i]):
        os.makedirs(boxpath + '/plots/' + events[i])



for i in range(len(correctedfilepath)):#in this case event paths are all sac files
    base = path.basename(correctedfilepath[i])
    event = correctedfilepath[i].split('/')[-2]
    network = base.split('_')[0]
    station = base.split('_')[1]
    full_channel = base.split('_')[2]
        
    uncorrfile = uncorr + event + '/' + base
    corrfile = correctedfilepath[i]
    cutfile = cut + event + '/' + base
    
    #read in traces
    stream = read(uncorrfile)
    tr = stream[0]
    data1 = tr.data
#    start = date2num(tr.stats.starttime.datetime)
#    end = date2num(tr.stats.endtime.datetime)
#    times = np.linspace(start, end, tr.stats.npts)
    
    stream = read(corrfile)
    tr = stream[0]
    data2 = tr.data
    mag = str(tr.stats.sac.mag)    

    fig = plt.figure(figsize = (25,20))
    fig.text(0.04, 0.5, 'Velocity amplitude', va='center', rotation='vertical', fontsize = 15)
    title = base + ' magnitude ' + mag
    fig.suptitle(title, fontsize = 20)

    plt.subplot(3,1,1)
    plt.plot(data1, color='black')
    plt.title('uncorrected', fontsize = 16)
    plt.ylabel('counts', fontsize = 16)
#    plt.xlim(times[0], times[-1])
    
#    stream = read(corrfile)
#    tr = stream[0]
#    data2 = tr.data
#    mag = tr.stats.sac.mag
#    print(mag)
    
    plt.subplot(3,1,3)
    plt.plot(data2, color='black')
    plt.title('corrected', fontsize = 16)
    plt.ylabel('m/s', fontsize = 16)
    plt.xlim(0, len(data2))
    plt.xlabel('samples', fontsize = 16)
    
    stream = read(cutfile)
    tr = stream[0]
    data3 = tr.data
    
    plt.subplot(3,1,2)
    plt.plot(data3, color='black')
    plt.title('cut', fontsize = 16)
    plt.ylabel('counts', fontsize = 16)
    plt.xlim(0, len(data3))

    
#    resp = data3/data2
#    print(min(resp), max(resp))
#    plt.subplot(4,1,4)
#    plt.plot(resp, color='black')
#    plt.title('response')
#    plt.xlim(0, len(data3))
#    plt.xlabel('sample number', fontsize = 15)
    
    print 'saving image: ' + boxpath + '/plots/' + event + '/' + base.split('.')[0] + '.png'
    plt.savefig(boxpath + '/plots/' + event + '/' + base.split('.')[0] + '.png')
    plt.close()
    plt.show()
    
    
    
    
    