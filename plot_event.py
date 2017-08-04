#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 31 15:57:15 2017

@author: escuser
"""

import matplotlib.pyplot as plt
import glob
import os
import os.path as path
from obspy import read


#box = 'Riverside_FRD_RDM'
#box = 'Imperial_Valley_PFO_TPFO_PMD'
#box = 'Imperial_Valley_SWS_ERR'
#box = 'Salton_Trough_SWS_ERR'

box = 'all_paths_subset'


boxpath = '/Users/escuser/project/boxes/' + box


#uncorrectedfilepath = glob.glob(boxpath + '/uncorrected/Event_*/*.SAC')#full path for only specified channel

event_dirs = glob.glob(boxpath + '/uncorrected/Event_*')
#events = []
#for i in range(len(event_dirs)):
#    events.append(path.basename(event_dirs[i]))
#make event directories within corrected local data
#make a directory for each event
#for i in range(len(events)):
#    if not path.exists(boxpath + '/event_plots/' + events[i]):
#        os.makedirs(boxpath + '/event_plots/' + events[i])

#
#for each event, plot 
for i in range(len(event_dirs)):#in this case event paths are all sac files
    station_files = glob.glob(event_dirs[i] + '/*HHE*.SAC')
    if len(station_files) > 0:
        print(station_files)
        event = station_files[0].split('/')[-2]
        num_subplots = len(station_files)
        fig = plt.figure(figsize = (25,25))
        fig.text(0.04, 0.5, 'Velocity amplitude', va='center', rotation='vertical', fontsize = 20)
        title = event_dirs[i] + ' HHE'
        plt.suptitle(title, fontsize = 20)
    
        for j in range(len(station_files)):
            base = path.basename(station_files[j])
#           event = station_files[j].split('/')[-2]
            network = base.split('_')[0]
            station = base.split('_')[1]
            full_channel = base.split('_')[2]
         
            stream = read(station_files[j])
            tr = stream[0]
            data = tr.data
            mag = str(tr.stats.sac.mag)
            plt.subplot(num_subplots, 1, j+1)
            plt.plot(data, color='black')
            plt.title(network + ' ' +  station + ' ' + mag, fontsize = 20)

        print 'saving image: ' + boxpath + '/event_plots/' + event + '.png'
        plt.savefig(boxpath + '/event_plots/' + event + '.png')
        plt.close()