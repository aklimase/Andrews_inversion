#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 12 11:36:07 2017

@author: escuser
uses the catalog from a certain box and makes a file of all events in the box and stations they were recorded on
"""
import numpy
from numpy import genfromtxt
import os

##################################################
box = '6'

path='/Users/escuser/project/catalog'
events_file='/Users/escuser/project/catalog/all_' + box + '.txt'
pick_files='/Users/escuser/project/catalog/all_' + box

#get info from event text file, all events in box
events=genfromtxt(events_file,usecols=5,dtype='S')
origin_time=genfromtxt(events_file,usecols=4,dtype='S')

#get list of all events with .pick files from the all directory
pick_list = []
for root, dirs, files in os.walk(pick_files):
    for file in files:
        if file.endswith('.pick'):
            pick_list.append(os.path.splitext(file)[0])

#if event is also in the pick list, keep it
#intersection of events in box and events with .picks for our stations
pick_event_list = []
pick_origin_time = []
intersection = set(pick_list).intersection(events)
for i in range(len(events)):
    if events[i] in intersection:
        pick_event_list.append(events[i])
        pick_origin_time.append(origin_time[i])
print(len(pick_event_list))

####################################################
#make a new event list file
events_file=path+'/'+'event_list_' + box + '.txt'
f=open(events_file,'w')
f.write('#event_id, origin time, station\n')
            

#for each event in the given box
#print the event id, origin time, and stations from pick files
for i in range(len(pick_event_list)):
    f.write(pick_event_list[i] +'\t' + pick_origin_time[i]+'\t')

    #open pick file
    pick = open(pick_files +'/'+ pick_event_list[i] + '.pick')
    stations=genfromtxt(pick, usecols=2, dtype='S')#list of all stations
    stations = numpy.atleast_1d(stations)#because throws error when only one element in list with genfromtxt
    
    if(len(stations)>1):
        station_list = list(set(stations))
    else:
        station_list = stations

    for j in range(len(station_list)):
        f.write(station_list[j] + '\t')
    f.write('\n')
    
    
f.close()