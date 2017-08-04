#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  2 09:26:23 2017

@author: escuser

before running delete everything in all_paths_subset/rawdata
"""

import random
import numpy as np
import obspy
import glob
import os.path as path
import shutil


lines = file('/net/anzanas.wr.usgs.gov/data/users/alexis/all_paths_M2.5_USGS_Catalog.txt').read().splitlines()
lines = lines[1:] #get rid of header
print(lines[0])
n = len(lines)
print(n)
print(0.1*n)
n_sample = int(np.round(0.1*n))

#list of line indexes
sample = random.sample(xrange(0, n-1), n_sample)

#get event info from those lines
print(sample)

sub_catalog = [lines[i] for i in sample]
print(len(sub_catalog))

#
time_cat = np.genfromtxt(sub_catalog, comments = '#', delimiter = '|', dtype = None, usecols = [1])

##need to truncate milliseconds

time_cat = [obspy.core.utcdatetime.UTCDateTime(x.split('.')[0]) for x in time_cat]


for i in range(0, len(time_cat)):
    yyyy = str(time_cat[i].year).zfill(4)
    month = str(time_cat[i].month).zfill(2)
    day = str(time_cat[i].day).zfill(2)
    hh = str(time_cat[i].hour).zfill(2)
    mm = str(time_cat[i].minute).zfill(2)
    ss = str(time_cat[i].second).zfill(2)
    base = 'Event_' + yyyy + '_' + month + '_' + day + '_' + hh + '_' + mm + '_' + ss
    src = '/net/anzanas.wr.usgs.gov/data/users/alexis/ANZA_boxes/all_paths/rawdata/' + base
    dst =  '/net/anzanas.wr.usgs.gov/data/users/alexis/ANZA_boxes/all_paths_subset/rawdata/' + base
    shutil.copytree(src, dst)

    
