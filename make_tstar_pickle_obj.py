#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 15 12:57:33 2018

@author: temp

takes in the file with the t* values for all paths and assigns an event location for each event name
makes a t star path object and writes to a pickle file

"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import obspy
import glob
import dread
from scipy.stats import pearsonr
import statsmodels.stats.power
import matplotlib.patches as mpatches
import cPickle as pickle
from obj import Tobj


plt.rcParams["font.family"] = "serif"
mpl.rcParams.update({'font.size': 22})

topdir = '/Volumes/USGS_Data/project/'

## Define the path to open:
#rpath = topdir + 'event_boxes/mixedcoeff_v3anza2013_pga_noVs30_5coeff_a4_-1.2_pga__ncoeff5_Mc_8.5_VR_99.6_a4_-1.2_robj.pckl'
#
## Open the file to read ('r'), save to rfile:
#rfile = open(rpath,'r')
#
## Load the object inside of rfile into a variable (this could take a while, the object is big):
#robj = pickle.load(rfile)
#
## Close the file:
#rfile.close()

#define an object for the t* paths
#class Tobj(object):
#    evid = ""
#    tstar = 0
#    dist = 0
#    station = ""
#    azimuth = 0
#    evlat = 0
#    evlon = 0
#    evdep = 0
#
#    # The class "constructor" - It's actually an initializer 
#    def __init__(self, evid, evlat, evlon, evdep, tstar, dist, station, azimuth):
#        self.evid = evid
#        self.tstar = tstar
#        self.dist = dist
#        self.station = station
#        self.azimuth = azimuth
##
def make_tobj(evid, evlat, evlon, evdep, tstar, dist, station, azimuth):
    tobj = Tobj(evid, evlat, evlon, evdep, tstar, dist, station, azimuth)
    return tobj

##path t* catalog
infile = topdir + 'boxes/all_paths/t*plots/all_paths_tstar.out'
data = np.genfromtxt(infile, delimiter='\t', names = True, dtype = None, encoding = None)

slist = [data['network_station'][i].split('_')[1] for i in range(len(data['network_station']))]

tobj = make_tobj(data['eventid'], 0, 0, 0, data['tstar'], data['distance_km'], slist, data['azimuth'])


#station_set = set(tobj.station)

evlat_cat = []
evlon_cat = []
evdep_cat = [] 
#
catalog = topdir + '/all_paths_M2.5_USGS_Catalog.txt'
##
time_cat = np.genfromtxt(catalog, comments = '#', delimiter = '|', dtype = None, usecols = [1])
###need to truncate milliseconds
f = np.genfromtxt(catalog, comments = '#', delimiter = '|', dtype = None, usecols = [2,3,4,10])
lat, lon, depth, mag = f.T
time_cat = [obspy.core.utcdatetime.UTCDateTime(x.split('.')[0]) for x in time_cat]


for i in range(len(tobj.evid)):
    base = tobj.evid[i] 
    (yyyy, month, day, hh, mm, ss) = [int(s) for s in base.split('_')]
    time = obspy.core.utcdatetime.UTCDateTime(yyyy, month, day, hh, mm, ss)
    s = tobj.station[i]

    if time in time_cat:
        print(time)
        ind = time_cat.index(time)
        evlat_cat.append(lat[ind])
        evlon_cat.append(lon[ind])
        evdep_cat.append(depth[ind])

tobj.evlat = evlat_cat    
tobj.evlon = evlon_cat
tobj.evdep = evdep_cat


output = open(topdir + 'boxes/all_paths/path_tstar_evlocations.pckl', 'wb')
pickle.dump(tobj, output)
output.close()


tfile = open(topdir + 'boxes/all_paths/path_tstar_evlocations.pckl', 'r')
tobj2 = pickle.load(tfile)





