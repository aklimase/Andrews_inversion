#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue May  8 15:09:52 2018

@author: temp

read in source parameters file (evid and stress drop)
read in residuals pickle object
match events based on event location
plot event residual vs. stress drop
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
#import obspy
#import glob
import dread
from scipy.stats import pearsonr
import statsmodels.stats.power
#import matplotlib.patches as mpatches
import cPickle as pickle
from temp_obj import Tobj

####set the maxiumum distance for finding matches here, km
match_dist = 20

plt.rcParams["font.family"] = "serif"
mpl.rcParams.update({'font.size': 22})
topdir = '/Volumes/USGS_Data/project/'

# residuals
rpath = topdir + 'event_boxes/mixedcoeff_v3anza2013_pga_noVs30_5coeff_a4_-1.2_pga__ncoeff5_Mc_8.5_VR_99.6_a4_-1.2_robj.pckl'
rfile = open(rpath,'r')
robj = pickle.load(rfile)
rfile.close()

r_evlat = robj.elat
r_evlon = robj.elon
r_evdep = robj.edepth
r_Eresidual = robj.E_terms

# stress drops
source_file =  topdir + 'source_params/secondo_2013_source.out'
data = np.genfromtxt(source_file ,delimiter='\t', dtype = None, names = True, encoding = None)

evid = data['evid']
log_s = data['log_stress_drop]


