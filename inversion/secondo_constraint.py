#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  8 15:14:47 2017

make a constraint file and put in the box directory
then constrain each event and station

if constraint is an event, divide event and multipy station
if constraint is station, multiply event and divide station

rename the outfile for the event or station

all in cm/s
"""

import numpy as np
import glob
import os.path as path


top_dir = '/Volumes/USGS_Data/project'
box = 'all_paths'
boxpath = top_dir + '/boxes/' + box

#######################################################
secondo_dir = 'secondo_rebin3'
constraint_file =  boxpath + '/constraint_rebin3_2010_05_25_19_49_51.out'
outfile_path = boxpath + '/secondo_rebin3_constrained_2010_05_25_19_49_51'
#########################################################

con = np.genfromtxt(constraint_file)
cf_spec = con.T[1] #second col
print(cf_spec)
secondo_ev =  glob.glob(boxpath + '/' + secondo_dir + '/2*.out')
secondo_stn = glob.glob(boxpath + '/'+ secondo_dir + '/[!2]*.out')
freq_list = con.T[0] 


##not in log space anymore
for i in range(len(secondo_ev)):#for each event
    #make each row into an array
    event = np.genfromtxt(secondo_ev[i])
    eventid = path.basename(secondo_ev[i]).split('.')[0]
    ######
    amp = event.T[1]/cf_spec
    outfile = open(outfile_path + '/' + eventid + '.out', 'w')
    out = (np.array([freq_list, amp])).T
    outfile.write('#freq_bins \t vel_spec_NE_m \n')
    np.savetxt(outfile, out, fmt=['%E', '%E'], delimiter='\t')
    outfile.close()
    
    
for i in range(len(secondo_stn)):#for each station
    #make each row into an array
    stn = np.genfromtxt(secondo_stn[i])
    stnid = path.basename(secondo_stn[i]).split('.')[0]
    ####
    amp = stn.T[1]*cf_spec
    outfile = open(outfile_path + '/' + stnid + '.out', 'w')
    out = (np.array([freq_list, amp])).T
    outfile.write('#freq_bins \t vel_spec_NE_m \n')
    np.savetxt(outfile, out, fmt=['%E', '%E'], delimiter='\t')
    outfile.close()