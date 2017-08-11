#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  8 15:14:47 2017

@author: escuser
"""

import numpy as np
import matplotlib.pyplot as plt
import glob
import os.path as path

constraint_file =  '/Volumes/USGS_Data/project/boxes/all_paths/constraint_2011_08_18_20_26_14.out'
con = np.genfromtxt(constraint_file)
cf_spec = np.exp(con.T[1]) #second col
print(cf_spec)
secondo_ev =  glob.glob('/Volumes/USGS_Data/project/boxes/all_paths/secondo/2*.out')
secondo_stn = glob.glob('/Volumes/USGS_Data/project/boxes/all_paths/secondo/[!2]*.out')
freq_list = con.T[0] 
outfile_path = '/Volumes/USGS_Data/project/boxes/all_paths/secondo_constrained_Brune'

for i in range(len(secondo_ev)):#for each event
    #make each row into an array
    #constrained with station so mult by station spec
    event = np.genfromtxt(secondo_ev[i])
    eventid = path.basename(secondo_ev[i]).split('.')[0]
    amp = event.T[1]/cf_spec
    outfile = open(outfile_path + '/' + eventid + '.out', 'w')
    out = (np.array([freq_list, amp])).T
    outfile.write('#freq_bins \t vel_spec_NE_cm \n')
    np.savetxt(outfile, out, fmt=['%E', '%E'], delimiter='\t')
    outfile.close()
    
#    plt.figure(figsize = (14,8))
#    plt.loglog(freq_list, amp, color='cornflowerblue')
#    plt.title('event: ' + eventid)
#    plt.xlabel('frequency (Hz)')
#    plt.ylabel('velocity spectrum cm/s')
#    plt.grid()
#    plt.savefig(outfile_path + '/' + eventid + '.png')
#    plt.close()
    
for i in range(len(secondo_stn)):#for each station
    #make each row into an array
    stn = np.genfromtxt(secondo_stn[i])
    stnid = path.basename(secondo_stn[i]).split('.')[0]
    amp = stn.T[1]*cf_spec
    outfile = open(outfile_path + '/' + stnid + '.out', 'w')
    out = (np.array([freq_list, amp])).T
    outfile.write('#freq_bins \t vel_spec_NE_cm \n')
    np.savetxt(outfile, out, fmt=['%E', '%E'], delimiter='\t')
    outfile.close()
    
    
#    plt.figure(figsize = (14,8))
#    plt.loglog(freq_list, amp, color='cornflowerblue')
#    plt.title('station: ' + stnid)
##    plt.ylim(0.1, 10)
#    plt.xlabel('frequency (Hz)')
#    plt.ylabel('velocity spectrum cm/s')
#    plt.grid()
#    plt.savefig(outfile_path + '/' + stnid + '.png')
#    plt.close()