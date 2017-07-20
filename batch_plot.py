#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 28 14:46:27 2017

@author: escuser

inputs: path to raw data, cut data, and corrected data

reads in all three sac files for a record and makes a plot showing all three 

outputs: generates plot to the screen, can change to save images as .png s
"""

#plot raw data, cut data, and corrected data
import matplotlib.pyplot as plt
import glob
import os.path as path
from obspy import read


#boxpath = '/Users/escuser/Documents/Alexis_Data/cut_sac_files'
box = 'Riverside_FRD_RDM'

boxpath = '/Users/escuser/project/boxes/' + box
event_dirs = glob.glob(boxpath + '/corrected/Event_*')

uncorrpaths = glob.glob(boxpath + '/uncorrected/Event_*/*.SAC')
corrfilepaths = glob.glob(boxpath + '/corrected/Event_*/*.SAC')
cutfilepaths = glob.glob(boxpath + '/cutdata_s/Event_*/*.SAC')

#cutfilepath = boxpath + '/cutdata_s'#full path for only specified channel
#icorrfilepath = boxpath + '/icorrdata'#full path for only specified channel


for i in range(0,10):#in this case event paths are all sac files
    base = path.basename(uncorrpaths[i])
    network = base.split('_')[0]
    station = base.split('_')[1]
    full_channel = base.split('_')[2]
    
    uncorrfile = uncorrpaths[i]
    corrfile = corrfilepaths[i]
    cutfile = cutfilepaths[i]

    fig = plt.figure(figsize = (25,20))
    fig.text(0.04, 0.5, 'Velocity amplitude (m/s)', va='center', rotation='vertical', fontsize = 15)
    fig.suptitle(base, fontsize = 20)
    #read in traces
    stream = read(uncorrfile)
    tr = stream[0]
    data = tr.data
    
    plt.subplot(3,1,1)
    plt.plot(data, color='black')
    plt.title('uncorrected')
    plt.xlim(0, len(data))
    
    stream = read(corrfile)
    tr = stream[0]
    data = tr.data
    
    plt.subplot(3,1,2)
    plt.plot(data, color='black')
    plt.title('corrected')
    plt.xlim(0, len(data))
    
    stream = read(cutfile)
    tr = stream[0]
    data = tr.data
    
    plt.subplot(3,1,3)
    plt.plot(data, color='black')
    plt.title('cut')
    plt.xlim(0, len(data))
    plt.xlabel('sample number', fontsize = 15)
    

    plt.savefig(boxpath + '/plots/' + base.split('.')[0] + '.png')
