#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 29 11:36:09 2017

@author: escuser

inputs: reads in the text files of NE average spectra from record\_spectra

method: reads in record spectra and corrects for distances using dread function, 
uses the Andrews 86 method and svd to compute station and event spectra from record spectra, 
plots the event and station spectra, then if a constraint station or event file is specified, 
the output files are rewritten for the constraint and plots are generated again, 
if no constraint, comment out that part of the code

outputs: writes spectra text files for each event and station in event_site_spectra directory
"""
import glob
import os.path as path
import numpy as np
import matplotlib.pyplot as plt
from obspy import read    
import dread


#read in the cut and corrected spectra = records
#records are in cm/s
#should be binned frequencies and amplitudes
#make list of records and the corresponding events and stations

#box = 'Imperial_Valley_PFO_TPFO_PMD'
#box = 'Imperial_Valley_SWS_ERR'
#box = 'Riverside_FRD_RDM'
#box = 'Salton_Trough_SWS_ERR'

#boxpath = '/Users/escuser/project/boxes/' + box

#boxpath = '/Users/escuser/Documents/Alexis_Data/cut_sac_files'
#record_path = glob.glob(boxpath + '/record_spectra/*.AZ.*.out')
#print 'Number of records: ', len(record_path)

#record_path = glob.glob('/Users/escuser/project/boxes/Imperial_Valley_PFO_TPFO_PMD/record_spectra/Event_*/*.out')#full path for only specified channel
#record_path.extend(glob.glob('/Users/escuser/project/boxes/Imperial_Valley_SWS_ERR/record_spectra/Event_*/*.out'))
#record_path.extend(glob.glob('/Users/escuser/project/boxes/Riverside_FRD_RDM/record_spectra/Event_*/*.out'))
#record_path.extend(glob.glob('/Users/escuser/project/boxes/Salton_Trough_SWS_ERR/record_spectra/Event_*/*.out'))

record_path = glob.glob('/Users/alexisklimasewski/Documents/USGS/Riverside_FRD_RDM/record_spectra/Event_*/*.out')
print(record_path)

print 'Number of records: ', len(record_path)
#read in all files to find networks and stations
stationlist = []
stn_lat = []
stn_lon = []
eventidlist = []
event_lat = []
event_lon = []
event_depth = []
record_freq = []
record_spec = []

##############################################################################
## read in the uncut sac files and get the distances between source and station

for i in range(len(record_path)):
    record = (record_path[i].split('/')[-1])
    ##############################################################################3
    box =  (record_path[i].split('/')[5])
#    print(box)
    base = path.basename(record)
    network, station, channel, loc = base.split('_')[0:4]
    yyyy, month, day, hh, mm, ss = base.split('_')[4:]
    ss = ss.split('.')[0]
    eventid = yyyy + '_' + month + '_' + day + '_' + hh + '_' + mm + '_' + ss
    #read in uncorrected data for header info
    #correct for distance
#    raw_file = boxpath + '/uncorrected/Event_'+ eventid + '/' + network + '_' + station + '_HHN_' + loc + '_' + eventid + '.SAC'
    ################################
#    raw_file = '/Users/escuser/project/boxes/' + box + '/uncorrected/Event_'+ eventid + '/' + network + '_' + station + '_HHN_' + loc + '_' + eventid + '.SAC'
    raw_file = '/Users/alexisklimasewski/Documents/USGS/Riverside_FRD_RDM/uncorrected/Event_'+ eventid + '/' + network + '_' + station + '_HHN_' + loc + '_' + eventid + '.SAC'
    stream = read(raw_file)
    tr = stream[0]
    
    evlon =  tr.stats.sac.evlo #deg
    evlat =  tr.stats.sac.evla #deg
    evdepth = tr.stats.sac.evdp #km
    stlon = tr.stats.sac.stlo #deg
    stlat = tr.stats.sac.stla #deg
    stdepth = tr.stats.sac.stdp #km
    
    #find distance between event and station
    dist =  dread.compute_rrup(evlon, evlat, evdepth, stlon, stlat, stdepth) #in km
    #km to cm
    dist = dist*100000

    #read in file
    data = np.genfromtxt(record_path[i], dtype = float, comments = '#', delimiter = None, usecols = (0,1))#only read in first two cols
    record_freq.append(data[:,0])
    #data is NE spectra
    #square for power spectra
    record_spec.append((data[:,1]*dist)**2)
    #if network and station not part of the list yet add
    if station not in stationlist:
        stationlist.append(station)
        stn_lat.append(stlat)
        stn_lon.append(stlon)
#        print station, stlat, stlon
    if eventid not in eventidlist:
        eventidlist.append(eventid)
        event_lat.append(evlat)
        event_lon.append(evlon)
        event_depth.append(evdepth)
#        print eventid, evlat, evlon, evdepth
        

        
#F_bins = len(record_freq[0])#number of bins out data has
freq_list = record_freq[0] 
F_bins = len(freq_list)
print(F_bins)

#rows = len(record_path)
rows = len(record_path) #testing first 10
print(rows)

index_matrix = [[0 for j in range(3)] for i in range(rows)]

#for i in range(len(records)):
for i in range(rows):
    record = record_path[i].split('/')[-1]
    base = path.basename(record)
    network, station, channel, loc = base.split('_')[0:4]
    yyyy, month, day, hh, mm, ss = base.split('_')[4:]
    ss = ss.split('.')[0]
    eventid = yyyy + '_' + month + '_' + day + '_' + hh + '_' + mm + '_' + ss
    #make a tuple of record, event, station so indices can be assigned
    index_matrix[i] = [base, eventidlist.index(eventid), stationlist.index(station)]

print(eventidlist[0])
print(stationlist)

I = len(eventidlist)#events
J = len(stationlist)#stations
K = len(record_path)#records
K = rows

print 'Number of events: ', I, ' Number of stations: ', J
print 'Number of rows (records): ', K, ' Number of cols (events+stations): ', I+J

#make the G matrix of 1s and 0s and R matrix of records
G1 = np.zeros((K,I))
G2 = np.zeros((K,J))
R = np.zeros((K,1))

for k in range(K):#for all records
    G1[k][index_matrix[k][1]] = 1 #record row, eventid col
    G2[k][index_matrix[k][2]] = 1 #record row, station col

G = np.concatenate((G1,G2), axis = 1)

print(G)

R = np.zeros((K,F_bins))
#populate R matrix
for k in range(K):#for all records
    #each row is a record, and col is a frequency band
    #set row equal to the that spectral array
    #here we take the log
    R[k:,] = np.log10(record_spec[k])
    
m1 = np.zeros((I+J, F_bins))

#do the inversion for each freq
for f in range(F_bins):
    d = R[:,f]#record for given frequency col
    dT = d.T
    G_inv = np.linalg.pinv(G)
    m1[:,f] = np.dot(G_inv,dT)

#now we need to turn log power spectra into just the velocity spectra
#should still be cm/s
#####################################################################
m1 = 10**(m1)
m1 = np.sqrt(m1)

print(m1.shape)

#now split m into an event matrix and a station matrix
event = m1[0:I,:] #take first I rows
station = m1[I:I+J,:]
print(event.shape, station.shape)

#now event has the dim of events x frequency bins

#write an output file for each event and station
#outfile_path = boxpath + '/secondo'
outfile_path = '/Users/escuser/project/boxes/secondo_all'


for i in range(I):#for each event
    #make each row into an array
    amp = event[i,:]
    outfile = open(outfile_path + '/' + eventidlist[i] + '.out', 'w')
    out = (np.array([freq_list, amp])).T
    outfile.write('#freq_bins \t vel_spec_NE_cm \n')
    np.savetxt(outfile, out, fmt=['%E', '%E'], delimiter='\t')
    outfile.close()
    
#    plt.figure(figsize = (10,8))
#    plt.loglog(freq_list, amp, color='cornflowerblue')
#    plt.title('event: ' + eventidlist[i])
#    plt.xlabel('frequency (Hz)')
#    plt.ylabel('velocity spectrum cm/s')
#    plt.show()
#    plt.savefig(outfile_path + '/' + eventidlist[i] + '.png')

#    
for i in range(J):#for each station
    #make each row into an array
    amp = station[i,:]
    outfile = open(outfile_path + '/' + stationlist[i] + '.out', 'w')
    out = (np.array([freq_list, amp])).T
    outfile.write('#freq_bins \t vel_spec_NE_cm \n')
    np.savetxt(outfile, out, fmt=['%E', '%E'], delimiter='\t')
    outfile.close()
    
    
#    plt.figure(figsize = (10,8))
#    plt.loglog(freq_list, amp, color='cornflowerblue')
#    plt.title('station: ' + stationlist[i])
#    plt.xlabel('frequency (Hz)')
#    plt.ylabel('velocity spectrum cm/s')
#    plt.show()
#    plt.savefig(outfile_path + '/' + stationlist[i] + '.png')


########################################################################
constraint = 'Event_'

#long period spectral level
u0 = 1 #displacement, units of cm
ml = 2.5 #moment magnitude
M0 = 10**(1.5*(0.754*ml + 11.584))# Moment, dyne cm, calculated from Ross 2016
sig = 5 #stress drop mega pascals, 1 Pa = 1 dyne/cm2
B = 3500 #cm/s
fc = B*((sig)/(8.47*M0))**3
Brune = (2.*np.pi*freq_list*u0)/(1+(freq_list/fc)**2)

print(Brune)


constraint_file = '/Users/escuser/project/boxes/secondo_all/' + constraint + '.out'
data = np.genfromtxt(constraint_file)
con_spec = data.T[1] #second col
#
#for i in range(I):#for each event
#    #make each row into an array
#    #constrained with station so mult by station spec
#    amp = event[i,:]*con_spec
#    outfile = open(outfile_path + '/' + eventidlist[i] + '.out', 'w')
#    out = (np.array([freq_list, amp])).T
#    outfile.write('#freq_bins \t vel_spec_NE_cm \n')
#    np.savetxt(outfile, out, fmt=['%E', '%E'], delimiter='\t')
#    outfile.close()
#    
##    plt.figure(figsize = (14,8))
##    plt.loglog(freq_list, amp, color='cornflowerblue')
##    plt.title('event: ' + eventidlist[i])
##    plt.xlabel('frequency (Hz)')
##    plt.ylabel('velocity spectrum cm/s')
##    plt.grid()
##    plt.savefig(outfile_path + '/' + eventidlist[i] + '.png')
##    
#for i in range(J):#for each station
#    #make each row into an array
#    amp = (station[i,:])/con_spec
#    print(amp)
#    outfile = open(outfile_path + '/' + stationlist[i] + '.out', 'w')
#    out = (np.array([freq_list, amp])).T
#    outfile.write('#freq_bins \t vel_spec_NE_cm \n')
#    np.savetxt(outfile, out, fmt=['%E', '%E'], delimiter='\t')
#    outfile.close()
#    
#    
##    plt.figure(figsize = (14,8))
##    plt.loglog(freq_list, amp, color='cornflowerblue')
##    plt.title('station: ' + stationlist[i])
##    plt.ylim(0.1, 10)
##    plt.xlabel('frequency (Hz)')
##    plt.ylabel('velocity spectrum cm/s')
##    plt.grid()
##    plt.savefig(outfile_path + '/' + stationlist[i] + '.png')
#
#
#
































