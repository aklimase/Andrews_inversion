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
from obspy import read    
import dread
import matplotlib.pyplot as plt
import time

#read in the cut and corrected spectra = records
#records are in cm/s
#should be binned frequencies and amplitudes
#make list of records and the corresponding events and stations

top_dir = '/Volumes/USGS_Data/project'

##############################################################################
#read in Joe's event list
#Joe_ev = top_dir + '/catalogs/Joe_evlist.txt'
#events = np.genfromtxt(Joe_ev, dtype = None, comments = '#', delimiter = None, names = True, encoding = None)
#ev_list = events['evid']

box = 'all_paths'

boxpath = top_dir + '/boxes/' + box
record_path = glob.glob(boxpath + '/record_spectra_meters/Event_2011*/*.out')
#

#record_path =[]
#for ev in ev_list:
#    record_path.extend(glob.glob(boxpath + '/record_spectra_meters/Event_' + ev + '/*.out'))
###         if not path.basename(fn).startswith('CI_ERR')])
###############################################################################

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
record_dist = []

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
    raw_file = boxpath + '/uncorrected/Event_'+ eventid + '/' + network + '_' + station + '_HHN_' + loc + '_' + eventid + '.SAC'
    stream = read(raw_file)
    tr = stream[0]
    
#    print('doing record: ' +  base)
    
    evlon =  tr.stats.sac.evlo #deg
    evlat =  tr.stats.sac.evla #deg
    evdepth = tr.stats.sac.evdp #km
    stlon = tr.stats.sac.stlo #deg
    stlat = tr.stats.sac.stla #deg
    stdepth = tr.stats.sac.stdp #km
    
    #find distance between event and station
    dist =  dread.compute_rrup(evlon, evlat, evdepth, stlon, stlat, stdepth) #in km
    #km to cm
    dist = dist*1000.
    record_dist.append(dist)
    #read in file
    data = np.genfromtxt(record_path[i], dtype = float, comments = '#', delimiter = None, usecols = (0,1))#only read in first two cols
    record_freq.append(data[:,0])
    #data is NE spectra
    #square for power spectra
    record_spec.append((data[:,1]*dist)**2.0)
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

print(freq_list)

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

C = np.ones((K,1))
G = np.concatenate((G1,G2, C), axis = 1)

print(G)

R = np.zeros((K,F_bins))
#populate R matrix
for k in range(K):#for all records
    #each row is a record, and col is a frequency band
    #set row equal to the that spectral array
    #here we take the log
    R[k:,] = np.log10(record_spec[k])
    ################################## solve for a distance dependent C
    C[k,0] = record_dist[k]
    
m1 = np.zeros((I+J+1, F_bins))#add in 1 for C

#do the inversion for each freq
for f in range(F_bins):
    t1 = time.time()
    d = R[:,f]#record for given frequency col
    dT = d.T
    print 'inverting for frequency: ', f, freq_list[f]
    G_inv = np.linalg.pinv(G, rcond=1e-13)
#    G_inv = scipy.linalg.pinv(G)
    m1[:,f] = np.dot(G_inv,dT)
    t2 = time.time()
    print 'time for inversion: (min) ', round((t2-t1)/60., 4)

#now we need to turn log power spectra into just the velocity spectra
#should still be cm/s
print(m1.shape)

#now split m into an event matrix and a station matrix
event = m1[0:I,:] #take first I rows
station = m1[I:I+J,:]
C = m1[I+J:I+J+1,:]
print(event.shape, station.shape, C.shape)

#now event has the dim of events x frequency bins

#write an output file for each event and station
outfile_path = boxpath + '/secondo_2011_Q'


for i in range(I):#for each event
    #make each row into an array
    amp = np.sqrt(np.power(10.0, event[i,:]))
    outfile = open(outfile_path + '/' + eventidlist[i] + '.out', 'w')
    out = (np.array([freq_list, amp])).T
    outfile.write('#freq_bins \t vel_spec_NE_m \n')
    np.savetxt(outfile, out, fmt=['%E', '%E'], delimiter='\t')
    outfile.close()

#    
for i in range(J):#for each station
    #make each row into an array
    amp = station[i,:]
#    print(stationlist[i], amp)
    amp = np.sqrt(np.power(10.0, station[i,:]))
    outfile = open(outfile_path + '/' + stationlist[i] + '.out', 'w')
    out = (np.array([freq_list, amp])).T
    outfile.write('#freq_bins \t vel_spec_NE_m \n')
    np.savetxt(outfile, out, fmt=['%E', '%E'], delimiter='\t')
    outfile.close()
    


#for Q
Beta = 3500. #m/s
#make each row into an array
#Q_f = np.log10(np.e)*((-2.0*np.pi*freq_list)/(Beta))*(1./C[0,:])
print(C[0])

####
outfile = open(outfile_path + '/' + 'C.out', 'w')
out = (np.array([freq_list, C[0]])).T
outfile.write('#freq_bins \t Bias \n')
np.savetxt(outfile, out, fmt=['%E', '%E'], delimiter='\t')
outfile.close()

###plot C
data = np.genfromtxt(outfile_path + '/' + 'C.out', dtype = float, comments = '#', delimiter = None, usecols = (0,1), encoding = None)#only read in first two cols
freq = data[:,0][30:48]
print freq
##these record spectra are in cm
C = (data[:,1])[30:48]


for i in range(len(freq)):
    print (np.log10(np.e)*(-2*np.pi*freq[i]/Beta))/C[i]
#
#G = np.zeros((len(freq), 2))
#d = np.zeros((len(freq), 1))
#
#Beta = 3500. #m/s
#
#for j in range(len(freq)):
#            #print(Amplitude[j], np.log(Amplitude[j]))
#    d[j][0] = (C[j])#lin to log
#    G[j][0] = (-1*np.pi*(freq[j]))/Beta
#    G[j][1] = 1
#    
#G_inv = np.linalg.pinv(G, rcond=1e-10)
#m = np.dot(G_inv,d)
#Q = 1./float(m[0])
#C0 = np.exp(float(m[1]))
#
#
#fig = plt.figure(figsize = (10,7))
#plt.ylabel('velocity amplitude (cm)', fontsize = 16)
#plt.xlim(0.5,70)
##plt.ylim(0.001,100)
#plt.semilogx(freq, C)
#
##decay = C[0]*np.exp(-np.pi*(freq)/(Q*Beta))
##plt.semilogx(freq, decay)
#
#plt.grid()
#plt.xlabel('frequency (hz)', fontsize = 16)
#plt.title('Constant bias')
#plt.tick_params(axis='both', which='major', labelsize=15)
#plt.tick_params(axis='both', which='both', length = 5, width = 1)
#plt.show()
#























