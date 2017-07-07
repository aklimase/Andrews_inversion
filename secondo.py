#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 29 11:36:09 2017

@author: escuser

takes in the power horizontal velocity spectra
then takes the log
after the inversion, raise (1/2) 10^(log(spectra^2))
"""
import glob
import os.path as path
import numpy as np
import matplotlib.pyplot as plt


#read in the cut and corrected spectra = records
#should be binned frequencies and amplitudes
#make list of records and the corresponding events and stations

boxpath = '/Users/escuser/Documents/Alexis_Data/cut_sac_files'
record_path = glob.glob(boxpath + '/record_spectra/*.AZ.*.out')
print 'Number of records: ', len(record_path)

#read in all files to find networks and stations
stationlist = []
eventidlist = []
record_freq = []
record_spec = []

#for i in range(len(records)):
for i in range(len(record_path)):
    record = (record_path[i].split('/')[7])
    base = path.basename(record)
    eventid, network, station, channel, extn = base.split('.')
    #read in file
    data = np.genfromtxt(record_path[i], dtype = float, comments = '#', delimiter = None, usecols = (0,1))#only read in first two cols
    record_freq.append(data[:,0])
    #data is NE power spectra
    record_spec.append(data[:,1])
    #make a tuple of record, event, station so indices can be assigned
    #if network and station not part of the list yet add
    if station not in stationlist:
        stationlist.append(station)
    if eventid not in eventidlist:
        eventidlist.append(eventid)
        
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
    record = (record_path[i].split('/')[7])
    base = path.basename(record)
    eventid, network, station, channel, extn = base.split('.')
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
m2 = np.zeros((I+J, F_bins))

#do the inversion for each freq
for f in range(F_bins):
    d = R[:,f]#record for given frequency col
    dT = d.T
    G_inv = np.linalg.pinv(G)
    m1[:,f] = np.dot(G_inv,dT)
    
#    #now do the svd for comparison
#    u,s,v = np.linalg.svd(G)
#    c = np.dot(u.T,d.reshape((rows,1)))
#    S = np.pad(np.diag(s), ((0, K - (I+J)),(0,0)), mode = 'constant')#resize s so that it has dim of rec x stations + events
#    S = np.linalg.pinv(S)
##    S = 1./s
#    w = np.dot(S,c)
#    x = np.dot(v.T,w)
#    m2[:,f] = x.flatten()
    
#    m = np.linalg.lstsq


#now we need to turn log power spectra into just the velocity spectra
#####################################################################
m1 = (0.5)*10**(m1)

print(m1.shape)

#now split m into an event matrix and a station matrix
event = m1[0:I,:] #take first I rows
station = m1[I:I+J,:]
print(event.shape, station.shape)

#now event has the dim of events x frequency bins

#write an output file for each event and station
outfile_path = boxpath + '/event_site_spectra'

print(freq_list)

for i in range(I):#for each event
    #make each row into an array
    amp = event[i,:]
    outfile = open(outfile_path + '/' + eventidlist[i] + '.out', 'w')
    out = (np.array([freq_list, amp])).T
    outfile.write('#freq_bins \t binned_vel_spectra \n')
    np.savetxt(outfile, out, fmt=['%E', '%E'], delimiter='\t')
    outfile.close()
    
    plt.figure(figsize = (10,8))
    plt.loglog(freq_list, amp, color='cornflowerblue')
    plt.title('event: ' + eventidlist[i])
    plt.xlabel('frequency (Hz)')
    plt.ylabel('velocity spectrum')
    plt.show()
    
for i in range(J):#for each station
    #make each row into an array
    amp = station[i,:]
    outfile = open(outfile_path + '/' + stationlist[i] + '.out', 'w')
    out = (np.array([freq_list, amp])).T
    outfile.write('#freq_bins \t binned_vel_spectra \n')
    np.savetxt(outfile, out, fmt=['%E', '%E'], delimiter='\t')
    outfile.close()
    
    
    plt.figure(figsize = (10,8))
    plt.loglog(freq_list, amp, color='cornflowerblue')
    plt.title('station: ' + stationlist[i])
    plt.xlabel('frequency (Hz)')
    plt.ylabel('velocity spectrum')
    plt.show()








































