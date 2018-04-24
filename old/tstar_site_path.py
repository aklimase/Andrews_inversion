#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  3 13:03:33 2017

@author: alexisklimasewski
"""

import numpy as np
import glob
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.patches as mpatches
import os.path as path
from obspy import read  
import dread

plt.rcParams["font.family"] = "serif"

mpl.rcParams.update({'font.size': 22})

#indices of frequency limits for inversion
f1 = 40
f2 = 48

top_dir = '/Volumes/USGS_Data/project'

boxpath = top_dir + '/boxes/'
#box = 'Imperial_Valley_SWS_ERR'
box_list = [('Imperial_Valley', 'PFO_TPFO_PMD'), ('Imperial_Valley', 'SWS_ERR'), ('Riverside', 'FRD_RDM'), ('Salton_Trough', 'SWS_ERR'), ('all_paths', '')]
box = box_list[1][0] + '_' + box_list[1][1]
residual_dir =  top_dir + '/boxes/all_paths/residual_spectra/' + box + '/'
site_dir = top_dir + '/boxes/all_paths/secondo_constrained_Brune/'

#list of files in the directory
record_list = glob.glob(residual_dir + '*.out')


I = len(record_list)
J = 50

residual = np.zeros((I, J))
decay_fit = np.zeros((I, f2-f1))

print 'number of records: ', I

G = np.zeros((I*J, I))
d = np.zeros((I*J, 1))

fig = plt.figure(figsize = (20,15))
plt.xscale('log')
plt.yscale('log')
plt.xlim(3, 40)
plt.ylim(0.001,100)
plt.ylabel(r'$\mathrm{residual_site spectra}$')
plt.xlabel(r'$\mathrm{frequency(Hz)}$')
plt.tick_params(axis='both', which='major', labelsize=15)
plt.tick_params(axis='both', which='both', length = 8, width = 2)

C = np.zeros((I))

for i in range(I):
    data = np.genfromtxt(record_list[i], comments = '#', dtype = float)
    freq = data.T[0][f1:f2]
    Amplitude = np.exp(data.T[1][f1:f2]) #data from log to lin
    residual[i:,] = np.exp(data.T[1])#residual in matrix, lin
    
    ########################################################################
    #now find the station where recorded and add site spectra on to the path
    station = record_list[i].split('/')[-1].split('_')[1]
    #station file path, open and add to amplitude
    
    
    A0 = Amplitude[0]#lin
    C[i] = A0#lin
    plt.plot(freq, Amplitude, color = 'blue', alpha= 0.2)#y = linear A
    plt.hold(True)
    for j in range(len(freq)):
        #print(Amplitude[j], np.log(Amplitude[j]))
        d[i*J+j][0] = np.log(Amplitude[j]/A0)#lin to log
        G[i*J+j][i] = -1*np.pi*j


G_inv = np.linalg.pinv(G, rcond=1e-10)

m = np.dot(G_inv,d)
print 'tstar array: ', m
print 'frequency list for inversion: ', freq

tstar = sorted(m, key=float)
tstar = np.concatenate(tstar)

for i in range(I):
    t = tstar[i]
    data = np.genfromtxt(record_list[i], comments = '#', dtype = float)
    freq = data.T[0][f1:f2]#linear A
    decay = C[i]*np.exp(-np.pi*t*freq)
    decay_fit[i:,] = decay
    plt.plot(freq, decay, color = 'black', alpha =0.2)
    plt.hold(True)

##############plot residual mean on top
mean = np.mean(residual, axis = 0)[f1:f2]
std = np.std(residual, axis = 0)[f1:f2]
plt.errorbar(freq, mean, yerr = std, zorder = 2000, color = 'red', elinewidth=2, capsize = 10, markeredgewidth=2, fmt='o', markersize=8)
plt.plot(freq, mean, color = 'red')
################plot t* decay mean on top

#remove outliers first 5 lowest tstar values
mean = np.mean(decay_fit[10:], axis = 0)
std = np.std(decay_fit[10:], axis = 0)
plt.errorbar(freq, mean, yerr = std, zorder = 2000, color = 'green', elinewidth=2, capsize = 10, markeredgewidth=2, fmt='o', markersize=8)
plt.plot(freq, mean, color = 'green')
##################

tstar_mean = np.mean(tstar)
print 'mean tstar:', tstar_mean 
#number of negative tstars
negtstar = []
for i in range(len(tstar)):
    if tstar[i]<0:
        negtstar.append(tstar[i])
frac_negtstar = 1.0*len(negtstar)/len(tstar)
print 'percentage of negative tstar values:', frac_negtstar*100

blue_patch = mpatches.Patch(color='blue', label='residual')
black_patch = mpatches.Patch(color='black', label='t* fit')
red_patch = mpatches.Patch(color='red', label='residual average')
green_patch = mpatches.Patch(color='green', label='t* fit average')
plt.legend(handles=[blue_patch, black_patch, red_patch, green_patch], prop={'size': 16})
############

plt.show()
#plt.savefig(top_dir + '/boxes/all_paths/residual_plots/t*_' + box + '.png')

################ tstar histogram
avg = np.mean(tstar)
sig = np.std(tstar)

fig = plt.figure(figsize = (20,15))
# the histogram of the data
plt.hist(tstar, 50, facecolor='green', alpha=0.75)

plt.xlabel(r'$\mathrm{t^*\ value}$')
plt.ylabel(r'$\mathrm{Number\ of\ Paths}$')
plt.title('Histogram of t* ' + r'$\mathrm{\mu}$' + ' =' + str(round(avg, 3)))
plt.grid(True)

plt.show()
plt.savefig(top_dir + '/boxes/all_paths/residual_plots/t*_histo_' + box + '.png')

#plot tstar vs distance
fig = plt.figure(figsize = (20,15))
plt.title(r'$\mathrm{residual\ t*\ vs.\ path\ distance}$')
#plt.xscale('log')
#plt.yscale('log')
#plt.xlim(3, 40)
#plt.ylim(0.001,100)
plt.ylabel(r'$\mathrm{t*}$')
plt.xlabel(r'$\mathrm{distance (km)}$')
plt.tick_params(axis='both', which='major', labelsize=15)
plt.tick_params(axis='both', which='both', length = 8, width = 2)

for i in range(I):
    record = (record_list[i].split('/')[-1])
    base = path.basename(record)
    network, station, channel, loc = base.split('_')[0:4]
    yyyy, month, day, hh, mm, ss = base.split('_')[4:]
    ss = ss.split('.')[0]
    eventid = yyyy + '_' + month + '_' + day + '_' + hh + '_' + mm + '_' + ss
    #read in uncorrected data for header info
    raw_file = boxpath + 'all_paths/uncorrected/Event_'+ eventid + '/' + network + '_' + station + '_HHN_' + loc + '_' + eventid + '.SAC'
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
    plt.scatter(dist, m[i])
plt.show()
#plt.savefig(top_dir + '/boxes/all_paths/residual_plots/t*vsdist_' + box + '.png')


