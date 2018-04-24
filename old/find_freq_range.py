#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 16 13:57:10 2017

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
import scipy

plt.rcParams["font.family"] = "serif"

mpl.rcParams.update({'font.size': 22})

#indices of frequency limits for inversion
#f1 = 38
#f2 = 46

top_dir = '/Volumes/USGS_Data/project'

boxpath = top_dir + '/boxes/'
box_list = [('Imperial_Valley', 'PFO_TPFO_PMD'), ('Imperial_Valley', 'SWS_ERR'), ('Riverside', 'FRD_RDM'), ('Salton_Trough', 'SWS_ERR'), ('all_paths', '')]
box = box_list[3][0] + '_' + box_list[3][1]
residual_dir =  top_dir + '/boxes/all_paths/residual_spectra/' + box + '/'

#list of files in the directory
record_list = glob.glob(residual_dir + '*.out')

#######################################
#find the best frequenies between
i1 = 34 #0.82 hz
i2 = 47 #42.2 hz

L2_norm_pair = []

data = np.genfromtxt(record_list[0], comments = '#', dtype = float)
freq_all = data.T[0]

pair_list = []
#make a list of all possible pairs of indices 24-48 at least 5 long
index = np.arange(i1,i2+1,1)
for ind1 in index:
    for ind2 in range(ind1+5, index[-1] +1):
        pair = (ind1, ind2)
        pair_list.append(pair)

for k in range(len(pair_list)):
    f1 = pair_list[k][0]
    f2 = pair_list[k][1]
    
    print 'doing indices: ', f1, f2

    #I = len(record_list)
    I=1
    J = 50
    
    residual = np.zeros((len(record_list), J))
    decay_fit = np.zeros((I, f2-f1))
    
    G = np.zeros((I*J, I))
    d = np.zeros((I*J, 1))
    
    C = np.zeros(len(record_list))
    tstar = []
    
    L2norm = 0
    
    for i in range(len(record_list)):
        data = np.genfromtxt(record_list[i], comments = '#', dtype = float)
        freq = data.T[0][f1:f2]
        Amplitude = data.T[1][f1:f2] #data from log to lin
        residual[i:,] = data.T[1]#residual in matrix, lin
        plt.hold(True)
        A0 = np.max(Amplitude)#lin
        C[i] = A0#lin
        plt.plot(data.T[0], np.log(data.T[1]), color = 'green', alpha= 0.2)#y = linear A
        plt.hold(True)
        for j in range(len(freq)):
            d[j][0] = np.log(Amplitude[j]/A0)#lin to log
            G[j][0] = -1*np.pi*(freq[j])
    
        G_inv = np.linalg.pinv(G, rcond=1e-10)
        m = float(np.dot(G_inv,d))
        tstar.append(m)
    
        t = m
        decay = C[i]*np.exp(-np.pi*t*(freq))
        decay_fit[i:,] = decay
        plt.plot(freq, np.log(decay), color = 'blue', alpha =0.2)
        plt.hold(True)
        ##average
        L2 = np.mean(np.sqrt((decay-Amplitude)**2)) #avg for all frequencies
        L2norm += L2 #sum all the L2 norms for the pair

    L2_norm_pair.append(L2norm)
    
#print(L2_norm_pair)
index = L2_norm_pair.index(min(L2_norm_pair))
print(len(L2_norm_pair))
print(len(pair_list))
print 'frequencies with smalles L2 norm: ', pair_list[index]