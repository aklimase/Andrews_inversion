#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  2 15:47:34 2017

@author: alexisklimasewski

input: reads in all of the site output files from secondo


"""

import numpy as np
import glob
import matplotlib.pyplot as plt
import os.path as path
import matplotlib.lines as mlines
#import mpl_defaults
#plot_defaults()


#indices of frequency limits for inversion
#f1 = 30 #30
#f2 = 48 #48

f1 = 28
f2 = 69

###########################

#############################
#directory of secondo event and site output files
#secondo_dir = 'secondo_constrained_2010_04_14_21_36_34'
secondo_dir = 'secondo_rebin3_constrained_2010_05_25_19_49_51'


top_dir = '/Volumes/USGS_Data/project/'
site_dir =  top_dir + 'boxes/all_paths/' + secondo_dir
#list of files in the directory
site_list = glob.glob(site_dir + '/[!2]*.out')

#############################
#directory of secondo event and site output files
#secondo_dir = 'secondo_meters_constrained_Brune'
#file name our output file of f*values
#outfilename = top_dir + 'tstar/site/secondo_constrained_2010_04_14_21_36_34/tstar_site.out'
#spec_levels_file =  top_dir + 'tstar/site/secondo_constrained_2010_04_14_21_36_34/spec_levels.out'

outfilename = top_dir + 'tstar/site/secondo_rebin3_constrained_2010_05_25_19_49_51/tstar_site.out'
spec_levels_file =  top_dir + 'tstar/site/secondo_rebin3_constrained_2010_05_25_19_49_51/spec_levels.out'


def calc_site_tstar(f1,f2):
    I = len(site_list)
    J = 75
    residual = np.zeros((len(site_list), J))
    
    site = np.zeros((I, J))
    decay_fit = np.zeros((I, f2-f1))
    
#    low_level = []
#    med_level = []
#    high_level = []
    
    print 'number of sites: ', I
    
    G = np.zeros((J, 2))
    d = np.zeros((J, 1))
    
    plt.figure(figsize = (20,15))
    plt.title(r'Site spectra')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlim(0.1, 70)
    plt.ylim(0.0001,5)
    plt.ylabel(r'site  spectra', fontsize = 26)
    plt.xlabel(r'frequency(Hz)', fontsize = 26)
    plt.tick_params(axis='both', which='major', labelsize=20)
    plt.tick_params(axis='both', which='both', length = 8, width = 2)
    plt.grid(linestyle='--', linewidth=1)
    
    line1 = mlines.Line2D([], [], color='red',linestyle = '-', label='site spectra')
    line2 = mlines.Line2D([], [], color='red', linestyle = '--', label='t* fit')
    legend1 = plt.legend(handles=[line1,line2],  loc='upper right', prop={'size': 20})
    plt.gca().add_artist(legend1)
    
#    C = np.zeros((I))
    tstar = []
    A = []
    L2norm = []
    
    cmap = plt.get_cmap('jet')
    colors = cmap(np.linspace(0, 1.0, I))
    
    ###################################
    avg_l = []
    avg_m = []
    avg_h = []
    site = []

    #for each site, calculate the 1-6, 6-14, and 14-36 Hz spectral levels
    for j in range(I):
        data = np.genfromtxt(site_list[j], comments = '#', dtype = float)
        freq = data.T[0]
        Amplitude = (data.T[1]) #data from log to lin

        l = []
        m = []
        h = []
    
        for i in range(len(freq)):
            if 14.0<freq[i]<=36.0:
                print freq[i]
                h.append(Amplitude[i])
            if 6.0<freq[i]<=14.0:
                m.append(Amplitude[i])
            if 1.0<freq[i]<=6.0:
                l.append(Amplitude[i])
        site.append(path.basename(site_list[j]).split('.')[0])
        avg_l.append(np.mean(l))
        avg_m.append(np.mean(m))
        avg_h.append(np.mean(h))
        print "{:.3f}".format(np.mean(l)), '\t&', "{:.3f}".format(np.mean(m)), '\t&', "{:.3f}".format(np.mean(h))
        
    outfile = open(spec_levels_file, 'w')
    out = (np.array([site, avg_l, avg_m, avg_h]).T)
    outfile.write('#site \t 1-6Hz(m) \t 6-14Hz(m) \t 14-36Hz(m) \n')
    np.savetxt(outfile, out, fmt='%s', delimiter='\t')
    outfile.close()


    
    #for each site
    for i in range(I):
        data = np.genfromtxt(site_list[i], comments = '#', dtype = float)
        freq = data.T[0][f1:f2]
        Amplitude = (data.T[1][f1:f2]) #data from log to lin
        
  
        residual[i:,] = (data.T[1])#residual in matrix, lin
        
        l = site_list[i].split('/')[-1].split('.')[0]
        
        plt.plot(data.T[0], data.T[1], color = colors[i], linestyle = '-', lw = 2, alpha= 1, label = l)#y = linear A
        for j in range(len(freq)):
            #print(Amplitude[j], np.log(Amplitude[j]))
            d[j][0] = np.log(Amplitude[j])#lin to ln
            G[j][0] = -1*np.pi*(freq[j])
            G[j][1] = 1
    
        G_inv = np.linalg.pinv(G, rcond=1e-10)
        m = np.dot(G_inv,d)
        t = float(m[0])
        A0 = np.exp(float(m[1]))
        tstar.append(float(t))
        A.append(A0)
    
    print'frequency range: ', (freq[0], freq[-1])
    
    
    for i in range(I):
        t = tstar[i]
        data = np.genfromtxt(site_list[i], comments = '#', dtype = float)
        freq = data.T[0][f1:f2]#linear A
        decay = A[i]*np.exp(-np.pi*t*(freq))
        decay_fit[i:,] = decay
        
        plt.plot(freq, decay, color = colors[i], linestyle = '--', alpha =1, lw = 2)
        
        data = np.genfromtxt(site_list[i], comments = '#', dtype = float)
        freq = data.T[0][f1:f2]
        Amplitude = (data.T[1][f1:f2]) #data from log to lin
    
        L2 = np.sqrt(np.sum((np.log10(decay)-np.log10(Amplitude))**2))#avg for all frequencies
        L2norm.append(L2)
    
    print'L2 norm average: ', np.mean(L2norm)
    
    plt.legend( loc='lower left', prop={'size': 20})
    plt.show()
#    plt.savefig('/Volumes/USGS_Data/project/tstar/site/t*_site_set1.png')
    
    for i in range(len(site_list)):
        print path.basename(site_list[i]).split('.')[0], '&', "{:.3f}".format(tstar[i]), '&', "{:.3f}".format(avg_l[i]), '\t&', "{:.3f}".format(avg_m[i]), '\t&', "{:.3f}".format(np.mean(avg_h[i])), ' \\\\'
        
    outfile = open(outfilename, 'w')
    out = (np.array([[path.basename(site_list[i]).split('.')[0] for i in range(len(site_list))], tstar, A]).T)
    outfile.write('#site \t tstar(s) \t A0(m) \n')
    np.savetxt(outfile, out, fmt='%s', delimiter='\t')
    outfile.close()
    
    return freq[0], freq[-1], np.mean(L2norm)


calc_site_tstar(f1,f2)