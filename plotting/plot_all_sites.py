#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 23 14:29:18 2018

@author: temp

plot up all of the site spectra on separate subplots
4x4 plots
"""

import matplotlib.pyplot as plt
import numpy as np
import glob
import matplotlib as mpl
#import mpl_defaults
#plot_defaults()


plt.style.use("classic")
mpl.rcParams['figure.subplot.hspace'] = 0.0
mpl.rcParams['figure.subplot.wspace'] = 0.0
mpl.rcParams['font.size'] = 20


constraint = '2010_05_25_19_49_51'

top_dir = '/Volumes/USGS_Data/project/'
#st_dir = top_dir + 'boxes/all_paths/old_inversion_runs/secondo_constrained_Brune'
st_dir = top_dir + 'boxes/all_paths/secondo_rebin3_constrained_' + constraint

station_files = glob.glob(st_dir + '/[!2]*.out')
#t_file = top_dir + 'tstar/site/old_inversion_cm/tstar_site.out'
t_file = top_dir + 'tstar/site/secondo_rebin3_constrained_' + constraint + '/tstar_site.out'



#f1 = 30
#f2 = 48

f1 = 28
f2 = 69

def main():
    st_list = [station_files[i].split('/')[7].split('.')[0] for i in range(len(station_files))]
    
    station_data = np.genfromtxt(st_dir + '/' + st_list[0] + '.out', dtype = float, comments = '#', delimiter = None, usecols = (0,1), encoding = None)
    f = station_data[:,0]
    
    data_t = np.genfromtxt(t_file, dtype = None, comments = '#', delimiter = None, encoding = None, names = True)
    t_st_list = data_t['site']
    t_list = data_t['tstars']
    A_list = data_t['A0m']
    
    nrows = 4
    ncols = 4
    fig, axes = plt.subplots(nrows, ncols, figsize = (20,12))
#    plt.tight_layout(pad=1, w_pad=1, h_pad=1.0)
    plt.subplots_adjust(wspace=0, hspace=0)
    fig.text(0.47, 0.05, 'frequency (Hz)', va='center')
    fig.text(0.06, 0.5, 'Velocity amplitude (m)', va='center', rotation = 'vertical')

    spec_list = []
    decay_list = []
    L2_norm_list = []


    for st in st_list:
        station_data = np.genfromtxt(st_dir + '/' + st + '.out', dtype = float, comments = '#', delimiter = None, usecols = (0,1), encoding = None)
        station_spec = station_data[:,1]
#        f = station_data[:,0]
        spec_list.append(station_spec)

        
        ind = list(t_st_list).index(st)
        freq_decay = station_data[:,0][f1:f2]#linear A
        decay_list.append(A_list[ind]*np.exp(-np.pi*t_list[ind]*freq_decay))
#        print np.log10(station_spec[f1:f2])
#        print np.log10(A_list[ind]*np.exp(-np.pi*t_list[ind]*freq_decay))
        L2_norm_list.append(np.mean((np.log10(station_spec[f1:f2]) - np.log10(A_list[ind]*np.exp(-np.pi*t_list[ind]*freq_decay)))**2.))

    print st_list
    print L2_norm_list
    #loop throught the 16 grid locations
    for i, ax in enumerate(fig.axes):
        print i
        plot_spectra(ax, x1 = f, y1 = spec_list[i], x2 = freq_decay, y2 = decay_list[i], name = st_list[i], t = "{:.3f}".format(t_list[i]), A = "{:.3f}".format(A_list[i]), L2 = "{:.3f}".format(L2_norm_list[i]))
      
#    plt.xlabel('frequency')
    
#    plt.savefig(top_dir + 'tstar/site/old_inversion_cm/tstar_allsites.png')
    plt.savefig(top_dir + 'tstar/site/secondo_rebin3_constrained_' + constraint + '/tstar_allsites.png')
    plt.show()

def plot_spectra(ax, x1, y1, x2, y2, name, t, A, L2):
    ax.loglog(x1,y1, c = 'C0')
    ax.loglog(x2,y2, ls = '--', c = 'C2')
#    ax.title.set_text(name)
    ax.grid(linestyle='--', linewidth=0.5, color = 'gray')
    ax.tick_params(axis='both', which='major', labelsize=18)
    ax.tick_params(axis='both', which='both', length = 5, width = 1)
    ax.text(0.55, 0.01, name)
    ax.text(0.55, 0.003, 'A0: ' + str(A))
    ax.text(0.55, 0.001, r'$\kappa$: ' + str(t))
    ax.text(4, 0.001, 'L2: ' + str(L2))
    ax.set_xlim(0.5,70)
    ax.set_ylim(0.0005,3)#1.4


main()

##make a spectra plotting function

