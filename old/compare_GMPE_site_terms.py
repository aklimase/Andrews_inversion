#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 27 14:47:14 2017

@author: alexisklimasewski

"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy.stats import pearsonr
import statsmodels
import statsmodels.stats.power
import glob
from matplotlib.gridspec import GridSpec
#import mpl_defaults
#plot_defaults()
mpl.rcParams['figure.subplot.right'] = 0.95

#####################################################
tstar_file = 'tstar_site.out'
GMPE_site_resid_file = 'site_terms.txt'
vs_30_file = 'v3anza2013_vs30.txt'
secondo_dir = 'secondo_meters_constrained_Brune'
spec_levels = 'spec_levels.out'
##outfile
#site_resid_amp = 'GMPE_site_resid_vs_site_amp'
#site_resid_tstar = 'GMPE_site_reid_vs_site_tstar'
######################################################

top_dir = '/Volumes/USGS_Data/project/'

sitefile = top_dir + 'tstar/site/' + tstar_file
#tuple with site name and tstar value
sitet = np.genfromtxt(sitefile, delimiter='\t', names = True, dtype = None)

#tuple with station, site residual, standard error
sitetermfile = top_dir + 'event_boxes/site/' + GMPE_site_resid_file
siteresid = np.genfromtxt(sitetermfile, delimiter='\t', names = True, dtype = None)

#tuple with station, lon, lat, vs30, method
vs30file = top_dir + 'event_boxes/' + vs_30_file
vs30 = np.genfromtxt(vs30file, delimiter='\t', names = True, dtype = None)

levels_file = top_dir + 'tstar/site/' + spec_levels
data = np.genfromtxt(vs30file, delimiter='\t', names = True, dtype = None)

##tupe with station, 1-6Hz amplification
spectra_paths = glob.glob(top_dir + 'boxes/all_paths/' + secondo_dir + '/[!2]*.out')
st_list = []
lowf_list = []
medf_list = []
highf_list = []
for st in spectra_paths:
    station_data = np.genfromtxt(st, dtype = float, comments = '#', delimiter = None, usecols = (0,1))
    station_spec = station_data[:,1]
    f_bins = station_data[:,0]
    l = []
    m = []
    h = []
    for i in range(len(f_bins)):
        if 15.0<f_bins[i]<36.0:
            h.append(station_spec[i])
        if 6.0<f_bins[i]<14.0:
            m.append(station_spec[i])
        if 1.0<f_bins[i]<6.0:
            l.append(station_spec[i])
    st_list.append(st.split('/')[7].split('.')[0])
    lowf_list.append(np.mean(l))
    medf_list.append(np.mean(m))
    highf_list.append(np.mean(h))
l_level = zip(st_list, lowf_list)
m_level = zip(st_list, medf_list)
h_level = zip(st_list, highf_list)


site = []
siteterm = []
sitevs30 = []
l = []
h = []
m = []

##go through and match with station name to get each spectral levels, site residual, and vs30
for i in range(len(l_level)):
    s = l_level[i][0]
    level1_6 = l_level[i][1]
    level6_14 = m_level[i][1]
    level14_36 = h_level[i][1]
    for j in range(len(siteresid)):
        if siteresid[j][0].split(':')[0] == s:
            r = siteresid[j][1]
            for k in range(len(vs30)):
                if vs30[k][0] == s:
                    v = vs30[k][3]
                    site.append(s)
                    l.append(level1_6)
                    m.append(level6_14)
                    h.append(level14_36)
                    siteterm.append(r)
                    sitevs30.append(v)
                    

def plot_spectral_levels(yvar, yname, color):
    #yvar either vs30 or GMPE site residual
    #cbar is True or False for coloring by vs30
    fig = plt.figure(figsize = (18,6))
    fig.subplots_adjust(wspace = 0)
    gs = GridSpec(1,3)
    
    cmap = mpl.cm.viridis
    norm = mpl.colors.Normalize(vmin=min(sitevs30), vmax=max(sitevs30))
    s_m = mpl.cm.ScalarMappable(cmap = cmap, norm=norm)
    s_m.set_array([])
    
    #1-6Hz subplot
    ax1 = plt.subplot(gs.new_subplotspec((0, 0), colspan=1))
    if color == True:
        plt.scatter(l, yvar, marker = 'o', s = 60, c = s_m.to_rgba(sitevs30), edgecolors = s_m.to_rgba(sitevs30))
    else:
        plt.scatter(l, yvar, marker = 'o', s = 60, c = 'blue')
    for i in range(len(l)):#default 2,5
       plt.annotate(site[i], xy = (l[i], yvar[i]), xytext = (2, 2), textcoords = 'offset points', fontsize = 20)
    plt.ylabel(yname)
    rval, pval = pearsonr(l, siteterm)
    power = statsmodels.stats.power.tt_solve_power(effect_size = rval, nobs = len(l), alpha = 0.05)
    label1 = 'Pearson R : ' + str(round(rval,4))
    label2 = 'power : ' + str(round(power,4))
    plt.annotate(label1, xy=(0.2, 0.16), xycoords='figure fraction', fontsize = 16)
    plt.annotate(label2, xy=(0.2, 0.13), xycoords='figure fraction', fontsize = 16)

    
    #6-14Hz subplot
    ax2 = plt.subplot(gs.new_subplotspec((0, 1), colspan=1), sharey = ax1)
    plt.setp(ax2.get_yticklabels(), visible = False)
    if color == True:
        plt.scatter(m, yvar, marker = 'o', s = 60,c = s_m.to_rgba(sitevs30), edgecolors = s_m.to_rgba(sitevs30))
    else:
         plt.scatter(m, yvar, marker = 'o', s = 60,c = 'blue')
    for i in range(len(m)):#default 2,5
       plt.annotate(site[i], xy = (m[i], yvar[i]), xytext = (2, 2), textcoords = 'offset points', fontsize = 20)
    rval, pval = pearsonr(m, siteterm)
    power = statsmodels.stats.power.tt_solve_power(effect_size = rval, nobs = len(m), alpha = 0.05)
    label1 = 'Pearson R : ' + str(round(rval,4))
    label2 = 'power : ' + str(round(power,4))
    plt.annotate(label1, xy=(0.5, 0.16), xycoords='figure fraction', fontsize = 16)
    plt.annotate(label2, xy=(0.5, 0.13), xycoords='figure fraction', fontsize = 16)
    
    #14-36Hz subplot
    ax3 = plt.subplot(gs.new_subplotspec((0, 2), colspan=1), sharey = ax1)
    plt.setp(ax3.get_yticklabels(), visible = False)
    if color == True:
        plt.scatter(h, yvar, marker = 'o', s = 60,c = s_m.to_rgba(sitevs30), edgecolors = s_m.to_rgba(sitevs30))
    else:
        plt.scatter(h, yvar, marker = 'o', s = 60,c = 'blue')
    for i in range(len(h)):#default 2,5
       plt.annotate(site[i], xy = (h[i], yvar[i]), xytext = (2, 2), textcoords = 'offset points', fontsize = 20)
    rval, pval = pearsonr(h, siteterm)
    power = statsmodels.stats.power.tt_solve_power(effect_size = rval, nobs = len(h), alpha = 0.05)
    label1 = 'Pearson R : ' + str(round(rval,4))
    label2 = 'power : ' + str(round(power,4))
    plt.annotate(label1, xy=(0.8, 0.16), xycoords='figure fraction', fontsize = 16)
    plt.annotate(label2, xy=(0.8, 0.13), xycoords='figure fraction', fontsize = 16)
    
    if color ==True:
        ##add color bar to right
        cbar = plt.colorbar(s_m)
        cbar.set_label(ur"Vs30 (m/s)", fontsize = 18)#"$azimuth$ (\u00b0)"
        cbar.ax.tick_params(labelsize = 18)
    
    plt.show()
    plt.savefig(top_dir + 'tstar/site/' + 'spectral_amplitude_vs_' + yname + '.png')


def plot_tstar(yvar, yname, color):
    #yvar name either siteterm or sitevs30
    site = []
    siteterm = []
    sitevs30 = []
    tstar = []
    #GMPE site residual vs. t*
    ##go through and match with station name to get each tstar, site residual, and vs30
    for i in range(len(sitet)):
        s = sitet[i][0]
        t = sitet[i][1]
        for j in range(len(siteresid)):
            if siteresid[j][0].split(':')[0] == s:
                r = siteresid[j][1]
                for k in range(len(vs30)):
                    if vs30[k][0] == s:
                        v = vs30[k][3]
                        site.append(s)
                        siteterm.append(r)
                        sitevs30.append(v)
                        tstar.append(t)

    cmap = mpl.cm.viridis
    norm = mpl.colors.Normalize(vmin=min(sitevs30), vmax=max(sitevs30))
    s_m = mpl.cm.ScalarMappable(cmap = cmap, norm=norm)
    s_m.set_array([])

                
    fig = plt.figure(figsize = (20,16))
    if color == True:
        plt.scatter(tstar, yvar, marker = 'o', s = 60,c = s_m.to_rgba(sitevs30), edgecolors = s_m.to_rgba(sitevs30))
    else:
        plt.scatter(tstar, yvar, marker = 'o', s=60, c = 'blue', edgecolors = 'blue')
    for i in range(len(tstar)):#default 2,5
       plt.annotate(site[i], xy = (tstar[i], yvar[i]), xytext = (2, 2), textcoords = 'offset points', fontsize = 20)
    plt.xlabel('Site t*')
    plt.ylabel(yname)
    rval, pval = pearsonr(tstar, siteterm)
    power = statsmodels.stats.power.tt_solve_power(effect_size = rval, nobs = len(tstar), alpha = 0.05)
    label1 = 'Pearson R : ' + str(round(rval,4))
    label2 = 'power : ' + str(round(power,4))
    plt.annotate(label1, xy=(0.1, 0.16), xycoords='figure fraction', fontsize = 18)
    plt.annotate(label2, xy=(0.1, 0.13), xycoords='figure fraction', fontsize = 18)
    
    if color == True:
        ##add color bar to right
        cbar = plt.colorbar(s_m)
        cbar.set_label(ur"Vs30 (m/s)", fontsize = 18)
        cbar.ax.tick_params(labelsize = 18)
    
    plt.show()
    plt.savefig(top_dir + 'tstar/site/tstar_vs_' + yname + '.png')


plot_tstar(yvar = siteterm, yname = 'GMPE site residual', color = True)
plot_tstar(yvar = sitevs30, yname = 'Site Vs30', color = False)

plot_spectral_levels(yvar = siteterm, yname= 'GMPE site residual', color = True)
plot_spectral_levels(yvar = sitevs30, yname = 'Site Vs30', color = False)
plot_spectral_levels(yvar = tstar, yname = 't*', color = False)






