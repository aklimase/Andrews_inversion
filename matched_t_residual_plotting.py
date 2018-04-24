#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 19 10:00:58 2018

@author: temp
"""

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib as mpl
from scipy.stats import pearsonr
import statsmodels.stats.power

#w,e,n,s
IV = [-115.75000 , -115.51000, 32.75000, 32.50000] #3a PFO, PMD and SWS, ERR
ST = [-115.75000, -115.20000, 33.40000, 33.00000]  #4 SWS,ERR
R = [-117.80000, -117.50000, 34.18000, 34.00000]#6 FRD,RDM


plt.rcParams["font.family"] = "serif"
mpl.rcParams.update({'font.size': 22})

sns.set_context("poster")
sns.set(rc={'axes.labelsize':20, 'xtick.labelsize':20, 'ytick.labelsize':20, "font.size":20})
sns.set(font_scale=2)

sns.set_style('whitegrid')

topdir = '/Volumes/USGS_Data/project/'

d = [20]#for each distance

for i in d:
    match_d = i
    
    f = topdir + 'boxes/all_paths/t*plots/pathterms/match_pathresid_' + str(match_d) +'km.txt'
    data = np.genfromtxt(f,delimiter='\t', dtype = None, names = True, encoding = None)
    t = data['tstar']
    r = data['path_residual']
    ev_lat = data['Lat_cat']
    ev_lon = data['Lon_cat']
    ev_dep = data['Depth_cat']
    station = data['station']
    
    
    rval, pval = pearsonr(t, r)
    h= sns.jointplot(x = np.asarray(t), y = np.asarray(r), kind="regression", space=0, size = 14, scatter_kws={'alpha':0.5})
    h.set_axis_labels('path t*', 'GMPE path residual', fontsize=20)
    
    # labels appear outside of plot area, so auto-adjust
    plt.tight_layout()
    plt.xlim(-0.06, 0.06)
    plt.ylim(-2,2.5)
    #plt.xlabel('path t*')
    #plt.ylabel('GMPE path residual')
    rval, pval = pearsonr(t, r)
    power = statsmodels.stats.power.tt_solve_power(effect_size = rval, nobs = len(t), alpha = 0.05)
    label2 = 'power ' +  ': ' + str(round(power,4))
    plt.annotate(label2, xy=(0.15, 0.75), xycoords='figure fraction', fontsize = 20)
    plt.savefig(topdir + 'boxes/all_paths/t*plots/pathterms/match_' + str(match_d) + 'km_scatter.png')
    #plt.xlim(-0.06, 0.06)
    #plt.ylim(-2,2.5)
    plt.show()
    
    
    h = sns.jointplot(x = np.asarray(t), y = np.asarray(r), kind="kde", color="g", space = 0, size = 14)
    h.set_axis_labels('path t*', 'GMPE path residual', fontsize=20)
    plt.tight_layout()
#    plt.xlim(-0.06, 0.06)
#    plt.ylim(-2,2.5)
    rval, pval = pearsonr(t, r)
    power = statsmodels.stats.power.tt_solve_power(effect_size = rval, nobs = len(t), alpha = 0.05)
    label2 = 'power ' +  ': ' + str(round(power,4))
    plt.annotate(label2, xy=(0.71, 0.75), xycoords='figure fraction', fontsize = 20)
    plt.savefig(topdir + 'boxes/all_paths/t*plots/pathterms/match_' + str(match_d) + 'km_contour.png')
    
    plt.show()
    

###############################################################################
##color boxes differently
#t_box = []
#r_box = []
#t_other = []
#r_other = []
#for i in range(len(t)):
#    if R[3] <= ev_lat[i] <= R[2] and R[0] <= ev_lon[i] <= R[1] and (station[i]=='FRD' or station[i] == 'RDM'):
#        t_box.append(t[i])
#        r_box.append(r[i])
#        print(ev_lat[i], ev_lon[i], station[i])
#    else:
#        t_other.append(t[i])
#        r_other.append(r[i])
#
#fig = plt.subplots(figsize = (15,15))
#plt.xlabel(r'$\mathrm{path \ t*}$')
#plt.ylabel(r'$\mathrm{GMPE \ path \ residual}$')
#plt.tick_params(axis='both', which='major', labelsize=18)
#plt.tick_params(axis='both', which='both', length = 8, width = 2)
#
#
#plt.scatter(t_other,r_other, color = 'blue', alpha = 0.4)
#plt.scatter(t_box,r_box, color = 'orange', alpha = 0.95)
#
#plt.title('GMPE path residual vs. path t*', fontsize = 22)
#
#plt.savefig(topdir + 'boxes/all_paths/t*plots/pathterms/match_' + str(match_d) + 'km_scatter_Riverside_FRD_RDM.png')
#
#plt.show()
