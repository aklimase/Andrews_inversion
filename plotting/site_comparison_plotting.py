#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 14 12:34:43 2018

@author: temp

plots site t* vs. Joe's site t*, Anderson 1991 site t*, and Kilb 2012 kappas
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy.stats import pearsonr
import statsmodels
import statsmodels.stats.power
import matplotlib.patches as mpatches
import glob
import mpl_defaults
plot_defaults()

#
#plt.rcParams["font.family"] = "serif"
#plt.rcParams['axes.axisbelow'] = True
#mpl.rcParams.update({'font.size': 22})

top_dir = '/Volumes/USGS_Data/project/'

#reading in data
######################
sitefile = top_dir + 'tstar/site/Joe_ev_tstar_site.out'
#tuple with site name and tstar value
sitet = np.genfromtxt(sitefile, delimiter='\t', names = True, dtype = None, encoding = None)

#tuple with station, site residual, standard error
sitetermfile = top_dir + 'event_boxes/site/site_terms.txt'
siteresid = np.genfromtxt(sitetermfile, delimiter='\t',names = True, dtype = None, encoding = None)

joesitefile = top_dir + 'event_boxes/site/Joe_anza_tstar.txt'
joesitet = np.genfromtxt(joesitefile, delimiter='\t', dtype = None, encoding = None)

jstat = (joesitet[1].split('[')[1]).split(']')
jstat = (jstat[0].replace("'","")).split(';')

###anderson kappa
andersonsitefile = top_dir + 'event_boxes/site/Anderson_1991_skappa.txt'
andsitek = np.genfromtxt(andersonsitefile, delimiter='\t',names = True, dtype = None, encoding = None)

#Kilb et al 2012 kappa
kilbfile = top_dir + 'event_boxes/site/Kilb_2012_table5.txt'
kilbk = np.genfromtxt(kilbfile, delimiter=' ',names = True, dtype = None, encoding = None)


for i in range(len(jstat)): #take LVA and LVA2 to be pretty much the same
    if jstat[i] == 'LVA':
        jstat[i] = 'LVA2'
        
for i in range(len(andsitek)): #take LVA and LVA2 to be pretty much the same
    if andsitek[i][1] == 'LVA':
        andsitek[i][1] = 'LVA2'
        

jtstar = (joesitet[3].split('[')[1]).split(']')
jtstar = jtstar[0].split(' ')
jsitet = zip(jstat, jtstar)

jtavgl = (joesitet[5].split('[')[1]).split(']')
jtavgl = jtavgl[0].split(' ')
jtavgl = zip(jstat, jtavgl)

jtavgh = (joesitet[7].split('[')[1]).split(']')
jtavgh = jtavgh[0].split(' ')
jtavgh = zip(jstat, jtavgh)
#######################################
spectra_paths = glob.glob(top_dir + 'boxes/all_paths/secondo_constrained_Brune/[!2]*.out')

st_list = []
l_list = []
h_list = []
for st in spectra_paths:
    station_data = np.genfromtxt(st, dtype = float, comments = '#', delimiter = None, usecols = (0,1))
    station_spec = station_data[:,1]
    f_bins = station_data[:,0]
    l = []
    h = []
    for i in range(len(f_bins)):
        if 1.0<f_bins[i]<6.0:
            l.append(station_spec[i])
        if 6.0<f_bins[i]<14.0:
            h.append(station_spec[i])
    st_list.append(st.split('/')[7].split('.')[0])
    l_list.append(np.mean(l))
    h_list.append(np.mean(h))
l_level = zip(st_list, l_list)
h_level = zip(st_list, h_list)

###############################################################################
#compare site tstar to Joe's site tstars
site = []
tstar = []
joe_tstar = []

for i in range(len(sitet)):
    s = sitet[i][0]
    t = sitet[i][1]
    for j in range(len(jsitet)):
        if jsitet[j][0] == s: #same station
            site.append(s)
            tstar.append(t)
            joe_tstar.append(float(jsitet[j][1]))
            
#relative change
median = np.median(joe_tstar)
j_rel = [(x-median)/median for x in joe_tstar]

median = np.median(tstar)
t_rel = [(x-median)/median for x in tstar]

j_var = zip(site, tstar, joe_tstar)
j_var_rel = zip(site, t_rel, j_rel)

#compare site spectral levels
site = []
t_l_level = []
joe_l_level = []

for i in range(len(l_level)):
    s = l_level[i][0]
    l = l_level[i][1]
    for j in range(len(jtavgl)):
        if jtavgl[j][0] == s: #same station
            site.append(s)
            t_l_level.append(l)
            joe_l_level.append(float(jtavgl[j][1]))
            
#relative change
median = np.median(t_l_level)
t_rel = [(x-median)/median for x in t_l_level]

median = np.median(joe_l_level)
j_rel = [(x-median)/median for x in joe_l_level]

l_level = zip(site, t_l_level, joe_l_level)
l_level_rel = zip(site, t_rel, j_rel)
            

#compare site spectral levels
site = []
t_h_level = []
joe_h_level = []

for i in range(len(h_level)):
    s = h_level[i][0]
    h = h_level[i][1]
    for j in range(len(jtavgh)):
        if jtavgh[j][0] == s: #same station
            site.append(s)
            t_h_level.append(h)
            joe_h_level.append(float(jtavgh[j][1]))
            
#relative change
median = np.median(t_h_level)
t_rel = [(x-median)/median for x in t_h_level]

median = np.median(joe_h_level)
j_rel = [(x-median)/median for x in joe_h_level]

h_level = zip(site, t_h_level, joe_h_level)
h_level_rel = zip(site, t_rel, j_rel)
            
#compare site tstar to Anderson 1991 site kappa
site = []
tstar = []
andersonk = []

for i in range(len(sitet)):
    s = sitet[i][0]
    t = sitet[i][1]
    for j in range(len(andsitek)):
        if andsitek[j][0] == s: #same station
            site.append(s)
            tstar.append(t)
            andersonk.append(andsitek[j][1]/1000)    
            
median = np.median(andersonk)
a_rel = [(x-median)/median for x in andersonk]

median = np.median(tstar)
t_rel = [(x-median)/median for x in tstar]         
            
anderson_var = zip(site, tstar, andersonk)
anderson_var_rel = zip(site, t_rel, a_rel)
            
site = []
tstar = []
channel = []
kah = []
kah_mean = []
ksfix = []
ksfix_mean = []

for i in range(len(sitet)):
    s = sitet[i][0]
    t = sitet[i][1]
    for j in range(len(kilbk)):
        if kilbk[j][0] == s: #same station
            site.append(s)
            tstar.append(t)
            channel.append(kilbk[j][1])
            kah.append((kilbk[j][2]/1000, kilbk[j+1][2]/1000))
            #average of E, N components
            kah_mean.append(np.mean([kilbk[j][2]/1000,kilbk[j+1][2]/1000]))
            ksfix.append((kilbk[j][3]/1000,kilbk[j+1][3]/1000))
            ksfix_mean.append(np.mean([kilbk[j][3]/1000,kilbk[j+1][3]/1000]))
            break   


median = np.median(kah_mean)
kah_rel = [(x-median)/median for x in kah_mean]

median = np.median(tstar)
t_rel = [(x-median)/median for x in tstar]         

median = np.median(ksfix_mean)
ksfix_rel = [(x-median)/median for x in ksfix_mean]

kilb_var_AH = zip(site, tstar, kah_mean)
kilb_var_sfix = zip(site, tstar, ksfix_mean)

kilb_rel_AH = zip(site, t_rel, kah_rel)
kilb_rel_sfix = zip(site, t_rel, ksfix_rel)
            
def plot_t(var_tuple, name):
    #var_tuple is a list of (station, tstar, othervariable)
    #name is name of the othervariable, string
    station = [str(l[0]) for l in var_tuple]
    x = [l[1] for l in var_tuple]
    y = [l[2] for l in var_tuple]
    
    rval, pval = pearsonr(x, y)
    power = statsmodels.stats.power.tt_solve_power(effect_size = rval, nobs = len(tstar), alpha = 0.05)

    label1 = 'Pearson R : ' + str(round(rval,4))
    label2 = 'power : ' + str(round(power,4))

    plt.figure(figsize = (15,15))

    plt.grid(ls = '--', lw = 1, zorder = 0)
    plt.scatter(x, y, s = 50)
#    plt.xlim(0,600)

    plt.xlabel(name + ' secondo events with Brune constraint')
    plt.ylabel(name + ' Joe\'s inversion')
    plt.tick_params(axis='both', which='major', labelsize=18)
    plt.tick_params(axis='both', which='both', length = 8, width = 2)

    for i in range(len(x)):
        if station[i] != 'KNW':
            plt.annotate(station[i], xy = (x[i], y[i]), xytext = (2, 5), textcoords = 'offset points', fontsize = 18)
        else:
            plt.annotate(station[i], xy = (x[i], y[i]), xytext = (-20, 5), textcoords = 'offset points', fontsize = 18)
    plt.annotate(label1, xy=(0.15, 0.83), xycoords='figure fraction', fontsize = 22)
    plt.annotate(label2, xy=(0.15, 0.80), xycoords='figure fraction', fontsize = 22)
    plt.savefig(top_dir + 'tstar/site/compare_Joe_ev_' + name + '.png' )
    plt.show()
#
#plot_t(j_var, 'Joe_site_tstar')
#plot_t(anderson_var, 'Anderson_1991_site_tstar')
#plot_t(kilb_var_AH, 'Kilb_Table_5_kappa_AH')
#plot_t(kilb_var_sfix, 'Kilb_Table_5_kappa_sfix')
#
###now lets look at relative change
#plot_t(j_var_rel, 'rel_diff_Joetstar')
#plot_t(anderson_var_rel, 'rel_diff_Andersontstar')
#plot_t(kilb_rel_AH, 'rel_diff_Kilb_Table_5_kappa_AH')
#plot_t(kilb_rel_sfix, 'rel_diff_Kilb_Table_5_kappa_sfix')
#
#plot_t(l_level_rel, 'rel_1-6_Hz_average_site_response')
#plot_t(h_level_rel, 'rel_6-14_Hz_average_site_response')
plot_t(l_level, '1-6_Hz_average_site_response')
plot_t(h_level, '6-14_Hz_average_site_response')
#
