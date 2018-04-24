#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 27 14:47:14 2017

@author: alexisklimasewski

Imperial Valley = 3a
Riverside = 6
Salton Trough = 4
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy.stats import pearsonr
import statsmodels
import statsmodels.stats.power

plt.rcParams["font.family"] = "serif"

mpl.rcParams.update({'font.size': 22})

top_dir = '/Volumes/USGS_Data/project/'


sitefile = top_dir + 'boxes/all_paths/t*plots/tstar_site.out'
#tuple with site name and tstar value
sitet = np.genfromtxt(sitefile, delimiter='\t', names = True, dtype = None)


sitetermfile = top_dir + 'event_boxes/site_terms.txt'
siteresid = np.genfromtxt(sitetermfile, delimiter='\t', names = True, dtype = None)

evboxes = ['PFO_3a', 'PMD_3a', 'SWS_3a', 'ERR_3a' ,'FRD_6', 'RDM_6', 'SWS_4', 'ERR_4']
path_residual_avg = []
pathint_avg = []
normpathint_avg = []
gradpathint_avg = []

#read in all evboxes files and compute the average path_residual, pathint, normpathint, gradpathint
for box in evboxes:
    infile =  top_dir + 'event_boxes/' + box + '_data.txt'
    data = np.genfromtxt(infile, delimiter='\t', names = True, dtype = None)
    path_residual = data['path_residual']
    pathint = data['pathint']
    normpathint = data['normpathint']
    gradpathint = data['gradpathint']
    
    path_residual_avg.append(np.mean(path_residual))
    pathint_avg.append(np.mean(pathint))
    normpathint_avg.append(np.mean(normpathint))
    gradpathint_avg.append(np.mean(gradpathint))
    print box, np.mean(path_residual)
    

#compare each path's average t* to path terms
mean_t_PFO_3a = 0.0031171819020574822
med_t_PFO_3a = 0.0043015897158499997
std_t_PFO_3a = 0.0088789492292223552

mean_t_PMD_3a = 0.0050330463401
med_t_PMD_3a = 0.00497124257326
std_t_PMD_3a = 0.00823585421476

mean_t_SWS_3a = 0.00556375667228
med_t_SWS_3a = 0.00664174677343
std_t_SWS_3a = 0.010412444393

mean_t_ERR_3a = 0.00291334405246
med_t_ERR_3a = 0.00321636143791
std_t_ERR_3a = 0.0115366611023

mean_t_FRD_6 = 0.0116251568606
med_t_FRD_6 = 0.0123565128715
std_t_FRD_6 = 0.00380773421202

mean_t_RDM_6 = 0.00331103325576
med_t_RDM_6 = 0.00323230300065
std_t_RDM_6 = 0.00162121395166

mean_t_SWS_4 = 0.00777983324141
med_t_SWS_4 = 0.00890941357616
std_t_SWS_4 = 0.0169305584537

mean_t_ERR_4 = 0.0191834133877
med_t_ERR_4 = 0.0194390018936
std_t_ERR_4 = 0.0162540751417

mean_t = [mean_t_PFO_3a, mean_t_PMD_3a, mean_t_SWS_3a, mean_t_ERR_3a, mean_t_FRD_6, mean_t_RDM_6, mean_t_SWS_4, mean_t_ERR_4]
x = mean_t
y = normpathint_avg


rval, pval = pearsonr(x, y)
power = statsmodels.stats.power.tt_solve_power(effect_size = rval, nobs = len(x), alpha = 0.05)

label1 = 'Pearson R : ' + str(round(rval,4))
label2 = 'power : ' + str(round(power,4))


fig = plt.figure(figsize = (20,15))
plt.scatter(x, y, s = 40)
plt.title('GMPE site residual vs. site t*')
plt.xlabel(r'$\mathrm{site \ t*}$')
plt.ylabel(r'$\mathrm{GMPE\ site \ residual}$')
plt.tick_params(axis='both', which='major', labelsize=18)
plt.tick_params(axis='both', which='both', length = 8, width = 2)

plt.annotate(label1, xy=(0.7, 0.85), xycoords='figure fraction', fontsize = 22)
plt.annotate(label2, xy=(0.7, 0.82), xycoords='figure fraction', fontsize = 22)

plt.show()
#plt.savefig(top_dir + 'boxes/all_paths/t*plots/pathterms/' + '')
