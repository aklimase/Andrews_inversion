#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 12 15:25:08 2018

@author: temp

reads in the object with Valerie's path residuals and the object with tstar values
station by station, finds the t* paths with events closest to residual paths within 5 km
plots the t* vs path residual for each t* path
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
#import obspy
#import glob
import dread
from scipy.stats import pearsonr
import statsmodels.stats.power
#import matplotlib.patches as mpatches
import cPickle as pickle
from temp_obj import Tobj

####set the maxiumum distance for finding matches here, km
match_dist = 20

plt.rcParams["font.family"] = "serif"
mpl.rcParams.update({'font.size': 22})

topdir = '/Volumes/USGS_Data/project/'

# Define the path to open:
rpath = topdir + 'event_boxes/mixedcoeff_v3anza2013_pga_noVs30_5coeff_a4_-1.2_pga__ncoeff5_Mc_8.5_VR_99.6_a4_-1.2_robj.pckl'

# Open the file to read ('r'), save to rfile:
rfile = open(rpath,'r')

# Load the object inside of rfile into a variable (this could take a while, the object is big):
robj = pickle.load(rfile)

# Close the file:
rfile.close()

tpath = topdir + 'boxes/all_paths/path_tstar_evlocations.pckl'

tfile = open(tpath, 'r')

tobj = pickle.load(tfile)
tfile.close()

sta_set = list(set(tobj.station))



outfile = open(topdir + 'boxes/all_paths/t*plots/pathterms/match_pathresid_' + str(match_dist) + 'km.txt', 'w')


#print '{0:23s} {1:10s} {2:12s} {3:12s} {4:13s} {5:12s} {6:12s} {7:12s} {8:12s} {9:12s} {10:12s}'.format('Event id cat', 'station', 'Lon cat', 'Lat cat', 'Depth cat', 'Lon evbox', 'Lat evbox', 'Depth evbox', 'distance', 'tstar', 'path residual')

#here are the matches
t = []
resid = []
s = []

#in order to use np loadtxt save to arrays to transform
v1 =[]
v2 = []
v3 = []
v4 = []
v5 = []
v6 = []
v7 = []
v8 = []


#for each station in both data sets
for k in range(len(sta_set)):
    print sta_set[k]
    t_index = [i for i, j in enumerate(tobj.station) if j == sta_set[k]]
    r_index = [i for i, j in enumerate(robj.sta) if j == sta_set[k]]
    print(len(t_index), len(r_index))
    
    #variables for stations of interest 
    t_evid = [tobj.evid[i] for i in t_index]
    t_evlat = [tobj.evlat[i] for i in t_index]
    t_evlon = [tobj.evlon[i] for i in t_index]
    t_evdep = [tobj.evdep[i] for i in t_index]
    t_tstar = [tobj.tstar[i] for i in t_index]
    
    r_evlat = [robj.elat[i] for i in r_index]
    r_evlon = [robj.elon[i] for i in r_index]
    r_evdep = [robj.edepth[i] for i in r_index]
    r_residual = [robj.path_terms[i] for i in r_index]
    
    #for each t* event, find the closest residual event
    for l in range(len(t_evlat)):
#    for l in range(10):
        d = []
        for m in range(len(r_evlat)):
            dist = dread.compute_rrup(t_evlon[l], t_evlat[l], t_evdep[l], r_evlon[m], r_evlat[m], -1*r_evdep[m]) #in km
#            print(t_evlon[l], t_evlat[l], t_evdep[l], r_evlon[m], r_evlat[m], -1*r_evdep[m])
            d.append(dist)
        if min(d) < match_dist:
            #it's a match!!
            ind = d.index(min(d))
#            print '{0:23s} {1:6s} {2:15f} {3:12f} {4:12f} {5:12f} {6:12f} {7:12f} {8:12f} {9:12f} {10:12f}'.format(t_evid[l], sta_set[k], t_evlon[l], t_evlat[l], t_evdep[l], r_evlon[ind], r_evlat[ind], r_evdep[ind], d[ind], t_tstar[l], r_residual[ind])
            v1.append(t_evid[l])
            s.append(sta_set[k])
            v2.append(t_evlon[l])
            v3.append(t_evlat[l])
            v4.append(t_evdep[l])
            v5.append(r_evlon[ind])
            v6.append(r_evlat[ind])
            v7.append(r_evdep[ind])
            v8.append(d[ind])
            t.append(t_tstar[l])
            resid.append(r_residual[ind])
            
#            outfile.write(str(t_evid[l]) + '\t'+ str(sta_set[k])+'\t'+ str(t_evlon[l])+'\t'+ str(t_evlat[l])+'\t'+ str(t_evdep[l])+'\t'+ str(r_evlon[ind])+'\t'+ str(r_evlat[ind])+'\t'+ str(r_evdep[ind])+'\t'+ str(round(d[ind],4))+'\t'+ str(round(t_tstar[l],4))+'\t'+ str(round(r_residual[ind],4)) + '\n')

out = (np.array([v1, s, v2, v3, v4, v5, v6, v7, v8, t, resid]).T)
outfile.write('#Event_id_cat \t station \t Lon_cat \t Lat_cat \t Depth_cat \t Lon_evbox \t Lat_evbox \t Depth_evbox \t distance \t tstar \t path_residual \n')
np.savetxt(outfile, out, fmt='%s', delimiter='\t')
outfile.close()
            
fig = plt.subplots(figsize = (15,15))
x = t 
y = resid
    
#plt.title(r'$\mathrm{residual\ t*\ vs.\ path\ distance}$')
plt.xlabel(r'$\mathrm{path \ t*}$')
plt.ylabel(r'$\mathrm{GMPE \ path \ residual}$')
plt.tick_params(axis='both', which='major', labelsize=18)
plt.tick_params(axis='both', which='both', length = 8, width = 2)

#colors = {'TPFO_3a':'red','PFO_3a':'red','PMD_3a':'red','SWS_3a':'green','ERR_3a':'green', 'FRD_6':'blue','RDM_6':'blue','SWS_4':'purple','ERR_4':'purple', 'TPFO_5':'orange','PFO_5':'orange','PMD_5':'orange', 'TRO_3a':'yellow', 'TRO_5':'yellow'}

#scatter = plt.scatter(x,y, color=[ colors[i] for i in box ], label = [colors.keys()[i] for i in range(len(colors.keys()))])

scatter = plt.scatter(x,y)

rval, pval = pearsonr(x, y)
power = statsmodels.stats.power.tt_solve_power(effect_size = rval, nobs = len(x), alpha = 0.05)
label1 = 'Pearsons R '  + ': ' + str(round(rval,4))
label2 = 'power ' +  ': ' + str(round(power,4))
plt.annotate(label1, xy=(0.15, 0.17), xycoords='figure fraction', fontsize = 18)
plt.annotate(label2, xy=(0.15, 0.13), xycoords='figure fraction', fontsize = 18)


#red_patch = mpatches.Patch(color='red', label='Imperial Valley PFO PMD')
#green_patch = mpatches.Patch(color='green', label='Imperial Valley SWS ERR')
#blue_patch = mpatches.Patch(color='blue', label='Riverside FRD RDM')
#purple_patch = mpatches.Patch(color='purple', label='Salton Trough SWS ERR')
#orange_patch = mpatches.Patch(color='orange', label='5 PFO PMD')
#yellow_patch = mpatches.Patch(color='yellow', label='5, 3a TRO')
#plt.legend(handles = [red_patch, green_patch, blue_patch, purple_patch, orange_patch , yellow_patch], fontsize = 15)
#
plt.title('GMPE path residual vs. path t*', fontsize = 22)

plt.savefig(topdir + 'boxes/all_paths/t*plots/pathterms/match_pathresid_' + str(match_dist) + 'km.png')

plt.show()
























