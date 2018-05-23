#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 16 13:53:33 2018

@author: temp
"""

import glob
import numpy as np
import obspy
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy.optimize import curve_fit
from uncertainties import ufloat
import uncertainties.umath as umath
from matplotlib.gridspec import GridSpec
plt.style.use("classic")
from obj import source_obj
import cPickle as pickle
from datetime import datetime

#plt.rcParams["font.family"] = "serif"
plt.rcParams['axes.axisbelow'] = True
mpl.rcParams.update({'font.size': 18})

#change these magnitude upper and lower bounds
#these are for Joe's ~3 eventds

beta = 3500. #3500m/s
U = 0.63#0.63
rho = 2750. #2750 kg/m^3

#read in all event spectra files
top_dir = '/Volumes/USGS_Data/project'
box = 'all_paths'
boxpath = top_dir + '/boxes/' + box
#event_spectra_dir = boxpath + '/secondo_rebin3_constrained_2013_11_30_11_36_35/'
#event_spectra_dir = boxpath + '/secondo_rebin3_constrained_2010_05_25_19_49_51/'
event_spectra_dir = boxpath + '/secondo_rebin3_constrained_2013_11_30_11_36_35/'


event_spectra = glob.glob(event_spectra_dir + '[2]*.out')


#find events in catalog that are in mag range
catalog = top_dir + '/catalogs/all_paths_M2.5_USGS_Catalog.txt'
cat = np.genfromtxt(catalog, comments = '#', delimiter = '|', dtype = None, usecols = [1,10])
event = []
magl = []
for i in range(len(cat)):
    m = cat[i][1]
    magl.append(cat[i][1])
    time = obspy.core.utcdatetime.UTCDateTime(cat[i][0].split('.')[0])
    ev = str(time.year).zfill(4) + '_' + str(time.month).zfill(2) + '_' + str(time.day).zfill(2) + '_' + str(time.hour).zfill(2) + '_' + str(time.minute).zfill(2) + '_' + str(time.second).zfill(2)
    event.append(ev)


#define Brune velocity spectra
def Brune(f, fc, omega0):
     return (2.*np.pi*f*omega0)/(1+(f/fc)**2.)
#def Brune(f, M0, rho, beta, stressdrop, omega0):
#     return (2.*np.pi*f*M0*0.63/(4.*np.pi*rho*beta**3))/(1+(f/(beta*(stressdrop/(8.47*M0))**(1./3)))**2.)

s = []
s_w = []
s_sig = []
Moment = []
Moment_sig = []
ml_calc = []
ml_calc_sig = []
omega0_calc = []
omega0_calc_sig =[]
fc_calc = []
fc_calc_sig = []
evid = []
magl_cat = []
error = 0

#for each event, use the curve fit to find the moment and stress drop from corner frequency and omega0
for i in range(len(event)):
#for i in range(20,40):
#    print event[i]
    f = event_spectra_dir + event[i] + '.out'
#    f = event_spectra[i]
    data = np.genfromtxt(f, dtype = float, comments = '#', delimiter = None, usecols = (0,1))#only read in first two cols
    freq = data[:,0][20:70]
    spec = data[:,1][20:70]
    
    try:
        popt, pcov = curve_fit(Brune, freq, spec, bounds=((0,0), (100., 1000.)))
    except RuntimeError:
        print "Cannot converge on an answer for event: " + event[i]
        error +=1
        continue
        
    fc_fit, omega0_fit = popt
    fc_std, omega0_std = np.sqrt(np.diag(pcov))
    
    fc = ufloat(fc_fit, fc_std)
    omega0 = ufloat(omega0_fit, omega0_std)
    
    M0 = (omega0*4.*np.pi*rho*beta**3.)/U
    stressdrop = ((fc/beta)**3.)*(8.47*M0)
    stressdrop_log = umath.log10(stressdrop/1e6)

#    print stressdrop/1e6
    
#    stressdropup = (((fc_fit+fc_std)/beta)**3.)*(8.47*M0)
#    stressdroplb = 
    
    M = (2./3.)*(umath.log10(M0) - 9.05)
    if M > 3:
        ml = M
    else:
        ml = (M - 0.884)/.667
    

    s.append(stressdrop_log.nominal_value)
    s_w.append(1./stressdrop_log.std_dev**2)
    if stressdrop_log.std_dev > 0.001:
        s_sig.append(stressdrop_log.std_dev)
    else:
        s_sig.append(0.001)
    Moment.append(M0.nominal_value)
    Moment_sig.append(M0.std_dev)
    ml_calc.append(ml.nominal_value)
    ml_calc_sig.append(ml.std_dev)
    omega0_calc.append(omega0.nominal_value)
    omega0_calc_sig.append(omega0.std_dev)
    fc_calc.append(fc.nominal_value)
    fc_calc_sig.append(fc.std_dev)
    evid.append(event[i])
    magl_cat.append(magl[i])
    
#    fig = plt.figure(figsize = (16,16))
#    plt.ylabel('Velocity amplitude (m)', fontsize = 16)
#    plt.xlim(0.5,70)
#    plt.ylim(.01,8)
#    plt.loglog(data[:,0] , data[:,1], color = 'green', label = 'event spectra')
#    plt.grid()
#    plt.loglog(freq, Brune(freq, fc_fit, omega0_fit), color = 'blue', label = 'Brune spectra')
#    plt.legend(loc = 'upper left', fontsize = 16)
#    plt.xlabel('Frequency (Hz)', fontsize = 16)
#    plt.title(event_spectra[i].split('/')[7].split('.')[0])
#    plt.tick_params(axis='both', which='major', labelsize=15)
#    plt.tick_params(axis='both', which='both', length = 5, width = 1)
#    plt.text(1, .035, 'ml cat: ' + "{:.3f}".format(magl[i]) + '   ml calc: ' + "{:.3f}".format(ml), fontsize = 16)
#    plt.text(1, .02, 'M0: ' + "{:.3e}".format(M0) + '   stress drop: ' + "{:.3f}".format(stressdrop/(10.**6.)) + ' MPa', fontsize = 16)
#    plt.text(1, .012, 'fc: ' + "{:.3f}".format(fc) + '   omega0: ' + "{:.3f}".format(omega0), fontsize = 16)
#    plt.show()

#    print event[i]
#    print 'ml cat: ' + "{:.3f}".format(magl[i]) + ' ml calc: ' + "{:.3f}".format(ml)
#    print 'fc fit: ' + "{:.3f}".format(fc_fit) + ' fc std: '+ "{:.3f}".format(fc_std) 
#    print 'omega0 fit: ' + "{:.3f}".format(omega0) #+ ' omega0 std: '+ "{:.3f}".format(omega0_std) 
#    print ' stress drop: ' + "{:.3f}".format(stressdrop/1e6) + ' MPa'

outfile = open(top_dir + '/source_params/secondo_2013_source.out', 'w')
#out = (np.array([np.array(evid), np.array(magl_cat), np.array(ml_calc), np.array(fc_calc), np.array(s), np.array(omega0), np.array(Moment)], dtype=['|S32', np.float64, np.float64, np.float64, np.float64, np.float64, np.float64]).T)
out = np.array([np.array(evid), np.array(magl_cat), np.array(ml_calc), np.array(ml_calc_sig), np.array(fc_calc),np.array(fc_calc_sig), np.array(s), np.array(s_sig), np.array(omega0_calc),np.array(omega0_calc_sig), np.array(Moment), np.array(Moment_sig)], dtype = object).T

outfile.write('#evid \t \t cat_ml \t fit_ml \t fit_ml_sig \t fit_fc \t fit_fc_sig \t log_stressdrop(MPa)_fit \t log_stressdrop_sig \t omega0_fit \t omega0_sig \t M0_fit \t M0_sig \n')
np.savetxt(outfile, out, fmt=['%s', '%.3f', '%.3f', '%.5f', '%.3f', '%.3f', '%.3f',' %.3f', '%.3f', '%.3f','%.3f', '%.3f'], delimiter='\t')
#np.savetxt(outfile, out, fmt= '%s', delimiter='\t')
#
outfile.close()


##############################
f = top_dir + '/source_params/secondo_2013_source.out'
data = np.genfromtxt(f, comments = '#', delimiter = None, dtype = None, names=True, encoding = None)
s = data['log_stressdropMPa_fit']

s_w = (1./data['log_stressdrop_sig']**2)

ml = data['fit_ml']
ml_w = (1./data['fit_ml_sig']**2)
cat_ml = data['cat_ml']

#weighted average
avg_w = np.sum([s_w[i]*s[i] for i in range(len(s))])/np.sum(s_w)
sig_avgw = 1/np.sqrt(np.sum(s_w))

#calculate ml and compare to actual ml
print 'stress drop mean, median, sigma: ', np.mean(s), np.median(s), np.std(s)
print 'weighted average :' + str(avg_w) + ' +- ' + str(sig_avgw)

label = 'mean : ' + "{:.3f}".format(10.**np.mean(s)) + ' */' + u"\u00F7" + "{:.3f}".format(np.std(s)) + '\n' + 'median : ' + "{:.3f}".format(10**np.median(s)) + '\n \n' + 'weighted mean : ' + "{:.3f}".format(10.**avg_w) + ' */' + u"\u00F7" + "{:.3f}".format(sig_avgw)

plt.figure(figsize = (14,14))
plt.tick_params(axis='both', which='major', labelsize=18)
plt.tick_params(axis='both', which='both', length = 5, width = 1)
plt.ylabel('events', fontsize = 28)
plt.xlabel('log stress drop (MPa)', fontsize = 28)
#plt.xlim(0,50)
plt.annotate(label, xy = (0.03,0.85), xycoords = 'axes fraction', fontsize = 24)
plt.hist(s, bins = 50)
plt.show()
plt.savefig(top_dir + '/source_params/stress_drop_histo.png')


gs = GridSpec(2,1)
fig = plt.figure(figsize = (14,16))
#fig.text(0.01, 0.5,'counts', va='center', rotation='vertical')
plt.subplot(gs.new_subplotspec((0, 0), colspan=1))
#weighted average
avg_w = np.sum([ml_w[i]*ml[i] for i in range(len(ml))])/np.sum(ml_w)
sig_avgw = 1/np.sqrt(np.sum(ml_w))

label = 'mean : ' + "{:.3f}".format(np.mean(ml)) + '\n' + 'median : ' + "{:.3f}".format(np.median(ml)) + '\n' + 'std dev : ' + "{:.3f}".format(np.std(ml)) + '\n \n' + 'weighted mean : ' + "{:.3f}".format(avg_w) + '\n' + 'weighted std dev : ' + "{:.3f}".format(sig_avgw)

plt.tick_params(axis='both', which='major', labelsize=18)
plt.tick_params(axis='both', which='both', length = 5, width = 1)
plt.ylabel('events', fontsize = 26)
plt.xlabel('magnitude', fontsize = 26)
plt.title('calculated magnitude', fontsize = 28)
plt.annotate(label, xy = (0.5,0.5), xycoords = 'axes fraction', fontsize = 24)
plt.hist(ml, bins = np.arange(2.5,6,0.25), color = 'green')

label = 'mean : ' + "{:.3f}".format(np.mean(cat_ml)) + '\n' + 'median : ' + "{:.3f}".format(np.median(cat_ml)) + '\n' + 'std dev : ' + "{:.3f}".format(np.std(cat_ml))

plt.subplot(gs.new_subplotspec((1, 0), colspan=1))
plt.tick_params(axis='both', which='major', labelsize=18)
plt.tick_params(axis='both', which='both', length = 5, width = 1)
plt.ylabel('events', fontsize = 26)
plt.xlabel('magnitude', fontsize = 26)
plt.title('catalog magnitude', fontsize = 28)
plt.annotate(label, xy = (0.5,0.5), xycoords = 'axes fraction', fontsize = 24)

plt.hist(cat_ml, bins = np.arange(2.5,6,0.25))

plt.show()
plt.savefig(top_dir + '/source_params/m_histo.png')

#save as a pickle object
f = top_dir + '/source_params/secondo_2013_source.out'
data = np.genfromtxt(f, comments = '#', delimiter = None, dtype = None, names=True, encoding = None)
s = data['log_stressdropMPa_fit']
s_sig = data['log_stressdrop_sig']
s_w = (1./data['log_stressdrop_sig']**2)
ml = data['fit_ml']
#ml_w = (1./data['fit_ml_sig']**2)
cat_ml = data['cat_ml']
evid = data['evid']
#weighted average
avg_w = np.sum([s_w[i]*s[i] for i in range(len(s))])/np.sum(s_w)
sig_avgw = 1/np.sqrt(np.sum(s_w))

#find events in catalog that are in mag range
catalog = top_dir + '/catalogs/all_paths_M2.5_USGS_Catalog.txt'
cat = np.genfromtxt(catalog, comments = '#', delimiter = '|', dtype = None, encoding = None, names = True)
event_cat = []
magl_cat = []
dep_cat = []
lat_cat = []
lon_cat = []
time_cat = []

for i in range(len(cat)):
    m = cat[i]['Magnitude']
    magl_cat.append(m)
    time = obspy.core.utcdatetime.UTCDateTime(cat[i]['Time'].split('.')[0])
    time_cat.append(time)
    ev = str(time.year).zfill(4) + '_' + str(time.month).zfill(2) + '_' + str(time.day).zfill(2) + '_' + str(time.hour).zfill(2) + '_' + str(time.minute).zfill(2) + '_' + str(time.second).zfill(2)
    event_cat.append(ev)
#    time_cat.append(datetime(time.year, time.month, time.day, time.hour, time.minute, time.second))
    dep_cat.append(cat[i]['Depthkm'])
    lat_cat.append(cat[i]['Latitude'])
    lon_cat.append(cat[i]['Longitude'])

x = 0
for i in range(len(event_cat)):
    for j in range(len(evid)):
        if event_cat[i] == evid[j]:
            if s_sig[j] < 0.5:
                x +=1 
ind = 0
#2d array
a = np.zeros(x, dtype={'names':('evid', 'm_cat', 'm_fit', 'lat', 'lon', 'depth', 'log_stressdrop', 'log_stressdrop_sig', 'time'), 'formats':('|S32', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', '|S32')})
#a = np.zeros(x, dtype={'names':('evid', 'm_cat', 'm_fit', 'lat', 'lon', 'depth', 'log_stressdrop', 'log_stressdrop_sig', 'time'), 'formats':('|S32', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'datetime64')})

#match catalog
for i in range(len(event_cat)):
    for j in range(len(evid)):
        if event_cat[i] == evid[j]:
            if s_sig[j] < 0.5:
                a[ind][0] = event_cat[i]
                a[ind][1] = magl_cat[i]
                a[ind][2] = ml[j]
                a[ind][3] = lat_cat[i]
                a[ind][4] = lon_cat[i]
                a[ind][5] = dep_cat[i]
                a[ind][6] = s[j]
                a[ind][7] = s_sig[j]
                a[ind][8] = time_cat[i]
                ind += 1

#reverse the array
a = a[::-1]

def make_source_obj(evid, m_cat, m_fit, lat, lon, depth, log_stressdrop, log_stressdrop_sig, time):
    s = source_obj(evid, m_cat, m_fit, lat, lon, depth, log_stressdrop, log_stressdrop_sig, time)
    return s

x = make_source_obj(a['evid'], a['m_cat'], a['m_fit'], a['lat'], a['lon'], a['depth'], a['log_stressdrop'], a['log_stressdrop_sig'], a['time'])
output = open(top_dir + '/source_params/Brune_fit_params_sig_lessthan_0.5.pckl', 'wb')
pickle.dump(x, output)
output.close()