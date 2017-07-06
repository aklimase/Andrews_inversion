#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 19 14:44:20 2017

@author: escuser
"""

#plotting corrected data
from obspy import read

eventid = '15001500'
network = 'AZ'
station = 'KNW'
full_channel = 'HHZ'


boxpath = '/Users/escuser/Documents/Alexis_Data/cut_sac_files_2'
icorr_path = boxpath + '/icorrdata'
sacfile = boxpath + '/rawdata/' + eventid + '.' + network + '.' + station + '.' + full_channel + '.sac'
icorr_sacfile = icorr_path + '/' + eventid + '.' + network + '.' + station + '.' + full_channel + '.sac'


stream1 = read(sacfile)
tr1 = stream1[0]
dt = tr1.stats.starttime
print 'Magnitude: ', tr1.stats.sac.get('mag')
print 'idep: ', tr1.stats.sac.get('idep')
#mag = tr1.stats.starttime
stream1.plot(size = (800, 400), starttime=dt + 30, endtime=dt + 70)


stream2 = read(icorr_sacfile)
tr2 = stream2[0]
dt = tr2.stats.starttime
print 'idep: ', tr2.stats.sac.get('idep')
#mag = tr2.stats.starttime
stream2.plot(size = (800, 400), starttime=dt + 30, endtime=dt + 70)

stream1 = read(sacfile + '0.1')
tr1 = stream1[0]
dt = tr1.stats.starttime
print 'idep: ', tr2.stats.sac.get('idep')
#mag = tr1.stats.starttime
stream1.plot(size = (800, 400), starttime=dt + 30, endtime=dt + 70)