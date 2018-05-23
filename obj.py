#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 18 13:58:36 2018

@author: temp
"""

class Tobj(object):
    evid = ""
    tstar = 0
    dist = 0
    station = ""
    azimuth = 0
    evlat = 0
    evlon = 0
    evdep = 0

    # The class "constructor" - It's actually an initializer 
    def __init__(self, evid, evlat, evlon, evdep, tstar, dist, station, azimuth):
        self.evid = evid
        self.tstar = tstar
        self.dist = dist
        self.station = station
        self.azimuth = azimuth

def make_tobj(evid, evlat, evlon, evdep, tstar, dist, station, azimuth):
    tobj = Tobj(evid, evlat, evlon, evdep, tstar, dist, station, azimuth)
    return tobj


class source_obj(object):
    evid = ""
    m_cat = 0
    m_fit = 0
    lat = 0
    lon = 0
    depth = 0
    log_stressdrop = 0
    log_stressdrop_sig = 0 
    time = ''
    
    def __init__(self, evid, m_cat, m_fit, lat, lon, depth, log_stressdrop, log_stressdrop_sig, time):
        self.evid = evid
        self.m_cat = m_cat
        self.m_fit = m_fit
        self.lat = lat
        self.lon = lon
        self.depth = depth
        self.log_stressdrop = log_stressdrop
        self.log_stressdrop_sig = log_stressdrop_sig
        self.time = time

def make_source_obj(evid, m_cat, m_fit, lat, lon, depth, log_stressdrop, log_stressdrop_sig, time):
    s = source_obj(evid, m_cat, m_fit, lat, lon, depth, log_stressdrop, log_stressdrop_sig, time)
    return s