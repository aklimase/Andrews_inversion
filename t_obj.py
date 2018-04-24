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
