#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu May  3 13:25:56 2018

@author: temp
"""

import glob
import numpy as np
import obspy
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.gridspec import GridSpec
from datetime import datetime
from datetime import timedelta
import dread
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.io.img_tiles as cimgt
import matplotlib.ticker as mtick
import cPickle as pickle
from shapely.geometry.polygon import Polygon
import matplotlib.path as mpltPath
import dread as dr
from obj import source_obj

plt.style.use("classic")
plt.rcParams['axes.axisbelow'] = True
mpl.rcParams.update({'font.size': 18})
mpl.rcParams['figure.subplot.left'] = 0.08
mpl.rcParams['figure.subplot.right'] = 0.95


#read in all event spectra files
top_dir = '/Volumes/USGS_Data/project'
box = 'all_paths'
boxpath = top_dir + '/boxes/' + box
event_spectra_dir = boxpath + '/secondo_rebin3_constrained_2013_11_30_11_36_35/'
event_spectra = glob.glob(event_spectra_dir + '[2]*.out')

# residuals
rpath = top_dir + '/event_boxes/mixedcoeff_v3anza2013_pga_noVs30_5coeff_a4_-1.2_pga__ncoeff5_Mc_8.5_VR_99.6_a4_-1.2_robj.pckl'
rfile = open(rpath,'r')
robj = pickle.load(rfile)
rfile.close()

r_evlat = robj.elat
r_evlon = robj.elon
r_evdep = robj.edepth
r_Eresidual = robj.E_residual

faultfile= top_dir + '/catalogs/Holocene_LatestPleistocene_118.0w_115.0w_32.3n_34.4n.pckl'
fault_segments=dr.read_obj_list(faultfile)

elMayortime = obspy.core.utcdatetime.UTCDateTime('2010-04-04T22:40:42')
elMayorlon = -115.295
elMayorlat = 32.286
elMayordepth = 10.0

sfile = open(top_dir + '/source_params/Brune_fit_params_sig_lessthan_0.5.pckl','r')
sobj = pickle.load(sfile)
sfile.close()

vertices_LagunaSalada = [(-116.1, 32.75), (-115.7, 32.75), (-115.3, 32.4), (-115.25, 32.3), (-115.6, 32.3), (-116.1, 32.75)]
vertices_CerroPrieto = [(-115.42, 32.6), (-115.15, 32.6), (-115.0, 32.3), (-115.25, 32.3), (-115.42, 32.6)]
vertices_ImperialFault = [(-115.85, 33.4), (-115.4,33.4),(-115.15,32.6),(-115.5,32.6),(-115.85, 33.4)]
vertices_Elsinor = [(-117.2, 33.5), (-116.95, 33.6), (-115.85, 32.75), (-116.3,32.75), (-117.2, 33.5)]
vertices_SanJacinto = [(-117.45, 34.1), (-117.15, 34.2), (-115.7, 33.05), (-115.95, 32.85), (-117.45, 34.07)]
vertices_SanAndreas = [(-117.35,34.4),(-116.2,34.4),(-116.2,33.8),(-116.7, 33.8), (-117.35,34.4)]
vertices_Riverside = [(-118, 34.2),(-117.5, 34.2),(-117.5, 33.75),(-118, 33.75)]
vertices_other1 = [(-118, 33.3),(-117.3,33.3),(-116.8, 32.2),(-118, 32.3), (-118, 33.3)]
#make regional object
vertex_list = [vertices_LagunaSalada, vertices_CerroPrieto, vertices_ImperialFault, vertices_Elsinor, vertices_SanJacinto, vertices_SanAndreas, vertices_Riverside, vertices_other1]
color = ['red', 'blue', 'orange', 'green', 'magenta', 'cyan', 'yellow', 'olive']
labels = ['Laguna Salada', 'Cerro Prieto', 'Imperial Fault', 'Elsinor', 'San Jacinto', 'San Andreas', 'Riverside', 'other 1']

def magnitude_map(mag_fit, mag_cat, lat, lon, save):
    '''
    Maps magnitudes for all event on map
    Input:
        mag_fit:    list/array of event magntiudes
        mag_cat:    list/array of event magntiudes
        lat:        list/array of event latitudes
        lon:        list/array of event longitudes
        save:       yes or no to save figure
    Output:
        plots a map
    '''
    
    import matplotlib.pyplot as plt
    import matplotlib as mpl
    from matplotlib.gridspec import GridSpec
    import cartopy.crs as ccrs
    import cartopy.io.img_tiles as cimgt
    
    gs = GridSpec(1,2)
    # make a color map
    cmap = mpl.cm.get_cmap('viridis')
    normalize = mpl.colors.Normalize(vmin=min(mag_fit), vmax=max(mag_fit))
    colors = [cmap(normalize(value)) for value in mag_fit]
    s_m = mpl.cm.ScalarMappable(cmap = cmap, norm=normalize)
    s_m.set_array([])
    
    fig = plt.figure(figsize = (20,10))
    
    plt.subplot(gs.new_subplotspec((0, 0), colspan=1))
    plt.title('Calculated magnitude', fontsize = 30)
    plt.tick_params(axis='both', which='major', labelsize=18)
    plt.tick_params(axis='both', which='both', length = 5, width = 1)
    plt.ylabel('latitude', fontsize = 28)
    plt.xlabel('longitude', fontsize = 28)
    plt.scatter(lon, lat, color = colors, s = 30)
    plt.xlim(-118.1, -115)
    plt.ylim(32.2,34.5)
    
    colors = [cmap(normalize(value)) for value in mag_cat]
    plt.subplot(gs.new_subplotspec((0, 1), colspan=1))
    plt.title('Catalog magnitude', fontsize = 30)
    plt.tick_params(axis='both', which='major', labelsize=18)
    plt.tick_params(axis='both', which='both', length = 5, width = 1)
    plt.ylabel('latitude', fontsize = 28)
    plt.xlabel('longitude', fontsize = 28)
    plt.scatter(lon, lat, color = colors, s = 30)
    #ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%.1f'))
    plt.xlim(-118.1, -115)
    plt.ylim(32.2,34.5)
    
    
    fig.subplots_adjust(right=0.9)
    cbar_ax = fig.add_axes([0.92, 0.15, 0.02, 0.75])
    cbar = plt.colorbar(s_m, cax=cbar_ax)
    cbar.set_label(ur"magnitude", fontsize = 22)#"$azimuth$ (\u00b0)"
    cbar.ax.tick_params(labelsize = 18)
    
    plt.show()
    if save == True:
        plt.savefig(top_dir + '/source_params/mag_map.png')

#magnitude_map(mag_fit = sobj.m_fit, mag_cat = sobj.m_cat, lat = sobj.lat, lon = sobj.lon, save = False)

def stressdrop_map(s, lat, lon, save, regions):
    '''
    Maps stress drops for all event on map
    Input:
        s:          list/array of log stressdrops
        lat:        list/array of event latitudes
        lon:        list/array of event longitudes
        save:       yes or no to save figure

    Output:
        plots a map
    '''

    import matplotlib.pyplot as plt
    import matplotlib as mpl
    import cartopy.crs as ccrs
    import cartopy.io.img_tiles as cimgt

    cmap = mpl.cm.get_cmap('gist_rainbow')
    normalize = mpl.colors.Normalize(vmin=min(s), vmax=max(s))
    colors = [cmap(normalize(value)) for value in s]
    s_m = mpl.cm.ScalarMappable(cmap = cmap, norm=normalize)
    s_m.set_array([])

    stamen_terrain = cimgt.StamenTerrain(desired_tile_form="L")
    fig = plt.figure(figsize = (18,16))
    ax = plt.axes(projection=stamen_terrain.crs)
    ax.set_extent([-118.01, -114.99, 32.29, 34.41])

    ax.add_image(stamen_terrain, 10, alpha = 0.5, cmap = 'gray')
    plt.tick_params(axis='both', which='major', labelsize=18)
    plt.tick_params(axis='both', which='both', length = 5, width = 1)

    xticks = [-118, -117, -116, -115]
    yticks = [32.3, 33, 33.5, 34, 34.4]
    ax.gridlines(xlocs=xticks, ylocs=yticks, draw_labels = True)

    for segment_i in range(len(fault_segments)):
        fault=fault_segments[segment_i]
        fault_z=np.zeros(len(fault))
        ax.plot(fault[:,0],fault[:,1],fault_z,color='k', transform=ccrs.Geodetic())

    plt.scatter(lon, lat, marker='o', color= colors, s=30, alpha=1.0, transform=ccrs.Geodetic())
#    geodetic_transform = ccrs.Geodetic()._as_mpl_transform(ax)

    #if regions is True, add polygons
    if regions == True:
        for i in range(len(vertex_list)):
            poly = Polygon(vertex_list[i])
            x,y = poly.exterior.xy
            plt.plot(x,y, color = color[i], transform=ccrs.Geodetic(), label = labels[i], lw = 2)

    plt.legend(loc=(0.1,-0.15),ncol=4)

    fig.subplots_adjust(right=0.88)
    cbar_ax = fig.add_axes([0.92, 0.15, 0.02, 0.75])
    cbar = plt.colorbar(s_m, cax=cbar_ax)
    cbar.set_label(ur"stressdrop", fontsize = 22)#"$azimuth$ (\u00b0)"
    cbar.ax.tick_params(labelsize = 18)

    plt.show()
    if save == True:
        plt.savefig(top_dir + '/source_params/stressdrop_map_regions.png')

#stressdrop_map(s = sobj.log_stressdrop, lat = sobj.lat, lon = sobj.lon, save = False, regions = True)

def stressdrop_vs_depth_mag(s, d, mag, sig, name, save):
    '''
    Maps stress drops vs. depths/mag for given events
    Input:
        s:          list/array of log stressdrops
        d:          list/array of event depths
        mag:        list/array of event magnitudes (catalog)
        name:       string name of polygon
        save:       True or False to save
    Output:
        plots 2 scatter plots of stress drop vs mag and stress drop versus depth
    '''
    import matplotlib.pyplot as plt
    import matplotlib as mpl    

    fig = plt.figure(figsize = (18,10))
    plt.subplot(111)
    
        #plt.scatter(a['time'], a['log_stressdrop'])
    plt.tick_params(axis='both', which='major', labelsize=18)
    plt.tick_params(axis='both', which='both', length = 5, width = 1)
    cmap = mpl.cm.get_cmap('viridis')
    normalize = mpl.colors.Normalize(vmin=min(mag), vmax=max(mag))
    colors = [cmap(normalize(value)) for value in mag]
    s_m = mpl.cm.ScalarMappable(cmap = cmap, norm=normalize)
    s_m.set_array([])
    
    for i in range(len(s)):
            plt.errorbar(x=s[i], y=d[i], xerr = sig[i], marker = '.', c = colors[i])
    plt.scatter(x=s, y=d, c = colors, edgecolors = colors)
    
    plt.ylabel('depth(km)', fontsize = 28)
    plt.xlabel('log(stress drop)', fontsize = 28)
    
    plt.gca().invert_yaxis()
    
    fig.subplots_adjust(right=0.9)
    cbar_ax = fig.add_axes([0.92, 0.15, 0.02, 0.75])
    cbar = plt.colorbar(s_m, cax=cbar_ax)
    cbar.set_label(ur"magnitude", fontsize = 22)#"$azimuth$ (\u00b0)"
    cbar.ax.tick_params(labelsize = 18)
    plt.show()
    plt.savefig(top_dir + '/source_params/polygon/stressdrop_vs_depth_' + name + '.png')
    
    fig = plt.figure(figsize = (18,10))
    plt.subplot(111)
        #plt.scatter(a['time'], a['log_stressdrop'])
    plt.tick_params(axis='both', which='major', labelsize=18)
    plt.tick_params(axis='both', which='both', length = 5, width = 1)
    cmap = mpl.cm.get_cmap('viridis')
    normalize = mpl.colors.Normalize(vmin=min(d), vmax=max(d))
    colors = [cmap(normalize(value)) for value in d]
    s_m = mpl.cm.ScalarMappable(cmap = cmap, norm=normalize)
    s_m.set_array([])
    
    for i in range(len(s)):
            plt.errorbar(x=s[i], y=mag[i], xerr = sig[i], marker = '.', c = colors[i])
    plt.scatter(x=s, y=mag, c = colors, edgecolors = colors)
    
    plt.ylabel('magnitude', fontsize = 28)
    plt.xlabel('log(stress drop)', fontsize = 28)

    fig.subplots_adjust(right=0.9)
    cbar_ax = fig.add_axes([0.92, 0.15, 0.02, 0.75])
    cbar = plt.colorbar(s_m, cax=cbar_ax)
    cbar.set_label(ur"depth(km)", fontsize = 22)#"$azimuth$ (\u00b0)"
    cbar.ax.tick_params(labelsize = 18)
    plt.show()
    
    if save == True:
        plt.savefig(top_dir + '/source_params/polygon/stressdrop_vs_mag_' + name + '.png')

#stressdrop_vs_depth_mag(s = sobj.log_stressdrop, d = sobj.depth, mag = sobj.m_cat, sig = sobj.log_stressdrop_sig, name = 'none', save = False)

def stressdrop_vs_time(ev, lon, lat, s, sig, m, name, x0, x1, save, days, regions):
    '''
    Maps stress drops for all event on map
    Input:
        t:          list/array of log stressdrops
        s:          list/array of event depths
        sig:        list/array of event magnitudes (catalog)
        m:          string name of polygon
        name:       name to saving plot
        x0:         utc datetime, first in range for plotting
        x1:         utc datetime, first in range for plotting
        save:       True or False to save
    Output:
        plots 2 scatter plots of stress drop vs mag and stress drop versus depth
    '''
    import matplotlib.pyplot as plt
    import matplotlib as mpl


    time = []
    #make evid into datetimes
    for i in range(len(ev)):
        x = (ev[i]).split('_')
        time.append(datetime(int(x[0]), int(x[1]), int(x[2]), int(x[3]), int(x[4]), int(x[5])))

    ###################
    #time bin data and find median
    data = np.column_stack((time, s))
    ## if data in first time bin
    res = []
    # end of first bin:
    binstart = x0
    #binstart = datetime(2010, 01, 01, 00, 00, 00)
    delta = timedelta(days = days)
    mid_bin = binstart + timedelta(days = days/2.)
    res.append([mid_bin, data[0][1]])
    # iterate through the data item
    for d in data:
#        mid_bin = binstart + timedelta(days = days/2.)
        # if the data item belongs to this bin, append it into the bin
        if d[0] < binstart + delta:
            res[-1].append(d[1])
            continue
        # otherwise, create new empty bins until this data fits into a bin
        binstart += delta
        mid_bin = binstart + timedelta(days = days/2.)
        while d[0] > binstart + delta:
            res.append([mid_bin, []])
            binstart += delta
        # create a bin with the data
#        mid_bin = binstart + timedelta(days = days/2.)
        res.append([mid_bin, d[1]])

    bin_med = []
    bin_t = []
    bin_sig = []
    for i in range(len(res)):
        bin_med.append(np.mean(res[i][1:-1]))
        bin_t.append(res[i][0])
        bin_sig.append(np.std(res[i][1:-1]))

############################################################
    if regions == True:
        fig = plt.figure(figsize = (18,10))
        plt.subplot(111)
        plt.tick_params(axis='both', which='major', labelsize=18)
        plt.tick_params(axis='both', which='both', length = 5, width = 1)
        plt.xlabel('time', fontsize = 28)
        plt.ylabel('log(stress drop)', fontsize = 28)
        plt.xlim(x0,x1)
        plt.errorbar(x = bin_t, y = bin_med, yerr = bin_sig, marker = '.', color = 'black', lw = 2,capsize=3, capthick = 2, zorder = 5000)
        for i in range(len(time)):
            plt.errorbar(x=time[i], y=s[i], yerr = sig[i], marker = '.', c = 'black')
        #plt.scatter(x=time, y=s, c = colors, edgecolors = colors)

#        color = ['red', 'blue', 'orange', 'green', 'magenta']
#        labels = ['Laguna Salada', 'Cerro Prieto', 'Imperial Fault', 'Elsinor', 'San Jacinto']
        for i in range(len(vertex_list)):
            print vertex_list[i]
            print labels[i]
            points = zip(lon,lat)
            path = mpltPath.Path(vertex_list[i])
            inside = path.contains_points(points)
            colors = color[i]
            t = []
            s = []
            sig = []
            for j in range(len(inside)):
                if inside[j] == True:
                    x = (sobj.evid[j]).split('_')
                    t.append(datetime(int(x[0]), int(x[1]), int(x[2]), int(x[3]), int(x[4]), int(x[5])))
                    s.append(sobj.log_stressdrop[j])
                    sig.append(sobj.log_stressdrop_sig[j])
            print len(t)
            for j in range(len(t)):
                plt.errorbar(x=t[j], y=s[j], yerr = sig[j], marker = '.', c = colors)
            plt.scatter(x=t, y=s, c = colors, edgecolors = colors, label = labels[i])
        #add color key
        plt.legend(loc='lower left',ncol=3, scatterpoints=1)
####################################################################3
    else:
        cmap = mpl.cm.get_cmap('viridis')
        normalize = mpl.colors.Normalize(vmin=min(m), vmax=max(m))
        colors = [cmap(normalize(value)) for value in m]
        s_m = mpl.cm.ScalarMappable(cmap = cmap, norm=normalize)
        s_m.set_array([])       

        fig = plt.figure(figsize = (18,10))
        plt.subplot(111)
        plt.tick_params(axis='both', which='major', labelsize=18)
        plt.tick_params(axis='both', which='both', length = 5, width = 1)
        for i in range(len(time)):
            plt.errorbar(x=time[i], y=s[i], yerr = sig[i], marker = '.', c = colors[i])
        plt.scatter(x=time, y=s, c = colors, edgecolors = colors)
    #plt.plot_date(x = bin_t, y = bin_med, ls = '-', color = 'orange', lw = 3, zorder = 5000)
        plt.errorbar(x = bin_t, y = bin_med, yerr = bin_sig, marker = '.', color = 'orange', lw = 2,capsize=3, capthick = 2, zorder = 5000)

        plt.xlabel('time', fontsize = 28)
        plt.ylabel('log(stress drop)', fontsize = 28)
        plt.xlim(x0,x1)


        fig.subplots_adjust(right=0.9)
        cbar_ax = fig.add_axes([0.92, 0.15, 0.02, 0.75])
        cbar = plt.colorbar(s_m, cax=cbar_ax)
        cbar.set_label(ur"magnitude", fontsize = 22)#"$azimuth$ (\u00b0)"
        cbar.ax.tick_params(labelsize = 18)

    plt.show()
    if save == True:
        plt.savefig(top_dir + '/source_params/polygon/stressdrop_vs_time_' + name + '.png')

stressdrop_vs_time(ev = sobj.evid, lon = sobj.lon, lat = sobj.lat, s = sobj.log_stressdrop, sig = sobj.log_stressdrop_sig, m = sobj.m_cat, name = 'Regions', x0 = datetime(2010, 1, 1, 0, 0, 0),  x1 = datetime(2017, 1, 1, 0, 0, 0), save = True, days = 60, regions = True)#

def points_in_poly(vertices, s, sig, name, lat, lon, m, d, x0, x1, save, days): 
    '''
    finds points in a polygon from vertices 
    plots polygon on stressdrop map
    plots stressdrop versus time and stressdrop vs mag/depth
    Input:
        vertices:   array of vertex points of polygon
        s:          list/array of event log stress drops
        sig:        list/array of event log stress drops sigma
        name:       string name of polygon
        lat:        list/array of event latitudes
        lon:        list/array of event longitudes
        m:          list/array of event magnitudes
        d:          list/array of event depths (km)
        x0:         utc datetime, first in range for plotting
        x1:         utc datetime, first in range for plotting
        save:       True or False to save
    Output:
        plots 2 scatter plots of stress drop vs mag and stress drop versus depth
    '''
    import matplotlib.pyplot as plt
    import matplotlib as mpl
    from shapely.geometry.polygon import Polygon
    import matplotlib.path as mpltPath
    import cartopy.crs as ccrs
    import cartopy.io.img_tiles as cimgt
    
    poly = Polygon(vertices)
    x,y = poly.exterior.xy
    
    cmap = mpl.cm.get_cmap('gist_rainbow')
    normalize = mpl.colors.Normalize(vmin=min(s), vmax=max(s))
    colors = [cmap(normalize(value)) for value in s]
    s_m = mpl.cm.ScalarMappable(cmap = cmap, norm=normalize)
    s_m.set_array([])
    
    stamen_terrain = cimgt.StamenTerrain(desired_tile_form="L")
    fig = plt.figure(figsize = (18,16))
    ax = plt.axes(projection=stamen_terrain.crs)
    ax.set_extent([-118.01, -114.99, 32.29, 34.41])
    
    ax.add_image(stamen_terrain, 10, alpha = 0.5, cmap = 'gray')
    plt.tick_params(axis='both', which='major', labelsize=18)
    plt.tick_params(axis='both', which='both', length = 5, width = 1)
    
    xticks = [-118, -117, -116, -115]
    yticks = [32.3, 33, 33.5, 34, 34.4]
    ax.gridlines(xlocs=xticks, ylocs=yticks, draw_labels = True)
    
    for segment_i in range(len(fault_segments)):
        fault=fault_segments[segment_i]
        fault_z=np.zeros(len(fault))
        ax.plot(fault[:,0],fault[:,1],fault_z,color='k', transform=ccrs.Geodetic())
    
    plt.scatter(lon, lat, marker='o', color= colors, s=25, alpha=1.0, transform=ccrs.Geodetic())
    plt.plot(x, y, color = 'blue', transform=ccrs.Geodetic())
#    geodetic_transform = ccrs.Geodetic()._as_mpl_transform(ax)
        
    fig.subplots_adjust(right=0.88)
    cbar_ax = fig.add_axes([0.92, 0.15, 0.02, 0.75])
    cbar = plt.colorbar(s_m, cax=cbar_ax)
    cbar.set_label(ur"stressdrop", fontsize = 22)#"$azimuth$ (\u00b0)"
    cbar.ax.tick_params(labelsize = 18)
    plt.savefig(top_dir + '/source_params/polygon/polygon_' + name + '.png')
    plt.show()
    
    points = zip(lon,lat)
    path = mpltPath.Path(vertices)
    inside = path.contains_points(points)
    
    ev = []
    t = []
    s = []
    m = []
    sig = []
    d = []
    
    for i in range(len(inside)):
        if inside[i] == True:
            x = (sobj.evid[i]).split('_')
            x = (datetime(int(x[0]), int(x[1]), int(x[2]), int(x[3]), int(x[4]), int(x[5])))
            #if in time range
            if (x >= x0) and (x <= x1):
                x = (sobj.evid[i]).split('_')
                ev.append(sobj.evid[i])
                t.append(datetime(int(x[0]), int(x[1]), int(x[2]), int(x[3]), int(x[4]), int(x[5])))
                s.append(sobj.log_stressdrop[i])
                m.append(sobj.m_cat[i])
                sig.append(sobj.log_stressdrop_sig[i])
                d.append(sobj.depth[i])
                
    stressdrop_vs_time(ev = ev, lat = lat, lon = lon, s = s, sig = sig, m = m, name = name, x0 = x0, x1 = x1, save = save, days = days, regions = False)
    stressdrop_vs_time(ev = ev, lat = lat, lon = lon, s = s, sig = sig, m = m, name = name + '_entire_window', x0 = datetime(2010, 01, 01, 0, 0, 0), x1 = datetime(2017, 01, 01, 0, 0, 0), save = save, days = 90, regions = False)
    stressdrop_vs_depth_mag(s = s, d = d, mag = m, sig = sig, name = name, save = save)


points_in_poly(vertices = vertices_ImperialFault, s = sobj.log_stressdrop, sig = sobj.log_stressdrop_sig, name = 'ImperialFault', lat = sobj.lat, lon = sobj.lon, m = sobj.m_cat, d = sobj.depth, x0 = datetime(2012, 8, 15, 0, 0, 0), x1 = datetime(2012, 9, 1, 0, 0, 0), save = True, days = 1)
points_in_poly(vertices = vertices_CerroPrieto, s = sobj.log_stressdrop, sig = sobj.log_stressdrop_sig, name = 'CerroPrieto', lat = sobj.lat, lon = sobj.lon, m = sobj.m_cat, d = sobj.depth, x0 = datetime(2010, 1, 1, 0, 0, 0), x1 = datetime(2017, 1, 1, 0, 0, 0), save = True, days = 120)
points_in_poly(vertices = vertices_LagunaSalada, s = sobj.log_stressdrop, sig = sobj.log_stressdrop_sig, name = 'LagunaSalada', lat = sobj.lat, lon = sobj.lon, m = sobj.m_cat, d = sobj.depth, x0 = datetime(2010, 3, 15, 0, 0, 0), x1 = datetime(2010, 5, 1, 0, 0, 0), save = True, days = 1)
points_in_poly(vertices = vertices_Elsinor, s = sobj.log_stressdrop, sig = sobj.log_stressdrop_sig, name = 'Elsinor', lat = sobj.lat, lon = sobj.lon, m = sobj.m_cat, d = sobj.depth, x0 = datetime(2010, 1, 1, 0, 0, 0), x1 = datetime(2017, 1, 1, 0, 0, 0), save = True, days = 120)
points_in_poly(vertices = vertices_SanJacinto, s = sobj.log_stressdrop, sig = sobj.log_stressdrop_sig, name = 'SanJacinto', lat = sobj.lat, lon = sobj.lon, m = sobj.m_cat, d = sobj.depth, x0 = datetime(2010, 1, 1, 0, 0, 0), x1 = datetime(2017, 1, 1, 0, 0, 0), save = True, days = 120)
points_in_poly(vertices = vertices_SanAndreas, s = sobj.log_stressdrop, sig = sobj.log_stressdrop_sig, name = 'SanAndreas', lat = sobj.lat, lon = sobj.lon, m = sobj.m_cat, d = sobj.depth, x0 = datetime(2010, 1, 1, 0, 0, 0), x1 = datetime(2017, 1, 1, 0, 0, 0), save = True, days = 120)
points_in_poly(vertices = vertices_Riverside, s = sobj.log_stressdrop, sig = sobj.log_stressdrop_sig, name = 'Riverside', lat = sobj.lat, lon = sobj.lon, m = sobj.m_cat, d = sobj.depth, x0 = datetime(2010, 1, 1, 0, 0, 0), x1 = datetime(2017, 1, 1, 0, 0, 0), save = True, days = 120)
points_in_poly(vertices = vertices_Riverside, s = sobj.log_stressdrop, sig = sobj.log_stressdrop_sig, name = 'other1', lat = sobj.lat, lon = sobj.lon, m = sobj.m_cat, d = sobj.depth, x0 = datetime(2010, 1, 1, 0, 0, 0), x1 = datetime(2017, 1, 1, 0, 0, 0), save = True, days = 120)



def match_residuals(evid, mcat, slon, slat, sd, s, rlon, rlat, rdepth, eresid):
    from scipy.stats import pearsonr
    import statsmodels
    import statsmodels.stats.power
    
    
    eventid = []
    stressdrop = []
    residual = []
    for i in range(len(evid)):
        d = []
        print (evid[i])
        for j in range(len(rlat)):
            dist = dread.compute_rrup(slon[i], slat[i], sd[i], rlon[j], rlat[j], -1*rdepth[j]) #in km
            d.append(dist)
        if min(d) < 20:
            ind = d.index(min(d))
            eventid.append(evid[i])
            stressdrop.append(s[i])
            residual.append(eresid[ind])
    cmap = mpl.cm.get_cmap('viridis')
    
    normalize = mpl.colors.Normalize(vmin=min(mcat), vmax=max(mcat))
    colors = [cmap(normalize(value)) for value in mcat]
    s_m = mpl.cm.ScalarMappable(cmap = cmap, norm=normalize)
    s_m.set_array([])
    
    fig = plt.figure(figsize = (16,14))
    plt.subplot(111)
    plt.tick_params(axis='both', which='major', labelsize=18)
    plt.tick_params(axis='both', which='both', length = 5, width = 1)
    plt.scatter(x= stressdrop, y=residual, c = colors, edgecolors = colors)
    plt.xlabel('log stress drop (MPa)', fontsize = 28)
    plt.ylabel('GMPE event residual', fontsize = 28)
    plt.gca().xaxis.set_major_formatter(mtick.FormatStrFormatter('%.2f'))
    
    
    rval, pval = pearsonr(x= stressdrop, y=residual)
    power = statsmodels.stats.power.tt_solve_power(effect_size = rval, nobs = len(stressdrop), alpha = 0.05)
    label1 = 'Pearson R : ' + "{:.4f}".format(rval)
    label3 = 'power : ' + "{:.4f}".format(power)
    label2 = 'pvalue: ' + "{:.4f}".format(pval)
    
    plt.annotate(label1 + '\n' + label2 + '\n' + label3, xy=(0.72, 0.02), xycoords='axes fraction', bbox=dict(facecolor='white', edgecolor='black', boxstyle='square,pad=0.2'))

    
    fig.subplots_adjust(right=0.9)
    cbar_ax = fig.add_axes([0.92, 0.15, 0.02, 0.75])
    cbar = plt.colorbar(s_m, cax=cbar_ax)
    cbar.set_label(ur"magnitude", fontsize = 22)#"$azimuth$ (\u00b0)"
    cbar.ax.tick_params(labelsize = 18)
    
    
    
    plt.show()
    plt.savefig(top_dir + '/source_params/residual_vs_stressdrop_match_events.png')
    
#match_residuals(evid = sobj.evid, mcat = sobj.m_cat, slon = sobj.lon, slat = sobj.lat, sd = sobj.depth, s = sobj.log_stressdrop, rlon = robj.elon, rlat = robj.elat, rdepth = robj.edepth, eresid = robj.E_residual)  
    
    
    
def make_map_images(evid, time, lon, lat, s):
    time = []
    for i in range(len(evid)):
        x = (evid[i]).split('_')
        time.append(datetime(int(x[0]), int(x[1]), int(x[2]), int(x[3]), int(x[4]), int(x[5])))
    #2010
    #data = np.column_stack((time[0:2100], a['lon'][0:2100], a['lat'][0:2100], a['log_stressdrop'][0:2100]))
    #elmayor, march20 - May20
#    data = np.column_stack((time[89:1550], a['lon'][89:1550], a['lat'][89:1550], a['log_stressdrop'][89:1550]))
#    data = np.column_stack((time[89:1550], lon[89:1550], lat[89:1550], s[89:1550]))
    data = np.column_stack((time, lon, lat, s))


    ## if data in first time bin
    res = []
    # end of first bin:
#    binstart = datetime(2010, 03, 20, 00, 00, 00)
    binstart = datetime(2010, 01, 01, 00, 00, 00)
    res.append([binstart, data[0]])
    delta = timedelta(hours = 6)
    
    # iterate through the data item
    for d in data:
        #print d[0]
        # if the data item belongs to this bin, append it into the bin
        if d[0] < binstart + delta:
            res[-1].append([d[0], d[1], d[2], d[3]])
            continue
    
        # otherwise, create new empty bins until this data fits into a bin
        binstart += delta
        while d[0] > binstart + delta:
            res.append([binstart, []])
            binstart += delta
    
        # create a bin with the data
        res.append([binstart, [d[0], d[1], d[2], d[3]]])

    x = [0]
    y = [0]
    label = []
    sdrop = [-1]
    
    cmap = mpl.cm.get_cmap('gist_rainbow')
    normalize = mpl.colors.Normalize(vmin=min(s), vmax= max(s))
    colors = [cmap(normalize(value)) for value in s]
    s_m = mpl.cm.ScalarMappable(cmap = cmap, norm=normalize)
    s_m.set_array([])
    
    stamen_terrain = cimgt.StamenTerrain(desired_tile_form="L")
    fig = plt.figure(figsize = (18,16.2))
    plt.tight_layout()
    ax = plt.axes(projection=stamen_terrain.crs)
    #entire region
    #ax.set_extent([-118.01, -114.99, 32.29, 34.41])
    #2010 region
    #ax.set_extent([-116.51, -114.99, 32.29, 33.21])
    #elMayor region
    ax.set_extent([-116.2, -115.1, 32.29, 33.0])
    
    c = [cmap(normalize(value)) for value in sdrop]
    
    ax.add_image(stamen_terrain, 10, alpha = 0.5, cmap = 'gray', zorder = 3)
    plt.tick_params(axis='both', which='major', labelsize=18)
    plt.tick_params(axis='both', which='both', length = 5, width = 1)
    plt.title(label, fontsize = 26)
    
    #entire region
#    xticks = [-118, -117, -116, -115]
#    yticks = [32.3, 33, 33.5, 34, 34.4]
    #2010
    xticks = [-117,-116.5, -116, -115.5, -115, -114]
    yticks = [32, 32.3, 32.5, 32.7, 32.9, 33.1, 33.3]
    ax.gridlines(xlocs=xticks, ylocs=yticks, draw_labels = True)
    
    for segment_i in range(len(fault_segments)):
        fault=fault_segments[segment_i]
        fault_z=np.zeros(len(fault))
        ax.plot(fault[:,0],fault[:,1],fault_z,color='k', transform=ccrs.Geodetic())
    
    ax.scatter(x, y, marker='o', color= c, s=50, alpha=1.0, transform=ccrs.Geodetic())
#    geodetic_transform = ccrs.Geodetic()._as_mpl_transform(ax)
    
    fig.subplots_adjust(right=0.88)
    cbar_ax = fig.add_axes([0.92, 0.1, 0.02, 0.8])
    cbar = plt.colorbar(s_m, cax=cbar_ax)
    cbar.set_label(ur"stressdrop", fontsize = 22)#"$azimuth$ (\u00b0)"
    cbar.ax.tick_params(labelsize = 18)
    
    #for each bin
    for i in range(len(res)):
        x = [0]
        y = [0]
        #label = []
        sdrop = []
        #print i
        binned_data = res[i][1:]
        label = res[i][0]
        if binned_data != [[]]:
            for j in range(len(binned_data)):
                x.append(binned_data[j][1])
                y.append(binned_data[j][2])
                #label.append(binned_data[i][0])
                sdrop.append(binned_data[j][3])
        c = [cmap(normalize(value)) for value in sdrop]       
        ax.set_title(label, fontsize = 26, y = 1.05)
        ax.scatter(x, y, marker='o', color= c, s=60, alpha=1.0, transform=ccrs.Geodetic())    
        #plt.show()
        plt.savefig(top_dir + '/source_params/map_images_delta6hr/img' + str(i).zfill(5) + '.png')

make_map_images(evid = sobj.evid, time = sobj.time, lon = sobj.lon, lat = sobj.lat, s = sobj.log_stressdrop)

def bin_residuals_stressdrops(dlat, dlon, stat):
    #bins the stressdrops and event residuals
    #takes in dlon, dlat for grid spacing
    import matplotlib.ticker as mtick
    from scipy.stats import pearsonr
    import statsmodels
    import statsmodels.stats.power
    
    #read in lat, lon, depth, event residual
    #robj.elon, robj.elat, robj.edepth, robj.E_residual
    #read in lat, lon, depth, stressdrop, stressdrop sigma
    #a['lon'], a['lat'], a['depth'], a['log_stressdrop'], a['log_stressdrop_sig']
    #define bin edges
    #assign bins
    from numpy import ones,r_,c_
    from scipy.stats import binned_statistic_dd
    
    lon_edges = np.arange(-118, -115, dlon)
    lat_edges = np.arange(32.3, 34.4, dlat)
    
    
    sample=c_[sobj.lon,sobj.lat]
    s = sobj.log_stressdrop
    bindims = [lon_edges, lat_edges]
    statistic_s,bin_edges,binnumber=binned_statistic_dd(sample,s,statistic=stat,bins=bindims)
    
    r = robj.E_residual
    sample = c_[robj.elon, robj.elat]
    
    bindims = [lon_edges, lat_edges]
    statistic_r,bin_edges,binnumber=binned_statistic_dd(sample,r,statistic=stat,bins=bindims)
    
    a1 = np.ndarray.flatten(statistic_s)
    a2 = np.ndarray.flatten(statistic_r)
    #statistics
    x = []
    y = []
    for i in range(len(a1)):
        if str(a1[i]) != 'nan' and str(a2[i]) != 'nan':
            x.append(a1[i])
            y.append(a2[i])
    rval, pval = pearsonr(x,y)
    power = statsmodels.stats.power.tt_solve_power(effect_size = rval, nobs = len(x), alpha = 0.05)
    
    plt.figure(figsize = (16,16.))
    plt.tight_layout()
    plt.scatter(np.ndarray.flatten(statistic_r), np.ndarray.flatten(statistic_s), s = 50)
    plt.ylabel('event residual', fontsize = 26)
    plt.xlabel('stress drop', fontsize = 26)
    
    label1 = 'Pearson R : ' + "{:.4f}".format(rval)
    label3 = 'power : ' + "{:.4f}".format(power)
    label2 = 'pvalue: ' + "{:.4f}".format(pval)
    
    plt.annotate(label1 + '\n' + label2 + '\n' + label3, xy=(0.72, 0.02), xycoords='axes fraction', bbox=dict(facecolor='white', edgecolor='black', boxstyle='square,pad=0.2'))

    plt.show()
    plt.savefig(top_dir + '/source_params/residual_vs_stressdrop_' + stat + '_binned_' +str(dlon) + '_' + str(dlat)+ '.png')

    xedges, yedges = bin_edges
    X,Y = np.meshgrid(xedges, yedges)
    fig, (ax,ax2)  = plt.subplots(figsize = (24,10), ncols=2)

    cmap = plt.get_cmap('coolwarm')
    norm = mpl.colors.Normalize(vmin=min(y),vmax = max(y))

    for segment_i in range(len(fault_segments)):
        fault=fault_segments[segment_i]
        fault_z=np.zeros(len(fault))
        ax.plot(fault[:,0],fault[:,1],fault_z,color='k')
        ax2.plot(fault[:,0],fault[:,1],fault_z,color='k')

    
    c1 = ax.pcolormesh(X,Y, statistic_r.T, cmap = cmap, norm = norm)
    cbar = fig.colorbar(c1, ax=ax)
    cbar.set_label(ur"event residual", fontsize = 22)
    ax.xaxis.set_major_formatter(mtick.FormatStrFormatter('%.1f'))
    ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1f'))
    ax.set_title('GMPE event residual')
    ax.set_xlim(-118.01, -114.99)
    ax.set_ylim(32.29, 34.41)

    cmap = plt.get_cmap('coolwarm')
    norm = mpl.colors.Normalize(vmin=min(x),vmax = max(x))
    
    c2 = ax2.pcolormesh(X,Y, statistic_s.T, cmap = cmap)
    cbar2 = fig.colorbar(c2, ax=ax2)
    cbar2.set_label(ur"log stressdrop", fontsize = 22)
    ax2.xaxis.set_major_formatter(mtick.FormatStrFormatter('%.1f'))
    ax2.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1f'))
    ax2.set_title('stress drop')
    ax2.set_xlim(-118.01, -114.99)
    ax2.set_ylim(32.29, 34.41)
    plt.show()
    plt.savefig(top_dir + '/source_params/compare_' + stat + '_residual_stressdrop_binned_' +str(dlon) + '_' + str(dlat)+ '.png')
    
#bin_residuals_stressdrops(dlon = 0.075, dlat = 0.05, stat = 'mean')
#bin_residuals_stressdrops(dlon = 0.075, dlat = 0.05, stat = 'median')
    
def dist_to_elMayor(vertices, lon, lat, m, s, sig, d):
    points = zip(sobj.lon, sobj.lat)
    path = mpltPath.Path(vertices)
    inside = path.contains_points(points)
    
    t = []
    s = []
    m = []
    sig = []
    d = []
    lat = []
    lon = []
    elMayor_dist = []
    
    for i in range(len(inside)):
        if inside[i] == True:
            x = (sobj.evid[i]).split('_')
            t.append(datetime(int(x[0]), int(x[1]), int(x[2]), int(x[3]), int(x[4]), int(x[5])))
            #t.append(points_t[i])
            s.append(sobj.log_stressdrop[i])
            m.append(sobj.m_cat[i])
            sig.append(sobj.log_stressdrop_sig[i])
            d.append(sobj.depth[i])
            lat.append(sobj.lat[i])
            lon.append(sobj.lon[i])
            dist = dread.compute_rrup(elMayorlon, elMayorlat,elMayordepth, sobj.lon[i], sobj.lat[i], -1*sobj.depth[i])
            elMayor_dist.append(dist)
    
    cmap = mpl.cm.get_cmap('viridis')
    normalize = mpl.colors.Normalize(vmin=min(m), vmax=max(m))
    colors = [cmap(normalize(value)) for value in m]
    s_m = mpl.cm.ScalarMappable(cmap = cmap, norm=normalize)
    s_m.set_array([])
    
    fig = plt.figure(figsize = (18,16))
    plt.scatter(elMayor_dist, s, marker='o', s=25, c = colors, edgecolors = colors)   
    for i in range(len(s)):
            plt.errorbar(x=elMayor_dist[i], y=s[i], yerr = sig[i], marker = '.', c = colors[i])
    plt.xlabel('distance to elMayor event (km)')
    plt.ylabel('log stressdrop')
    
    fig.subplots_adjust(right=0.88)
    cbar_ax = fig.add_axes([0.92, 0.1, 0.02, 0.8])
    cbar = plt.colorbar(s_m, cax=cbar_ax)
    cbar.set_label(ur"magnitude", fontsize = 22)#"$azimuth$ (\u00b0)"
    cbar.ax.tick_params(labelsize = 18)
    
    plt.show()
    plt.savefig(top_dir + '/source_params/dist_to_elMayor_polygon1.png')

#dist_to_elMayor(vertices = [(-116.1, 32.75),(-115.7, 32.75),(-115.3, 32.4),(-115.25, 32.3),(-115.6, 32.3),(-116.1, 32.75)], lon = sobj.lon, lat = sobj.lat, m = sobj.m_cat, s = sobj.log_stressdrop, sig = sobj.log_stressdrop_sig, d = sobj.depth) 

