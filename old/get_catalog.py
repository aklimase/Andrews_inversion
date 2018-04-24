# -*- coding: utf-8 -*-

from obspy.clients.fdsn import Client
from obspy import UTCDateTime
from numpy import zeros,genfromtxt
from matplotlib import pyplot as plt
import milne

# For downloading waveforms; if it chokes in the middle of the download, can
#   start here at waveform number hot_start
hot_start=None

####  Get a catalogue of events ####

#Time period for search (A few years)
tstart=UTCDateTime("1998-01-01")
tend=UTCDateTime("2016-12-31")

# geographic Box for search38.22°N 122.31°WCoordinates: 38.22°N 122.31°W[1]
lon_query=[-115.75000,-115.20000]
lat_query=[33.00000,33.40000]  
region_name = '4' 

minimum_magnitude = 2.5

out_file='/Users/escuser/project/catalog/all_4'
data_center='USGS'


#Run the thing and don't return an object
#  Here it is returning the catalog as an object, "catalog", so set
#   return_object=True, not False
#
#   Also want it in a JSON format, so say format='JSON', in all caps.  Available
#   formats are in the documentation for obspy.core.events.catalog class 

catalogue = milne.get_catalogue(tstart,tend,lon_query,lat_query,out_file=out_file,minmagnitude=minimum_magnitude,
        human_readable=True,data_center=data_center,return_object=True,format='QUAKEML')



## Plot the events, by color and magitude ##

# Load in data from the human readable file:
evdata = genfromtxt(out_file+'.txt',skip_header=1)
evcolor = evdata[:,3]
evsize = 10**(evdata[:,3])/100 + 5

eventfig = plt.figure()
plt.scatter(evdata[:,0],evdata[:,1],c=evcolor,s=evsize)

# Add colorbar and labels:
plt.colorbar()
plt.xlabel('Longitude')
plt.ylabel('Latitude')
plt.title('Events in region ' + region_name)


###
## If you want to download waveforms as well...
###


##Now download all the waveforms
#path_to_catalogue=out_file+'.QUAKEML'
#
#network='NC'
#station='C031'
#channels='HNZ'
#download_folder='/home/dmelgar/re1906/EEW/stations/NC.J056'
#tbefore=5.0  #Seconds before P-wave pick
#tafter=15.0  #Seconds after P-wave pick
#
#milne.get_waveforms_one_station(path_to_catalogue,station,network,channels,download_folder,
#        apply_gain=True,tbefore=5.0,tafter=15.0,data_center=data_center,quiet=True,hot_start=hot_start)