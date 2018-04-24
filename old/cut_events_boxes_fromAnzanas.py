####
# Make event files for all events in catalog, for all stations available
# VJS 6/2017
####

import waveforms as wf
import numpy as np
from obspy.core import UTCDateTime
from glob import glob

# Path to the event box catalog file:
catalogfile = '/net/anzanas.wr.usgs.gov/data/users/alexis/Imperial_Valley.txt'

# Main anzanas directory:
anzanas_dir = '/net/anzanas.wr.usgs.gov/data/ANZA_waveforms/'

# Acceptable channels:
accept_channel = '*H[H,N][N,Z,E]*'

# Cut event file directory:
cut_dir = '/net/anzanas.wr.usgs.gov/data/users/alexis/ANZA_boxes/Imperial_Valley_PFO_TPFO_PMD/'

#############################################################################
# Read in catalog file - in floats, and strings:
evcatalog_floats = np.genfromtxt(catalogfile,usecols=(0,1,2,3))
evcatalog_strings = np.genfromtxt(catalogfile,usecols=(4,5),dtype='S')

# Collect columns of data from floats::
evlon = evcatalog_floats[:,0]
evlat = evcatalog_floats[:,1]
evdepth = evcatalog_floats[:,2]
evmag = evcatalog_floats[:,3]

# Collect columns of string data, origin time and evend id's:
evorigin = evcatalog_strings[:,0]
evid = evcatalog_strings[:,1]

##########################
# Loop over the events
#   For every event, get the event info, and then find all files from that 
#   day and cut them.

for eventi in range(len(evid)):
    evlon_i = evlon[eventi]
    evlat_i = evlat[eventi]
    evdepth_i = evdepth[eventi]
    evmag_i = evmag[eventi]
    
    # Event origin time needs to be converted to a utcdatetime object:
    evorigin = UTCDateTime(evorigin[eventi])
    # Remove the ci from the event id
    evid = evid[eventi].split('ci')[1]
    
    # Now get the year and julian day:
    evyear_i = evorigin.year
    evjulday_i = evorigin.julday
    
    ## Now, cut every file.
    # Directory - 
    event_dir = anzanas_dir + evyear_i + '/' + evjulday_i + '/'
    
    # Glob all files here:
    globlist = glob(event_dir + '*H[H,N][N,Z,E]*')
    
    ## For each file in this list, cut the file and store to the new directory.
    # New event directory:
    event_dir_new = cut_dir + 'Event_' + evid + '/'
    
    #for each station and channel of an event, convert to sac and cut, put into event dir
    or sta_chan in range(len(globlist)):
        
    