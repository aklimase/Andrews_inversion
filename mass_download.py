'''
(1) Get a cataloge of events then (2) mass download waveforms within a certain 
distance of each one of those events with obspy's mass_downloader() framework

Written by DMM, 2016; Modified by VJS, 6/2017
'''

from os import path,makedirs,remove
from obspy.clients.fdsn import Client
from obspy.core import quakeml
from obspy import UTCDateTime
from numpy import zeros,genfromtxt
from matplotlib import pyplot as plt
import milne
from obspy.clients.fdsn.mass_downloader import CircularDomain, Restrictions, MassDownloader


################         What do you want to do?            ###################

# Download catalog?
get_catalog=False
# Make a human readable format file for the catalog?
make_events_file=False
# Make files for pick times, one file per event, with stations recording it?
make_pick_files=True
# Download waveforms?
mass_download=False

###############################################################################

# In case things choke up, restart here at event number:
hot_start=0

#What data center and where should I save it dofus?
data_center='USGS'

catalog_file='/Users/escuser/project/catalog/all_6.quakeml'
events_file='/Users/escuser/project/catalog/all_6.txt'
# This folder must exist:
pick_file_path='/Users/escuser/project/catalog/all_6'

#Which channels?
channels=["HNZ", "HLZ", "HHZ", "BHZ"]

# Which stations do you want to get events on?
stations = ['FRD','RDM']

################################################################################
        
        
#load the catalog

catalog=quakeml.readQuakeML(catalog_file)
print 'Read catalog'

# If you're making pick files...
if make_pick_files:
   
    # Read in the events from the human readable text file; define the pick client
    events=genfromtxt(events_file,usecols=5,dtype='S')
    pick_client = Client(data_center)
    
    # For every event in the catalog, loop through to find the event id...
    for kevent in range(hot_start,len(events)):
        
        id_event=events[kevent]
        
        #if 'ci' in id_event: #SoCla event, get picks, otehrwise ignore
        ## We don't need the above line since all are in teh socal catalog, so instead
        #    add ci to all events
        #id_event = 'ci' + id_event
         
        # Make a pick file:   
        print 'Creating pick file for event %d of %d' %(kevent,len(events))
        
        # Make the pick file path, open the file, and write the header. 
        pick_file=pick_file_path+'/'+id_event+'.pick'
        f=open(pick_file,'w')
        f.write('#pick time, network, station ,channel\n')
                
        # Replace the event id starting with 'ci' with nothing, if 'ci' is in it:
        #id_event=id_event.replace('ci','')
        
        print id_event
        
        # Get the event, station, and pick info from the client
        ev=pick_client.get_events(eventid=id_event,includearrivals=False)
        
        # Get just the picks
        picks=ev[0].picks
        
        # If the stations of interest are included for this event, this will change
        station_counter=0
        
        # For every pick in the lengt of picks for this, get the network, station, channel, etc:
        for kpick in range(len(picks)):
            pick=picks[kpick]
            time=pick.time
            net=pick.waveform_id.network_code
            sta=pick.waveform_id.station_code
            chan=pick.waveform_id.channel_code
            
            
            # Is this event recorded on the stations of interest?  If so, write 
            #   to file...
            if sta in stations:
                line='%s\t%s\t%s\t%s\n' % (time,net,sta,chan)
                f.write(line)
                
                # Add to the counter, since it has been recorded on stations of interest
                station_counter+=1
            
        f.close()
                
        # Is this file empty, or did it record on the stations of interest a write 
        #     to file?  If it didn't write to file/doesn't have stations of interest,
        #       remove the file.
        if station_counter==0:
            remove(pick_file)
               
        
        
