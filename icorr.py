#### 
# Test waveforms.py to download resp file and use it to remove response
# VJS 6/2017
###

import waveforms as wf
import os
import os.path as path
from obspy import read
import glob

############################################################################
#makes response with wf.download_response
#demeans, detrends and prefilters data
#removes instrument response

# Set parameters for resp file:
network = 'AZ'
station = 'BZN'
location = '*'
channel = 'BH*'

# Unit for time series?
tsunit = 'VEL'

# Start and end date/time for resp file:
starttime = '1998-01-01T00:00.000'
endtime = '2599-12-31T23:59.000'

# Pre-filter for instrument response correction:
#   The last number should be the nyquist frequency; second to last reduces
#   the filter towards the nyquist; first two are for the "water level"
# for BH Nyquist frequency is 20, for HH Nyquist frequency is 50
# for magnitude 4-5 the corner frequency is below 1 Hz so play around with LB
prefilt = (0.6,0.9,45,50)

# sacfile (uncorrected) path, respfile path, instr. corrected sac file path:
sacfile = '/Users/escuser/Documents/Alexis_Data/cut_sac_files_2/10701405.AZ.SND.HNZ.sac'
#respfile = 'test.resp'
respfile = network + '_' + station + '_' + 'HH*.resp'
icorr_sacfile = 'test.sac'

testfilebase = path.basename(sacfile)
network = testfilebase.split('.')[1]
station = testfilebase.split('.')[2]
channel = (testfilebase.split('.')[3])
print(channel)
#channel = 'BH*'

#first rmean and rtrend
stream = read(sacfile)
tr = stream[0]
tr.detrend(type = 'demean')#removes mean of data
tr.detrend(type = 'simple')#rtrend linear from first and last samples
tr.taper(max_percentage = 0.05, type ='hann')

#rewrite to a sac file
tr.write('temp.sac', format = 'sac')
sacfile = 'temp.sac'

############################################################################

# Get response file - downloads data from IRIS and saves to respfile:
wf.download_response(network,station,location,channel,starttime,endtime,respfile)


# Correct file - loads in uncorrected sac file, corrects w/ prefilter data,
#     and removes instrument response - saves corrected to icorr_sacfile:
wf.remove_response(sacfile,respfile,icorr_sacfile,prefilt,tsunit)



#############################################################################
#plot with obspy and compare

stream1 = read(sacfile)
tr1 = stream1[0]
dt = tr1.stats.starttime
stream1.plot(size = (800, 600), starttime=dt + 30, endtime=dt + 70)


stream2 = read(icorr_sacfile)
tr2 = stream2[0]
dt = tr2.stats.starttime
stream2.plot(size = (800, 600), starttime=dt + 30, endtime=dt + 70)
#############################################################################





### Could batch this - so could do it in two steps:
# 1.   Read in all files you have, and for every station, download a resp file
#       for all HH channels at once (so set location = '*', channel = 'HH*').
#      Set it to save all these to a resp file directory.
#
# 2.   Read in all uncorrected sac files in a directory, loop through and for each
#       one, look in the resp file directory and find the one that corresponds to 
#       your station/instrument in the loop (can use "split" to get this info, I 
#       can help if you want).  Then, use wf.remove_response inside the loop for 
#       each file, save to a directory of corrected sac files. 





boxpath = '/net/anzanas.wr.usgs.gov/data/users/alexis/ANZA_boxes/Riverside_FRD_RDM'
#events = os.listdir(boxpath + '/rawdata/')
#print(events)
eventpaths = glob.glob(boxpath + '/rawdata/*')#full path
events = []
for i in range(len(eventpaths)):
    events.append(eventpaths[i].split('/')[9])
print(events)

#read in all files to find networks and stations
networklist = []
stationlist = []

for i in range(len(events)):
    files = os.listdir(eventpaths[i])
    #get the path of every file in the event directory and split into network.staation, and channel
    for j in range(len(files)):
        base = path.basename(files[j])
        network = base.split('_')[0]
        station = base.split('_')[1]
        #if network and station not park of the list yet add
        if network not in networklist:
            networklist.append(network)
        if station not in stationlist:
            stationlist.append(station)

#make response files in response directory for each combo of networks and stations
response_path = '/net/anzanas.wr.usgs.gov/data/users/alexis/ANZA_boxes/response_files'

for network in networklist:
    for station in stationlist:
        respfile = response_path + '/' + network + '_' + station + '_' + 'HH*.resp'
        channel = 'HH*'
        wf.download_response(network,station,location,channel,starttime,endtime,respfile)
        respfile = response_path + '/' + network + '_' + station + '_' + 'BH*.resp'
        channel = 'BH*'
        wf.download_response(network,station,location,channel,starttime,endtime,respfile)
        
#read in all uncorrected sac files and loop through and for each look in response directory, then remove response and save to directory of corrected files

icorr_path = boxpath + '/icorrdata'
#make a directory for each event
for i in range(len(events)):
    if not os.path.exists(icorr_path + '/' + events[i]):
        os.makedirs(icorr_path + '/' + events[i])

#read in raw data, rmean and rtrend, find .resp file, correct and add to icorr dir
for i in range(len(events)):
    files = os.listdir(eventpaths[i])
    #get the path of every file in the event directory and split into network.station, and channel
    for j in range(len(files)):
        base = path.basename(files[j])
        network = base.split('_')[0]
        station = base.split('_')[1]
        full_channel = base.split('_')[2]
        channel = full_channel[:2] + '*'
        #find response file
        respfile = response_path + '/' + network + '_' + station + '_' + channel + '.resp'
        #first rmean and rtrend
        stream = read(files[j])
        tr = stream[0]
        tr.detrend(type = 'demean')#removes mean of data
        tr.detrend(type = 'simple')#rtrend linear from first and last samples
        tr.taper(max_percentage = 0.05, type ='hann')
        #rewrite to a sac file
        tr.write('temp.sac', format = 'sac')
        sacfile = 'temp.sac'
        icorr_sacfile = icorr_path + '/' + events[i] + '/' + network + '_' + station + '_' + full_channel + '.sac'
        wf.remove_response(sacfile,respfile,icorr_sacfile,prefilt,tsunit)#prefilt values for HH









