#Make the Anza faults file
#VJS 9/2016

import dread
import cPickle as pickle

#Name of the gmt multisegment file:
multisegpath='/Volumes/USGS_Data/project/python_codes/Holocene_LatestPleistocene.txt'

#Name of the pckl output file:
#pcklpath='/Users/vsahakian/anza/data/faults/Holocene_LatestPleistocene_117.5w_115.5w_33n_34n.txt'
pcklpath='/Volumes/USGS_Data/project/catalogs/Holocene_LatestPleistocene_118.0w_115.0w_32.3n_34.4n.pckl'


#Limits to output:
pathlimits=[[-118.0,-115.0],[32.3,34.4]]

#Convert:
allsegments=dread.multiseg2pckl(multisegpath,pcklpath,pathlimits)

#Write to the pickle file:
fout=open(pcklpath,'w')
for segment_i in range(len(allsegments)):
    pickle.dump(allsegments[segment_i],fout)
fout.close()