#### Compre gridded path terms with a material model
## VJS 12/2017

import dread as dr
import cPickle as pickle
import res_analysis as ra
import numpy as np
import matplotlib.pyplot as plt


####################################################################
#Change this parameter depending on where you run:
#0=desktop
#1=mac

what_home=0

if what_home==0:
    #Desktop:
    HOME='/katmai'
elif what_home==1:
    #Mac:
    HOME='/Users/vsahakian'

#############################
########### Paths ###########
#############################

home=HOME+'/anza/models/residuals/'
run_name = 'mixedregr_v3anza2013_pga_5coeff_a4_-1.20_Mc_8.5_res4_noVs30/'
residualpath=home + run_name + 'mixedcoeff_v3anza2013_pga_noVs30_5coeff_a4_-1.2_pga__ncoeff5_Mc_8.5_VR_99.6_a4_-1.2_robj_FILTERED_raydat.pckl'
materialpath=HOME+'/anza/data/pckl/FangVs.pckl'
#gobjpath=home+run_name+'/'+'mixedcoeff_v3anza2013_pga_noVs30_5coeff_a4_-1.2_pga__ncoeff5_Mc_8.5_VR_99.6_a4_-1.2_robj_FILTERED_mean_grid.pckl'
gobjpath_mean = '/home/vsahakian/Desktop/mixedcoeff_v3anza2013_pga_noVs30_5coeff_a4_-1.2_pga__ncoeff5_Mc_8.5_VR_99.6_a4_-1.2_robj_FILTERED_mean_grid.pckl'
gobjpath_count = '/home/vsahakian/Desktop/mixedcoeff_v3anza2013_pga_noVs30_5coeff_a4_-1.2_pga__ncoeff5_Mc_8.5_VR_99.6_a4_-1.2_robj_FILTERED_count_grid.pckl'
gobjpath_median = '/home/vsahakian/Desktop/mixedcoeff_v3anza2013_pga_noVs30_5coeff_a4_-1.2_pga__ncoeff5_Mc_8.5_VR_99.6_a4_-1.2_robj_FILTERED_median_grid.pckl'

statistic = 'mean'

# Set raytpe to Vs, whoich is 1:
raytype = 1

####################################################################

###############################################
########### Step 1 - Read in models ###########
###############################################
print 'Reading in material model'
# Material model:
mfile = open(materialpath,'r')
mobj = pickle.load(mfile)
mfile.close()

print 'Reading in residual model'
# Residuals model:
rfile = open(residualpath,'r')
robj = pickle.load(rfile)
rfile.close()


###############################################
########### Step 2 - Get bin edges  ###########
###############################################

print 'Getting bin edges'

# Material object has nodes - that means that there will be len(nodes) + 2 edges:
x_numedges = len(mobj.x) + 1
y_numedges = len(mobj.y) + 1
z_numedges = len(mobj.z) + 1

# Get the differences between nodes:
x_diff = np.diff(mobj.x)
y_diff = np.diff(mobj.y)
z_diff = np.diff(mobj.z)

# Initiate empty bin edges:
x_binedges = np.zeros((x_numedges))
y_binedges = np.zeros((y_numedges))
z_binedges = np.zeros((z_numedges))


# Fill first of each:
x_binedges[0] = mobj.x[0] - (x_diff[0]/2.)
y_binedges[0] = mobj.y[0] - (y_diff[0]/2.)
z_binedges[0] = mobj.z[0] - (z_diff[0]/2.)

#
for i_node in range(len(mobj.x)):
    if i_node < len(mobj.x) - 1:
        x_binedges[i_node + 1] = mobj.x[i_node] + (x_diff[i_node]/2.)
    else:
        x_binedges[i_node + 1] = mobj.x[i_node] + (x_diff[i_node - 1]/2.)
    
for i_node in range(len(mobj.y)):
    if i_node < len(mobj.y) - 1:
        y_binedges[i_node + 1] = mobj.y[i_node] + (y_diff[i_node]/2.)
    else:
        y_binedges[i_node + 1] = mobj.y[i_node] + (y_diff[i_node - 1]/2.)

for i_node in range(len(mobj.z)):
    if i_node < len(mobj.z) - 1:
        z_binedges[i_node + 1] = mobj.z[i_node] + (z_diff[i_node]/2.)
    else:
        z_binedges[i_node + 1] = mobj.z[i_node] + (z_diff[i_node - 1]/2.)

# Make binedges arrays:
all_binedges = [x_binedges,y_binedges,z_binedges]


###############################################
########### Step 3 - Compute grid   ###########
###############################################
# Use all rays:
ind_subset = 'all'

#print 'Gridding for mean...'
#
#gridded_obj_mean=ra.grid_path_term(robj,all_binedges,raytype,statistic,rpath_type='object',subset=ind_subset)
#
#print 'Saving mean grid...'
## Save grid to file:
#gfile = open(gobjpath_mean,'w')
#pickle.dump(gridded_obj_mean,gfile)
#gfile.close()

####################
print 'Gridding for counts...'
# Also get counts:
statistic = 'count'
gridded_obj_count = ra.grid_path_term(robj,all_binedges,raytype,statistic,rpath_type='object',subset=ind_subset)

print 'Saving count grid...'
# Save grid to file:
gfile = open(gobjpath_count,'w')
pickle.dump(gridded_obj_count,gfile)
gfile.close()

####################
print 'Gridding for median...'
# Also get median:
statistic = 'median'
gridded_obj_median = ra.grid_path_term(robj,all_binedges,raytype,statistic,rpath_type='object',subset=ind_subset)

print 'Saving median grid...'
# Save grid to file:
gfile = open(gobjpath_median,'w')
pickle.dump(gridded_obj_median,gfile)
gfile.close()

###############################################
###########     Step 4 - Reshape    ###########
###############################################

#print 'Reshaping...'
## First reshape the statistic:
#array_length = np.shape(gridded_obj_mean.statistic)[0] * np.shape(gridded_obj_mean.statistic)[1] * np.shape(gridded_obj_mean.statistic)[2]
#mean_reshape = np.reshape(gridded_obj_mean.statistic,array_length)
#
## Now reshape the velocity model 
##  For every depth slice, convert the x by y to a vector; then concatenate all at the end.
#for i_depth in range(np.shape(mobj.materials)[0]):
#    i_array_length = np.shape(mobj.materials)[1] * np.shape(mobj.materials)[2]
#    i_slice = np.reshape(mobj.materials[i_depth], [1,i_array_length])
#
#    # If it's the first slice, make the velocity model equal to it:
#    if i_depth == 0:
#        velmodel_reshape = i_slice
#    else:
#        velmodel_reshape = np.append(velmodel_reshape,i_slice)
#
## Also use np.sparse

###############################################
########### Step 5 - Correlate bins  ##########
###############################################


###############################################
###########     Step 6 - Plots      ###########
###############################################



