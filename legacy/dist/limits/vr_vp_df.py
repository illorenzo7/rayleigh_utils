# Author: Loren Matilsky
# Created: 05/06/2018
# This script calculates distribution limits for a joint pdf of v_r and v_phi
# from Meridional slices, including only the downflows (v_r < 0.)
# By default, the last 100 meridional slices are used.
# Extrema in the values in v_r and v_phi in latitude are found at each radius. 
# The output is (vr_mins, vp_maxes, vp_mins), each with shape (nr,), 
# saved in the file
# [dirname]_vr_vp_dist_limits_df_[first iter]_[last iter].npy

# Import relevant modules
import numpy as np
import sys, os
from diagnostic_reading import Meridional_Slice
from common import *

# Get the name of the run directory
dirname = sys.argv[1]
# Get the stripped name to use in file naming
dirname_stripped = strip_dirname(dirname)

# Find the relevant place to store the data, and create the directory if it
# doesn't already exist
datadir = dirname + '/data/'
if (not os.path.isdir(datadir)):
    os.makedirs(datadir)

# Directory where the Rayleigh data is kept
data_type = 'Meridional_Slices'
radatadir = dirname + '/' + data_type + '/'

file_list, int_file_list, nfiles = get_file_lists(radatadir)

# Get grid info (if it's not computed already using grid_info.py,
# this will fail)
rr,tt,cost,sint,rr_depth,ri,ro,d = np.load(datadir + 'grid_info.npy')
nr = len(rr)
nt = len(tt)

# Read in CLAs (if any) to change default file range to use (by default, 
# this will be the last 100 meridional slices)
args = sys.argv[2:]
nargs = len(args)

if (nargs == 0):
    index_first, index_last = nfiles - 100, nfiles
else:
    index_first, index_last = get_desired_range(int_file_list, args)

# Set the savename by the directory, what we are saving, and first and last
# iteration files for the average
savename = dirname_stripped + '_vr_vp_dist_limits_df_' +\
    file_list[index_first] + '_' + file_list[index_last - 1] + '.npy'
savefile = datadir + savename    

# Initialize mins and maxes so that the first data point read will become
# both the min/max

# vr_max will always be 0. since we are considering downflows only.
vr_mins = np.ones(nr)*np.inf

vp_maxes = -np.ones(nr)*np.inf
vp_mins = np.ones(nr)*np.inf

# Loop over the desired file range to set the extrema
for i in range(index_first, index_last):
    print ('considering Meridional_Slices/%s for the limits ...'\
            %file_list[i])
    mer = Meridional_Slice(radatadir + file_list[i], '')
    
    niter = mer.niter
    nphi = mer.nphi

    # Loop through each phi-value in the Meridional_Slice and each
    # time in it.
    for pindex in range(nphi):
        for tindex in range(mer.niter):

            # At each t and phi, find the v_r and v_phi values in the 
            # meridional plane; label these arrays with suffix "_loc" 
            vr_loc = mer.vals[pindex,:,:,mer.lut[1],tindex]/100. # report
                # velocities in m/s instead of cm/s
            vp_loc = mer.vals[pindex,:,:,mer.lut[3],tindex]/100.
            
            # Get the indices of the downflows, and only use those when 
            # determining the max/min in the meridional plane
            where_df = np.where(vr_loc < 0)

            vr_df = np.zeros_like(vr_loc)
            vr_df[where_df] = vr_loc[where_df]

            vp_df = np.zeros_like(vp_loc)
            vp_df[where_df] = vp_loc[where_df]

            # Compute the "local" (as in, local to this particular t and phi)
            # min/max in v_r and v_phi for the downflows, treating each radius
            # separately. If they exceed/are less than the previous min/max
            # values, update the maxes/mins.
            vr_mins_loc = np.min(vr_df, axis=0)

            vp_maxes_loc = np.max(vp_df, axis=0)
            vp_mins_loc = np.min(vp_df, axis=0)

            for ir in range(nr):
                if (vp_maxes_loc[ir] > vp_maxes[ir]):
                    vp_maxes[ir] = vp_maxes_loc[ir]

                if (vr_mins_loc[ir] < vr_mins[ir]):
                    vr_mins[ir] = vr_mins_loc[ir]
                if (vp_mins_loc[ir] < vp_mins[ir]):
                    vp_mins[ir] = vp_mins_loc[ir]
                    
print ('Saving file at ' + savefile + ' ...')
np.save(savefile, (vr_mins, vp_maxes, vp_mins))
