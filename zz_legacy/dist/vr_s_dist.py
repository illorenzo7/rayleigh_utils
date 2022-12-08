# Author: Loren Matilsky
# Created: 02/25/2018
# This script computes the 2-D distribution for v_r and S' at each point in 
# the nt x nr meridional plane. 
# Requires the data produced by dist/bininfo/vr_s.py stored in 
# vrs_binedges.npy, vrs_binspacing.npy
# Outputs the distribution as (nt, nr, nbins_vr, nbins_s) (big file!) at 
# [dirname]_vr_s_dist_full_[first iter]_[last iter].npy

# Import relevant modules
import numpy as np
import sys, os
from diagnostic_reading import Meridional_Slice
from common import *
    get_widest_range_file

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

# Get grid info (if it's not computed already using grid_info.py, this will fail)
rr,tt,cost,sint,rr_depth,ri,ro,d = np.load(datadir + 'grid_info.npy')
nr = len(rr)
nt = len(tt)

# Get spherically averaged mean entropy to subtract off 
s_mean_file = get_widest_range_file(datadir, 's_spherical_mean')
s_spherical_mean = np.load(datadir + s_mean_file)
s_spherical_mean_2d = s_spherical_mean.reshape((1, nr))

# Read in CLAs (if any) to change default file range to use (by default, this
# will be the last 100 meridional slices)
args = sys.argv[2:]
nargs = len(args)
if (nargs == 0):
    index_first, index_last = nfiles - 100, nfiles
else:
    index_first, index_last = get_desired_range(int_file_list, args)
print ('Making distribution for files in the range %s to %s' 
       %(file_list[index_first], file_list[index_last - 1]))
   
# Set the savename by the directory, what we are saving, and first and last
# iteration files for the average
savename = dirname_stripped + '_vr_s_dist_full_' + file_list[index_first] + '_' +\
    file_list[index_last - 1] + '.npy'
savefile = datadir + savename  

# Get the bin information from the output of dist/bininfo/vr_s.py
vr_binedges, s_binedges = np.load(datadir + 'vrs_binedges.npy')
vr_space, s_space = np.load(datadir + 'vrs_binspacing.npy')

vr_maxes = np.max(vr_binedges, axis=1)
vr_mins = np.min(vr_binedges, axis=1)
s_maxes = np.max(s_binedges, axis=1)
s_mins = np.min(s_binedges, axis=1)

nbins_vr = len(vr_binedges[0, :]) - 1
nbins_s = len(s_binedges[0, :]) - 1

# Keep a logfile to monitor program's progress
logname = 'logfile_' + savename.split('.')[0]
print ('Saving log at ' + datadir + logname + ' ...')
logfile = open(datadir + logname, 'w') # For the FIRST TIME (outside the long
    # "for loop,") open the file to write. Within the loop, open it to append.
    # The purpose of closing and reopening the file is so the progress of the
    # script (which can take over a day) can be checked remotely.
opened_only_once = True

# Array to hold the distribution
vr_s_dist = np.zeros((nt, nr, nbins_vr, nbins_s))

for i in range(index_first, index_last):
    # Log the progress we are making!
    if (not opened_only_once):
        logfile = open(datadir + logname, 'a')
    else:
        opened_only_once = False
    logfile.write('adding %s/Meridional_Slices/%s to the distribution...\n'\
                  %(dirname.split('/')[-1],file_list[i]))
    logfile.close()    
    # Print the message for good measure
    print('adding %s/Meridional_Slices/%s to the distribution...\n'\
          %(dirname.split('/')[-1],file_list[i]))
    
    mer = Meridional_Slice(radatadir + file_list[i], '')  
    niter = mer.niter
    nphi = mer.nphi

    # Loop through each phi-value the Meridional_Slice and each
    # time in it.
    for pindex in range(nphi):
        for tindex in range(mer.niter):
            # Compute the vr and vp values for this phi-value and time
            # convert data from (cm/s) --> (m/s)
            vr_mer = mer.vals[pindex,:,:,mer.lut[1],tindex]/100.
            s_mer = mer.vals[pindex,:,:,mer.lut[64],tindex] -\
                s_spherical_mean_2d

            # Loop through the current slice in theta-r space,
            # binning each value of vr and s accordingly.
            
            for ir in range(nr):

                # Get max bin values and spacing
                minvr_loc = vr_mins[ir]
                maxvr_loc = vr_maxes[ir]
                mins_loc = s_mins[ir]
                maxs_loc = s_maxes[ir]

                vr_space_loc = vr_space[ir]
                s_space_loc = s_space[ir]

                for it in range(nt):
                    # compute the local vr/vp for  this cell in the slice
                    vr_loc = vr_mer[it,ir]
                    s_loc = s_mer[it,ir]
                    # bin vr, putting all values outside [minvr,maxvr] 
                    # in the
                    # closest bin in range (the outermost bins)
                    if (vr_loc < minvr_loc + vr_space_loc):
                        bin_index_vr = 0
                    elif (vr_loc > maxvr_loc - vr_space_loc):
                        bin_index_vr = nbins_vr - 1
                    else:
                        bin_index_vr = np.floor((vr_loc - minvr_loc)\
                                /vr_space_loc)
                        bin_index_vr = int(bin_index_vr)

                    # bin vp according to the same process
                    if (s_loc < mins_loc + s_space_loc):
                        bin_index_s = 0
                    elif (s_loc > maxs_loc - s_space_loc):
                        bin_index_s = nbins_s - 1
                    else:
                        bin_index_s = np.floor((s_loc - mins_loc)\
                                /s_space_loc)
                        bin_index_s = int(bin_index_s)
                   
                    vr_s_dist[it,ir,bin_index_vr,bin_index_s] += 1.

print ('Saving the vr_vp distribution at ' + savefile + ' ...')
np.save(savefile, vr_s_dist)