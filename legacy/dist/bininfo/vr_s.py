# Author: Loren Matilsky
# Created: 05/06/2018
# This script computes the edges of the bins to be used by dist/vr_s_df.py
# in computing the joint v_r and S' distribution. The bins are equally spaced,
# and determined from the datafile 
# [dirname]_vr_s_dist_limits_full_[first iter]_[last iter].npy, produced by
# the script dist/limits/vr_s.py. nbins_vr/nbins_s (set to 100 each today)
# set the number of bins to be used in each direction.
# The script outputs 
# (vr_binedges, s_binedges) in vrs_binedges.npy (shapes (nr, nbins_var + 1]))
# (vr_bincenters, s_bincenters) in vrs_bincenters.npy (shapes (nr, nbins_var))
# and (vr_space, s_space) in vrs_binspacing.npy (shapes (nr,))

# Import relevant modules
import numpy as np
import sys, os
from common import *

# Get the name of the run directory
dirname = sys.argv[1]

# Find the relevant place to store the data, and create the directory if it
# doesn't already exist
datadir = dirname + '/data/'
if (not os.path.isdir(datadir)):
    os.makedirs(datadir)

# Get grid info (if it's not computed already using grid_info.py, this will fail)
rr,tt,cost,sint,rr_depth,ri,ro,d = np.load(datadir + 'grid_info.npy')
nr = len(rr)
nt = len(tt)

# make many bins for vr and vp (set to 100 on 05/06/2018)
nbins_vr = 100
nbins_s = 100

# Read in the bin limits (function of radius; multiply by 1.05 to be safe)
# Note the limits are in m/s
limits_file = get_widest_range_file(datadir, 'vr_s_dist_limits_full')
print ('Getting dist limits from ' + datadir + limits_file)
vr_maxes, vr_mins, s_maxes, s_mins = np.load(datadir + limits_file)
vr_maxes *= 1.05
vr_mins *= 1.05
s_maxes *= 1.05
s_mins *= 1.05

# Compute binning structure for each radius
vr_binedges = np.zeros((nr, nbins_vr + 1))
vr_bincenters = np.zeros((nr, nbins_vr))

s_binedges = np.zeros((nr, nbins_vr + 1))
s_bincenters = np.zeros((nr, nbins_vr))

vr_space = np.zeros(nr)
s_space = np.zeros(nr)

for ir in range(nr):
    minvr = vr_mins[ir]
    maxvr = vr_maxes[ir]
    mins = s_mins[ir]
    maxs = s_maxes[ir]

    vr_binedges[ir,:] = np.linspace(minvr, maxvr, nbins_vr+1)
    vr_bincenters[ir,:] = 0.5*(vr_binedges[ir, :-1] + vr_binedges[ir, 1:])
    vr_space[ir] = (maxvr - minvr)/nbins_vr

    s_binedges[ir,:] = np.linspace(mins, maxs, nbins_s+1)
    s_bincenters[ir,:] = 0.5*(s_binedges[ir, :-1] + s_binedges[ir, 1:])
    s_space[ir] = (maxs - mins)/nbins_s

# Save binspacing and binedges to use later in plots
np.save(datadir + 'vrs_binedges.npy', (vr_binedges, s_binedges))
np.save(datadir + 'vrs_bincenters.npy', (vr_bincenters, s_bincenters))
np.save(datadir + 'vrs_binspacing.npy', (vr_space, s_space))