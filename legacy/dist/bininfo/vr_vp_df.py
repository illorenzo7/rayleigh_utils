# Author: Loren Matilsky
# Created: 05/06/2018
# This script computes the edges of the bins to be used by dist/vr_vp_df.py
# in computing the joint v_r and v_phi distribution for the downflows only.
# The bins are equally spaced, and determined from the datafile 
# [dirname]_vr_vp_dist_limits_df_[first iter]_[last iter].npy, produced by
# the script dist/limits/vr_vp_df.py. nbins_vr/nbins_vp (set to 100 each today)
# set the number of bins to be used in each direction.
# The script outputs 
# (vr_binedges, vp_binedges) in vrvp_df_binedges.npy (shapes (nr, nbins_var + 1]))
# (vr_bincenters, vp_bincenters) in vrvp_df_bincenters.npy (shapes (nr, nbins_var))
# and (vr_space, vp_space) in vrvp_df_binspacing.npy (shapes (nr,))

# Import relevant modules
import numpy as np
import sys, os
from common import get_widest_range_file

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

# make many bins for vr and vp
nbins_vr = 100
nbins_vp = 100

# Read in the bin limits (function of radius; multiply by 1.05 to be safe)
# Note the limits are in m/s
limits_file = get_widest_range_file(datadir, 'vr_vp_dist_limits_df')
print ('Getting dist limits from ' + datadir + limits_file)
vr_mins, vp_maxes, vp_mins = np.load(datadir + limits_file)
vp_maxes *= 1.05
vr_mins *= 1.05
vp_mins *= 1.05

# Compute binning structure for each radius
vr_binedges = np.zeros((nr, nbins_vr + 1))
vr_bincenters = np.zeros((nr, nbins_vr))

vp_binedges = np.zeros((nr, nbins_vp + 1))
vp_bincenters = np.zeros((nr, nbins_vp))

vr_space = np.zeros(nr)
vp_space = np.zeros(nr)

for ir in range(nr):
    minvr = vr_mins[ir]
    maxvr = 0
    minvp = vp_mins[ir]
    maxvp = vp_maxes[ir]

    vr_binedges[ir,:] = np.linspace(minvr, maxvr, nbins_vr+1)
    vr_bincenters[ir,:] = 0.5*(vr_binedges[ir, :-1] + vr_binedges[ir, 1:])
    vr_space[ir] = (maxvr - minvr)/nbins_vr

    vp_binedges[ir,:] = np.linspace(minvp, maxvp, nbins_vp+1)
    vp_bincenters[ir,:] = 0.5*(vp_binedges[ir, :-1] + vp_binedges[ir, 1:])
    vp_space[ir] = (maxvp - minvp)/nbins_vp

# Save binspacing and binedges to use later in plots
np.save(datadir + 'vrvp_df_binedges.npy', (vr_binedges, vp_binedges))
np.save(datadir + 'vrvp_df_bincenters.npy', (vr_bincenters, vp_bincenters))
np.save(datadir + 'vrvp_df_binspacing.npy', (vr_space, vp_space))