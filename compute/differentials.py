###########################
# Author: Loren Matilsky  #
# Revised on: 05/03/2018  #
#############################################################################
# Computes the differentials (dphi, dt, dr) associated with a Rayleigh grid #
#############################################################################
# Notes: uses basic grid information to create makeshift differentials      #
# (not necessarily the correct ones; must ask Nick about this.) Assumes     #
# grid info has already been calculated through "grid_info.py."             #
#############################################################################

# Import relevant modules
import numpy as np
import sys, os

# Get directory name for the Rayleigh run
dirname = sys.argv[1]

# Get the data directory
datadir = dirname + '/data/'
if (not os.path.isdir(datadir)):
    os.makedirs(datadir)

# Get basic grid info
rr,tt,cost,sint,rr_depth,ri,ro,d = np.load(datadir + 'grid_info.npy')
nr = len(rr); nt = len(tt)

# Initialize differentials
dr = np.zeros_like(rr)
dt = np.zeros_like(tt)

# Theta grid (excludes the endpoints).
# Goes "backward" (from pi to 0).
# Differentials are computed as the distance between midpoints of the theta
# grid, to try to make each theta collocation ponit "representative" of its
# associated differential. 
# The first and last differential must extend to pi and 0, respectively, since
# theta does not include the endpoints.
tt_midpoints = (tt[1:] + tt[:-1])/2
dt[0] = np.pi - tt_midpoints[0]
dt[1:-1] = tt_midpoints[:-1] - tt_midpoints[1:] 
dt[-1] = tt_midpoints[-1] - 0

# radial grids (includes the endpoints)
# Also goes "backward" (from r_o to r_i)
rr_midpoints = (rr[1:] + rr[:-1])/2
dr[0] = ro - rr_midpoints[0]
dr[1:-1] = rr_midpoints[:-1] - rr_midpoints[1:]
dr[-1] = rr_midpoints[-1] - ri

# For completeness, we calculate the dphi (uniform) differential. 
# This is just 2*pi/nphi
nphi = 2*nt
dphi = np.array([2*np.pi/nphi])

savefile = datadir + 'differentials.npy'
np.save(savefile, (dphi, dt, dr))