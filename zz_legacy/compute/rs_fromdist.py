# Author: Loren Matilsky
# Created: 05/08/2018
# This script computes three components of the Reynolds stress tensor due 
# to the upflows and downflows separately for a Rayleigh run (stored in dir-
# name.) The script uses the run's full joint distribution of v_r and v_phi
# (produced by script dist/vr_vp_dist.py) and computes the 
# Reynolds stresses at each point in the meridional plane. 
# Output is saved in the same format  as in the output from rs.py, but only for
# the v_r/v_phi terms:
# (vr2_p, vp2_p, vrvp_p,     vr2_m, vp2_m, vrvp_m,     fplus, fminus)
# in the file [dirname]_rs_fromdist_[first iter]_[last iter].npy
# Units are in (cm/s)^2
# Quantities associated with v_phi have v_phi's deviation from the azimuthal
# mean (i.e., the differential rotation has been subtracted out)

# Import relevant modules
import numpy as np
import os, sys
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

# Get grid info (if it's not computed already using grid_info.py, this will fail)
rr,tt,cost,sint,rr_depth,ri,ro,d = np.load(datadir + 'grid_info.npy')
nr = len(rr)
nt = len(tt)

# Read in distribution
dist_file = get_widest_range_file(datadir, 'vr_vp_dist_full')
print ('About to read in distribution from ' + datadir + dist_file + ' ...')
dist = np.load(datadir + dist_file)

# Get the iterations for the file names of the saved data
iter1, iter2 = get_iters_from_file(dist_file)
savename = dirname_stripped + '_rs_fromdist_' +\
    str(iter1).zfill(8) + '_' + str(iter2).zfill(8) + '.npy'

# Get vavg file and data
vavg_file = get_widest_range_file(datadir, 'vavg')
print ('getting vavg from %s...' %(datadir + vavg_file))
vr_av, vt_av, vp_av = np.load(datadir + vavg_file)

# Filling factors for upflows and downflows
fplus = np.zeros_like(vp_av)
fminus = np.zeros_like(vp_av)

vr2_p = np.zeros_like(vp_av)
vr2_m = np.zeros_like(vp_av)

vp2_p = np.zeros_like(vp_av)
vp2_m = np.zeros_like(vp_av)


vrvp_p = np.zeros_like(vp_av)
vrvp_m = np.zeros_like(vp_av)
    
# Get bin info:
vr_bincenters, vp_bincenters = np.load(datadir + 'vrvp_bincenters.npy')
# convert back to cm/s
vr_bincenters *= 100
vp_bincenters *= 100

# Do a loop over the meridional plane, computing Reynolds stress at each 
# point
#for ir in range(20,21): # for debugging purposes
for ir in range(nr):
    # Get a 2D grid of the v_r/v_phi values associated with this particular
    # radius. Indexing='ij' forces the 2D values to match the way in which 
    # dist indexes its "counts" (i.e., the first index indicates the value of
    # v_r, the second indicates the value of v_phi)
    vrvals, vpvals = np.meshgrid(vr_bincenters[ir], vp_bincenters[ir],\
                                 indexing='ij')
                                 
    where_df, where_uf = np.where(vrvals < 0), np.where(vrvals > 0)
    
    
#    for it in range(384,385): # for debugging purposes
    for it in range(nt):
        dist_loc = dist[it, ir, :, :]
        prob = dist_loc/np.sum(dist_loc) # Make a probability distribution out
            # the "counts"
        prob_df, prob_uf = np.zeros_like(prob), np.zeros_like(prob)
        prob_df[where_df] = prob[where_df] #
        prob_uf[where_uf] = prob[where_uf] # probabilities "normalized" so 
            # that prob_df + prob_uf = prob. I.e., masking various quantities
            # by prob_df and prob_uf breaks up the contribution of that
            # quantity into contributions from upflows and downflows
        
        # Compute the various contributions from upflows and downflows to the
        # filling factors and Reynolds stress quantities, by averaging over
        # the distribution
        fplus[it, ir] = np.sum(prob_uf)
        fminus[it, ir] = np.sum(prob_df)
        
        vr_fluc = vrvals - vr_av[it, ir]
        vp_fluc = vpvals - vp_av[it, ir]
        
        vr2_p[it, ir] = np.sum(prob_uf*vr_fluc**2)
        vr2_m[it, ir] = np.sum(prob_df*vr_fluc**2)
        
        vp2_p[it, ir] = np.sum(prob_uf*vp_fluc**2)
        vp2_m[it, ir] = np.sum(prob_df*vp_fluc**2)
        
        vrvp_p[it, ir] = np.sum(prob_uf*vr_fluc*vp_fluc)
        vrvp_m[it, ir] = np.sum(prob_df*vr_fluc*vp_fluc)
        
savefile = datadir + savename
print ('Saving ' + savefile + ' ...')
np.save(savefile, (vr2_p, vp2_p, vrvp_p, vr2_m, vp2_m, vrvp_m, fplus, fminus))
