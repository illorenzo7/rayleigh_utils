# Author: Loren Matilsky
# Created: 05/06/2018
# This script takes the joint v_r and S' distribution at all latitudes (computed
# from dist/vr_s.py) and creates six subdistributions:
# 1_lowlat (+/- 15 degrees)
# 2_lowmidlat (15 to 30 degrees North and South)
# 3_midlat (30 to 45 degrees North and South)
# 4_midhighlat (45 to 60 degrees North and South)
# 5_highlat (60 to 75 degrees North and South)
# 6_superhighlat (75 to 90 degrees North and South)

import numpy as np
import sys, os
from common import strip_dirname, get_widest_range_file, get_iters_from_file

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
dist_file = get_widest_range_file(datadir, 'vr_s_dist_full')
print ('About to read in distribution from ' + datadir + dist_file + ' ...')
dist = np.load(datadir + dist_file)
# Get the iterations for the file names of the saved data
iter1, iter2 = get_iters_from_file(dist_file)

# Theta grid in latitude degrees
tt_lat = tt*180/np.pi - 90

# Average over the various six latitude intervals:

# 1, lowlat
it1, it2 = np.argmin(np.abs(tt_lat - 15)), np.argmin(np.abs(tt_lat + 15))
subdist = np.sum(dist[it1:it2, :, :, :], axis=0)
savename = dirname_stripped + '_vr_s_dist_full_1_lowlat_' +\
    str(iter1).zfill(8) + '_' + str(iter2).zfill(8) + '.npy'
print ('Saving ' + datadir + savename + ' ...')
np.save(datadir + savename, subdist)

# 2, lowmidlat
it1, it2 = np.argmin(np.abs(tt_lat - 30)), np.argmin(np.abs(tt_lat - 15))
it3, it4 = np.argmin(np.abs(tt_lat + 15)), np.argmin(np.abs(tt_lat + 30))
subdist = np.sum(dist[it1:it2, :, :, :], axis=0) +\
    np.sum(dist[it3:it4, :, :, :], axis=0)
savename = dirname_stripped + '_vr_s_dist_full_2_lowmidlat_' +\
    str(iter1).zfill(8) + '_' + str(iter2).zfill(8) + '.npy'
print ('Saving ' + datadir + savename + ' ...')
np.save(datadir + savename, subdist)

# 3, midlat
it1, it2 = np.argmin(np.abs(tt_lat - 45)), np.argmin(np.abs(tt_lat - 30))
it3, it4 = np.argmin(np.abs(tt_lat + 30)), np.argmin(np.abs(tt_lat + 45))
subdist = np.sum(dist[it1:it2, :, :, :], axis=0) +\
    np.sum(dist[it3:it4, :, :, :], axis=0)
savename = dirname_stripped + '_vr_s_dist_full_3_midlat_' +\
    str(iter1).zfill(8) + '_' + str(iter2).zfill(8) + '.npy'
print ('Saving ' + datadir + savename + ' ...')
np.save(datadir + savename, subdist)

# 4, midhighlat
it1, it2 = np.argmin(np.abs(tt_lat - 60)), np.argmin(np.abs(tt_lat - 45))
it3, it4 = np.argmin(np.abs(tt_lat + 45)), np.argmin(np.abs(tt_lat + 60))
subdist = np.sum(dist[it1:it2, :, :, :], axis=0) +\
    np.sum(dist[it3:it4, :, :, :], axis=0)
savename = dirname_stripped + '_vr_s_dist_full_4_midhighlat_' +\
    str(iter1).zfill(8) + '_' + str(iter2).zfill(8) + '.npy'
print ('Saving ' + datadir + savename + ' ...')    
np.save(datadir + savename, subdist)

# 5, highlat
it1, it2 = np.argmin(np.abs(tt_lat - 75)), np.argmin(np.abs(tt_lat - 60))
it3, it4 = np.argmin(np.abs(tt_lat + 60)), np.argmin(np.abs(tt_lat + 75))
subdist = np.sum(dist[it1:it2, :, :, :], axis=0) +\
    np.sum(dist[it3:it4, :, :, :], axis=0)
savename = dirname_stripped + '_vr_s_dist_full_5_highlat_' +\
    str(iter1).zfill(8) + '_' + str(iter2).zfill(8) + '.npy'
print ('Saving ' + datadir + savename + ' ...')    
np.save(datadir + savename, subdist)

# 6, superhighlat
it1, it2 = np.argmin(np.abs(tt_lat - 90)), np.argmin(np.abs(tt_lat - 75))
it3, it4 = np.argmin(np.abs(tt_lat + 75)), np.argmin(np.abs(tt_lat + 90))
subdist = np.sum(dist[it1:it2, :, :, :], axis=0) +\
    np.sum(dist[it3:it4, :, :, :], axis=0)
savename = dirname_stripped + '_vr_s_dist_full_6_superhighlat_' +\
    str(iter1).zfill(8) + '_' + str(iter2).zfill(8) + '.npy'
print ('Saving ' + datadir + savename + ' ...')   
np.save(datadir + savename, subdist)