# Author: Loren Matilsky
# Created: 05/14/2018
# This script fixes the error in the df distributions, in which all flows 
# with vr > 0 were binned in the "vr=0" bin. 

import numpy as np
import sys, os
from common import strip_dirname, get_widest_range_file

# Get the name of the run directory
dirname = sys.argv[1]
# Get the stripped name to use in file naming
dirname_stripped = strip_dirname(dirname)

# Find the relevant place to store the data, and create the directory if it
# doesn't already exist
datadir = dirname + '/data/'
if (not os.path.isdir(datadir)):
    os.makedirs(datadir)

# Do vr_s_df first
dist_file = get_widest_range_file(datadir, 'vr_s_dist_df')
print ('About to read in distribution from ' + datadir + dist_file + ' ...')
dist = np.load(datadir + dist_file)
print ('Fixing the "vr = 0" bin ...')
dist[:, :, -1, :] = dist[:, :, -2, :]
print ('Saving the modified distribution back at ' + datadir + dist_file + ' ...')
np.save(datadir + dist_file, dist)

# Then vr_vp
dist_file = get_widest_range_file(datadir, 'vr_vp_dist_df')
print ('About to read in distribution from ' + datadir + dist_file + ' ...')
dist = np.load(datadir + dist_file)
print ('Fixing the "vr = 0" bin ...')
dist[:, :, -1, :] = dist[:, :, -2, :]
print ('Saving the modified distribution back at ' + datadir + dist_file + ' ...')
np.save(datadir + dist_file, dist)