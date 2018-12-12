# Typical averaging routine for Rayleigh
# Created by: Loren Matilsky
# On: 04/18/2018
##################################################################
# This routine computes the average velocities over a period of time specified
# by the user from the AZ_Avgs data directory ("datadir") in a rayleigh run. 
# By default, the routine averages over the last 100 files of datadir, though
# the user can specify a different range like:
# -n 10 (last 10 files)
# -range iter1 iter2 (names of start and stop data files; iter2 can be "last")
# -centerrange iter0 nfiles (average about central file iter0 over nfiles)

# Import relevant modules
import numpy as np
import sys, os
sys.path.append(os.environ['rasource'] + '/post_processing')
from rayleigh_diagnostics import AZ_Avgs
from common import get_file_lists, get_desired_range, strip_dirname

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
data_type = 'AZ_Avgs'
radatadir = dirname + '/' + data_type + '/'

# Get all the file names in datadir and their integer counterparts
file_list, int_file_list, nfiles = get_file_lists(radatadir)

# Get grid info (if it's not computed already using grid_info.py, this will fail)
rr,tt,cost,sint,rr_depth,ri,ro,d = np.load(datadir + 'grid_info.npy')
nr = len(rr)
nt = len(tt)

# Read in CLAs
args = sys.argv[2:]
nargs = len(args)

if (nargs == 0):
    index_first, index_last = nfiles - 100, nfiles # By default average over the
                                            # last 100 files
else:
    index_first, index_last = get_desired_range(int_file_list, args)

# Set the savename by the directory, what we are saving, and first and last
# iteration files for the average
savename = dirname_stripped + '_vavg_' + file_list[index_first] + '_' +\
    file_list[index_last - 1] + '.npy'
savefile = datadir + savename    

# Initialize velocity averages
vr = np.zeros((nt,nr))
vt = np.zeros((nt,nr))
vp = np.zeros((nt,nr))

# Average over the relevant data range, summing everything and then dividing
#   by the number of "slices" added at the end
print ('Considering AZ_Avg files %08i through %08i for the average ...'\
       %(int_file_list[index_first], int_file_list[index_last - 1]))
count = 0
for i in range(index_first, index_last):
    print ('Adding ' + data_type + '/%s to the average ...' %file_list[i])
    az = AZ_Avgs(radatadir + file_list[i], '')
    local_ntimes = az.niter
    for j in range(local_ntimes):
        vr_loc = az.vals[:,:,az.lut[1],j]
        vt_loc = az.vals[:,:,az.lut[2],j]
        vp_loc = az.vals[:,:,az.lut[3],j]

        vr += vr_loc
        vt += vt_loc
        vp += vp_loc
        count += 1

vr /= count
vt /= count
vp /= count

# Save the avarage
print ('Saving file at ' + savefile + ' ...')
np.save(savefile, (vr,vt,vp))
