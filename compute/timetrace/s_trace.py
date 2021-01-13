# Author: Loren Matilsky
# Created: well before 05/06/2019
# NEEDS UPDATED DOCS AND TESTING

import numpy as np
import sys, os
from diagnostic_reading import ShellAverage
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
data_type = 'Shell_Avgs'
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
    index_first, index_last = nfiles - 100, nfiles # By default average over the                                           # last 100 files
else:
    index_first, index_last = get_desired_range(int_file_list, args)

# Set the savename by the directory, what we are saving, and first and last
# iteration files for the average
savename = dirname_stripped + '_s_vs_time_' + file_list[index_first] + '_' +\
    file_list[index_last - 1] + '.npy'
savefile = datadir + savename    

# Average over the relevant data range, summing everything and then dividing
#   by the number of "slices" at the end
s_vals = []
print ('Considering Shell_Avg files %08i through %08i for the average ...'\
       %(int_file_list[index_first], int_file_list[index_last - 1]))
count = 0
for ii in range(index_first, index_last):
    print ('Adding Shell_Avgs/%s to the average...' %file_list[ii])
    a = ShellAverage(radatadir + file_list[ii], '')
    local_ntimes = a.niter
    for jj in range(local_ntimes):
        entropy_loc = a.vals[:, 0, a.lut[64],jj]
        s_vals.append(np.mean(entropy_loc))

print ('Saving entropy_trace at ' + savefile + ' ...')
np.save(savefile, np.array(s_vals))
