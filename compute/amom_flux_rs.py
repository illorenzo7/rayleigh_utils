# Typical averaging routine for Rayleigh
# Created by: Loren Matilsky
# On: 05/24/2018
##################################################################
# This routine computes the angular momentum fluxes due to the Reynolds stress
#  over a period of time specified
# by the user from the AZ_Avgs data directory ("datadir") in a rayleigh run. 
# By default, the routine averages over the last 100 files of datadir, though
# the user can specify a different range like:
# -n 10 (last 10 files)
# -range iter1 iter2 (names of start and stop data files; iter2 can be "last")
# -centerrange iter0 nfiles (average about central file iter0 over nfiles)

# Import relevant modules
import numpy as np
import sys, os
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
use_old_lut = False
for i in range(nargs):
    arg = args[i]
    if (arg == '-old'):
        use_old_lut = True

if (nargs == 0 or (nargs == 1 and use_old_lut == True)):
    index_first, index_last = nfiles - 100, nfiles
else:
    index_first, index_last = get_desired_range(int_file_list, args)

if (use_old_lut):
    from diagnostic_reading import AzAverage
    reading_function = AzAverage
else:
    from rayleigh_diagnostics import AZ_Avgs
    reading_function = AZ_Avgs

# Set the savename by the directory, what we are saving, and first and last
# iteration files for the average
savename = dirname_stripped + '_amomflux_rs_' + file_list[index_first] + '_' +\
    file_list[index_last - 1] + '.npy'
savefile = datadir + savename    

f_rs_r = np.zeros((nt,nr))
f_rs_t = np.zeros((nt,nr))

# Get indices for the angular momentum flux
if (use_old_lut):
    index_f_rs_r = 149
    index_f_rs_t = 150
else:
    index_f_rs_r = 1807
    index_f_rs_t = 1808
print(index_f_rs_r, index_f_rs_t)
print ('Considering AZ_Avg files %08i through %08i for the average ...'\
       %(int_file_list[index_first], int_file_list[index_last - 1]))

count = 0
for i in range(index_first, index_last):
    print ('Adding ' + data_type + '/%s to the average ...' %file_list[i])
    a = reading_function(radatadir + file_list[i], '')
    local_ntimes = a.niter
    for j in range(local_ntimes):
        f_rs_r_loc = a.vals[:, :, a.lut[index_f_rs_r], j]
        f_rs_t_loc = a.vals[:, :, a.lut[index_f_rs_t], j]

        f_rs_r += f_rs_r_loc
        f_rs_t += f_rs_t_loc
        count += 1

f_rs_r /= count
f_rs_t /= count

# Save the avarage
print ('Saving file at ' + savefile + ' ...')
np.save(savefile, (f_rs_r, f_rs_t))