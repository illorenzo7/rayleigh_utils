# Routine to average Rayleigh Shell_Avgs data in time
# Created by: Loren Matilsky
# On: 11/12/2018
##################################################################
# This routine computes the average in time of the values in the Shell_Avgs data 
# for a particular simulation. 

# By default, the routine averages over the last 100 files of datadir, though
# the user can specify a different range in sevaral ways:
# -n 10 (last 10 files)
# -range iter1 iter2 (names of start and stop data files; iter2 can be "last")
# -centerrange iter0 nfiles (average about central file iter0 over nfiles)

# Import relevant modules
import numpy as np
import sys, os
sys.path.append(os.environ['rasource'] + '/post_processing')
sys.path.append(os.environ['co'])
from rayleigh_diagnostics import Shell_Avgs
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

radatadir = dirname + '/Shell_Avgs/'

# Get all the file names in datadir and their integer counterparts
file_list, int_file_list, nfiles = get_file_lists(radatadir)

# Read in CLAs
args = sys.argv[2:]
nargs = len(args)

if (nargs == 0):
    index_first, index_last = nfiles - 101, nfiles - 1  
    # By default average over the last 100 files
else:
    index_first, index_last = get_desired_range(int_file_list, args)

# Set the timeavg savename by the directory, what we are saving, and first and last
# iteration files for the average
savename = dirname_stripped + '_Shell_Avgs_' + file_list[index_first] + '_' +\
    file_list[index_last] + '.npy'
savefile = datadir + savename    

# Initialize empty "vals" array for the time average
sh0 = Shell_Avgs(radatadir + file_list[index_first], '')
vals = np.zeros_like(sh0.vals[:, 0, :, 0]) # 0 in second axis; get mean (not curtosis, etc.)

# Grid info
rr = sh0.radius
nr = sh0.nr
ri, ro = np.min(rr), np.max(rr)
d = ro - ri
rr_height = 

# Average over the relevant data range, summing everything and then dividing
#   by the number of "slices" added at the end
print ('Considering Shell_Avgs files %s through %s for the average ...'\
       %(file_list[index_first], file_list[index_last]))

count = 0
iter1, iter2 = int_file_list[index_first], int_file_list[index_last]

for i in range(index_first, index_last + 1):
    print ('Adding Shell_Avgs/%s to the average ...' %file_list[i])
    if i == index_first:
        sh = sh0
    else:   
        sh = Shell_Avgs(radatadir + file_list[i], '')

    local_ntimes = sh.niter
    for j in range(local_ntimes):
        vals += sh.vals[:, 0, :, j]
        count += 1

vals /= count
print ('Averaged over %i Shell_Avgs slice(s) ...' %count)

# Save the avarage
print ('Saving file at ' + savefile + ' ...')
np.save(savefile, {'vals': vals, 'lut': sh0.lut, 'count': count, 'iter1': iter1, 'iter2': iter2},\
       'rr': sh0.radius, 'nr': sh0.nr)