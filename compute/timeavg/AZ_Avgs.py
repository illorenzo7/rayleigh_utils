# Routine to average Rayleigh AZ_Avgs data in time
# Created by: Loren Matilsky
# On: 11/10/2018
##################################################################
# This routine computes the average in time of the values in the AZ_Avgs data 
# for a particular simulation. 

# By default, the routine averages over the last 100 files of datadir, though
# the user can specify a different range in sevaral ways:
# -n 10 (last 10 files)
# -range iter1 iter2 (names of start and stop data files; iter2 can be "last")
# -centerrange iter0 nfiles (average about central file iter0 over nfiles)

# Import relevant modules
import numpy as np
import sys, os
sys.path.append(os.environ['rapp'])
sys.path.append(os.environ['co'])
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

radatadir = dirname + '/AZ_Avgs/'

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
savename = dirname_stripped + '_AZ_Avgs_' + file_list[index_first] + '_' +\
    file_list[index_last] + '.npy'
savefile = datadir + savename    

# Initialize empty "vals" array for the time average
az0 = AZ_Avgs(radatadir + file_list[index_first], '')
nr = az0.nr
nt = az0.ntheta
qv = az0.qv
nqv_tot = len(qv)
qinds = np.arange(len(qv)) # initial indices of the qv in the AZ_Avgs 
    # vals array
counts = np.zeros_like(qv) # keep track of how many slices we averaged
        # over for each variable separately (variables may be added 
        # or removed mid-simulation because we are all human)
iters1 = np.zeros_like(qv) + az0.iters[0]
iters2 = np.zeros_like(qv) + az0.iters[0]

vals = np.zeros_like(az0.vals[:, :, :, 0])

# Average over the relevant data range, summing everything and then dividing
#   by the number of "slices" added at the end
print ('Considering AZ_Avgs files %s through %s for the average ...'\
       %(file_list[index_first], file_list[index_last]))
for i in range(index_first, index_last + 1):
    print ('Adding AZ_Avgs/%s to the average ...' %file_list[i])
    if i == index_first:
        az = az0
    else:   
        az = AZ_Avgs(radatadir + file_list[i], '')

    current_qv = az.qv

    if not np.array_equal(current_qv, qv): # the outputs have changed!
        # must adjust array dimensions and qinds
        old_qv = np.copy(qv)
        nqv_before = len(old_qv)
        nqv_now = len(current_qv)
        nqv_tot = len(np.unique(np.append(old_qv, current_qv)))
        print(nqv_before, nqv_now, nqv_tot)
        print ('old_qv: ', old_qv)
        print ('new_qv: ', current_qv)
        nappended = 0
        qinds = np.zeros(nqv_now, dtype=int)
        for j in range(nqv_now):
            qv_loc = current_qv[j]
            if qv_loc in old_qv:
                qinds[j] = np.argmin(np.abs(old_qv - qv_loc))
                print('old output: ', qv_loc, j, qinds[j])
            else:
                qv = np.append(qv, qv_loc)
                dummy = np.zeros((nt, nr, 1))
                vals = np.append(vals, dummy, 2)
                counts = np.append(counts, 0)
                iters1 = np.append(iters1, az.iters[0])
                iters2 = np.append(iters2, 0)

                qinds[j] = nqv_before + nappended
                nappended += 1
                print('new output: ', qv_loc, j, qinds[j])
                print (np.shape(vals))
#                print(np.shape(vals))
#    print(qv)
#    print(qinds)
    local_ntimes = az.niter
    for j in range(local_ntimes):
        vals[:, :, qinds] += az.vals[:, :, :, j]
        counts[qinds] += 1
        iters2[qinds] = az.iters[j]

vals /= counts.reshape((1, 1, nqv_tot))
#print ('Averaged over %i AZ_Avgs slice(s) ...' %count)

# Save the avarage
print ('Saving file at ' + savefile + ' ...')
np.save(savefile, (vals, qv, counts, iters1, iters2))
