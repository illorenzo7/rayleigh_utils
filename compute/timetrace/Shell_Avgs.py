# Routine to trace Rayleigh Shell_Avgs data in time
# Created by: Loren Matilsky
# On: 11/08/2019
##################################################################
# This routine computes the trace in time of the values in the Shell_Avgs data 
# for a particular simulation. 

# By default, the routine traces all files of datadir, though
# the user can specify a different range in sevaral ways, e.g.:
# -n 10 (last 10 files)
# -range iter1 iter2 (names of start and stop data files; can also
# be the strings "first" or "last")
# -centerrange iter0 nfiles (trace about central file iter0 over nfiles)

# Import relevant modules
import numpy as np
import pickle
import sys, os
sys.path.append(os.environ['rapp'])
sys.path.append(os.environ['raco'])
from rayleigh_diagnostics import Shell_Avgs
from common import *

# Get the name of the run directory
dirname = sys.argv[1]
# Get the stripped name to use in file naming
dirname_stripped = strip_dirname(dirname)

# Find the relevant place to store the data, and create the directory if it
# doesn't already exist
datadir = dirname + '/data/'
if not os.path.isdir(datadir):
    os.makedirs(datadir)

radatadir = dirname + '/Shell_Avgs/'

# Get all the file names in datadir and their integer counterparts
file_list, int_file_list, nfiles = get_file_lists(radatadir)

# Read in CLAs
args = sys.argv[2:]
nargs = len(args)

index_first, index_last = get_desired_range(int_file_list, args)

nrec = 1 # By default, only include 1 record per Shell_Avgs slice (otherwise
        # there will be a LOT of data!
for i in range(nargs):
    arg = args[i]
    if arg == '-nrec':
        nrec = int(args[i+1])

# Set the timetrace savename by the directory, what we are saving, and first and last
# iteration files for the trace
savename = dirname_stripped + '_trace_Shell_Avgs_' + file_list[index_first] + '_' +\
    file_list[index_last] + '.pkl'
savefile = datadir + savename    

# Initialize empty "vals" array for the timetrace
sh0 = Shell_Avgs(radatadir + file_list[index_first], '')

print ('Considering Shell_Avgs files %s through %s for the trace ...'\
       %(file_list[index_first], file_list[index_last]))

count = 0
iter1, iter2 = int_file_list[index_first], int_file_list[index_last]

vals = []
times = []
iters = []

for i in range(index_first, index_last + 1):
    print ('Adding Shell_Avgs/%s to the trace ...' %file_list[i])
    if i == index_first:
        sh = sh0
    else:   
        sh = Shell_Avgs(radatadir + file_list[i], '')

    nrec_tot = sh.niter
    nstep = nrec_tot // nrec
    if nstep == 0:
        nstep = 1
    for j in range(0, nrec_tot, nstep):
        vals.append(list(sh.vals[:, 0, :, j])) #  Each Shell_Avgs slice has shape (nr, nq)
                        # End "vals" array will have shape (ntimes, nr, nq)
        times.append(sh.time[j])
        iters.append(sh.iters[j])
        count += 1

vals = np.array(vals)
times = np.array(times)
iters = np.array(iters)

print ('Traced over %i Shell_Avgs slice(s) ...' %count)

# Save the avarage
print ('Saving file at ' + savefile + ' ...')
f = open(savefile, 'wb')
pickle.dump({'vals': vals, 'times': times, 'iters': iters, 'lut': sh0.lut, 'ntimes': count, 'iter1': iter1, 'iter2': iter2, 'rr': sh0.radius, 'nr': sh0.nr, 'qv': sh0.qv, 'nq': sh0.nq},\
        f, protocol=4)
f.close()
