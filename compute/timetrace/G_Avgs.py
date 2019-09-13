# Routine to trace Rayleigh G_Avgs data in time
# Created by: Loren Matilsky
# On: 12/18/2018
##################################################################
# This routine computes the trace in time of the values in the G_Avgs data 
# for a particular simulation. 

# By default, the routine traces over the last 100 files of datadir, though
# the user can specify a different range in sevaral ways:
# -n 10 (last 10 files)
# -range iter1 iter2 (names of start and stop data files; iter2 can be "last")
# -centerrange iter0 nfiles (trace about central file iter0 over nfiles)

# Import relevant modules
import numpy as np
import pickle
import sys, os
sys.path.append(os.environ['rapp'])
sys.path.append(os.environ['raco'])
from rayleigh_diagnostics import G_Avgs
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

radatadir = dirname + '/G_Avgs/'

# Get all the file names in datadir and their integer counterparts
file_list, int_file_list, nfiles = get_file_lists(radatadir)

# Read in CLAs
args = sys.argv[2:]
nargs = len(args)

if (nargs == 0):
    index_first, index_last = nfiles - 101, nfiles - 1  
    # By default trace over the last 100 files
else:
    index_first, index_last = get_desired_range(int_file_list, args)

# Set the timetrace savename by the directory, what we are saving, and first and last
# iteration files for the trace
savename = dirname_stripped + '_trace_G_Avgs_' + file_list[index_first] + '_' +\
    file_list[index_last] + '.pkl'
savefile = datadir + savename    

# Initialize empty "vals" array for the timetrace
g0 = G_Avgs(radatadir + file_list[index_first], '')
vals = np.zeros_like(g0.vals[0, :]) # G_Avgs has just numbers: one for each quantity

print ('Considering G_Avgs files %s through %s for the trace ...'\
       %(file_list[index_first], file_list[index_last]))

count = 0
iter1, iter2 = int_file_list[index_first], int_file_list[index_last]

vals = []
times = []
iters = []

for i in range(index_first, index_last + 1):
    print ('Adding G_Avgs/%s to the trace ...' %file_list[i])
    if i == index_first:
        g = g0
    else:   
        g = G_Avgs(radatadir + file_list[i], '')

    local_ntimes = g.niter
    for j in range(local_ntimes):
        vals.append(list(g.vals[j, :]))
        times.append(g.time[j])
        iters.append(g.iters[j])
        count += 1

vals = np.array(vals)
vals = np.transpose(vals)
times = np.array(times)
iters = np.array(iters)

print ('Traced over %i G_Avgs slice(s) ...' %count)

# Save the avarage
print ('Saving file at ' + savefile + ' ...')
f = open(savefile, 'wb')
pickle.dump({'vals': vals, 'times': times, 'iters': iters, 'lut': g0.lut, 'count': count, 'iter1': iter1, 'iter2': iter2, 'qv': g0.qv, 'nq': g0.nq},\
        f, protocol=4)
f.close()
