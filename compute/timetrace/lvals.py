# Routine to trace Rayleigh G_Avgs data in time
# Created by: Loren Matilsky
# On: 03/07/2019
############################################################################
# This routine computes the trace in time/spherical harmonic degree l of
# quantities in the  AZ_Avgs data for a particular simulation (typically 
# these are the fluid variables v-vector, b-vector, S, and P).
# 
# By default, the routine traces over the last  niter = 100 files available,# though the user can specify a different range in sevaral ways:
# -n 10 (last 10 files)
# -range iter1 iter2 (no.s for start/stop data files; iter2 can be "last")
# -centerrange iter0 nfiles (trace about central file iter0 over nfiles)
#
# The final datacube output ('vals') will have shape
# (niter, ntheta, ndepths, nq)

# Import relevant modules
import numpy as np
import pickle
import sys, os
sys.path.append(os.environ['rasource'] + '/post_processing')
sys.path.append(os.environ['raco'])
from rayleigh_diagnostics import Shell_Spectra
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

radatadir = dirname + '/Shell_Spectra/'

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

# Set the timetrace savename by the directory, what we are saving, 
# and first and last iteration files for the trace
savename = dirname_stripped + '_trace_lvals_' +\
        file_list[index_first] + '_' + file_list[index_last] + '.pkl'
savefile = datadir + savename    

# Read in first Shell_Spectra file for (spectral) grid info
spec0 = Shell_Spectra(radatadir + file_list[index_first], '')

# Create 1d array of (float) l-values 
lmax = spec0.lmax
lvals = np.arange(lmax + 1, dtype='float')

# Trace over the relevant data range
print ('Considering Shell_Spectra files %s through %s for the trace ...'\
       %(file_list[index_first], file_list[index_last]))

iter1, iter2 = int_file_list[index_first], int_file_list[index_last]

iter1, iter2 = int_file_list[index_first], int_file_list[index_last]

vals = []
times = []
iters = []

for i in range(index_first, index_last + 1):
    print ('Adding Shell_Spectra/%s to the trace ...' %file_list[i])
    if i == index_first:
        spec = spec0
    else:   
        spec = Shell_Spectra(radatadir + file_list[i], '')

    local_ntimes = spec.niter
    for j in range(local_ntimes):
        vals_loc = spec.lpower[:, :, :, j, :] 
                # [ilval, ir, iq, iiter, i(tot,mean,conv)]
        vals.append(vals_loc.tolist())
        times.append(spec.time[j])
        iters.append(spec.iters[j])

# Convert lists into arrays
vals = np.array(vals)
times = np.array(times)
iters = np.array(iters)

# Get final shape of "vals" array
niter = len(iters)

print ('Traced over %i Shell_Spectra slice(s) ...' %niter)

# Save the trace
print ('Saving file at ' + savefile + ' ...')
f = open(savefile, 'wb')
pickle.dump({'vals': vals, 'times': times, 'iters': iters, \
'iter1': iter1, 'iter2': iter2, 'niter': niter, 'lut': spec0.lut, \
'qv': spec0.qv, 'nq': spec0.nq, 'rvals': spec0.radius, 'rinds': spec0.inds,\
'nr': spec0.nr, 'lvals': lvals, 'nell': spec0.nell, 'lmax': lmax},\
    f, protocol=4)
f.close()
