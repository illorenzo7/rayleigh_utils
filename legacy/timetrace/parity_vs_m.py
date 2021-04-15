# Routine to trace Rayleigh G_Avgs data in time
# Created by: Loren Matilsky
# On: 09/29/2019
############################################################################
# This routine computes the trace in time/spherical harmonic degree m of
# the power in even/odd modes (for even: l even if m is even, l odd if m is odd;
# for odd: l odd if m is even, l even if m is odd) for quantities in the 
# Shell_Spectra data for a particular simulation (typically these are the fluid # variables v-vector, b-vector, S, and P).
# 

# The final datacube output ('vals') will have shape
# (ntimes, nm, ndepths, nq, 2)

# Import relevant modules
import numpy as np
import pickle
import sys, os
sys.path.append(os.environ['rapp'])
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
if not os.path.isdir(datadir):
    os.makedirs(datadir)

radatadir = dirname + '/Shell_Spectra/'

# Get all the file names in datadir and their integer counterparts
file_list, int_file_list, nfiles = get_file_lists(radatadir)

# Read in CLAs

args = sys.argv[2:]
nargs = len(args)

if nargs == 0:
    index_first, index_last = 0, nfiles - 1  
    # By default trace over all files
else:
    index_first, index_last = get_desired_range(int_file_list, args)

# Set the timetrace savename by the directory, what we are saving, 
# and first and last iteration files for the trace
savename = 'parity_vs_m-' +\
        file_list[index_first] + '_' + file_list[index_last] + '.pkl'
savefile = datadir + savename    

# Read in first Shell_Spectra file for (spectral) grid info
spec0 = Shell_Spectra(radatadir + file_list[index_first], '')

# Create 1d array of (float) l-values 
lmax = spec0.lmax
mmax = spec0.mmax
nm = spec0.nm
lvals = np.arange(lmax + 1, dtype='float')
mvals = np.arange(mmax + 1, dtype='float')
qv = spec0.qv
nq = spec0.nq
rvals = spec0.radius
nr = spec0.nr

# Trace over the relevant data range
print ('Considering Shell_Spectra files %s through %s for the trace ...'\
       %(file_list[index_first], file_list[index_last]))

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
        vals_loc = spec.vals[:, :, :, :, j]
        parity_loc = np.zeros((nm, nr, nq, 2)) 
        # last index: 0 for even, 1 for odd
        # conveniently, im = mval, il = lval
        for im in range(nm):
            if im % 2 == 0: # m is even: for even parity (index = 0), sum over
                            # even l; for odd parity (index = 1), sum over odd l
                parity_loc[im, :, :, 0] =\
                        np.sum(np.abs(vals_loc[::2, im, :, :])**2., axis=0)
                parity_loc[im, :, :, 1] =\
                        np.sum(np.abs(vals_loc[1::2, im, :, :])**2., axis=0)
            elif im % 2 == 1: # m is odd: for even parity (index = 0), sum over
                            # odd l; for odd parity (index = 1), sum over even l
                parity_loc[im, :, :, 0] =\
                        np.sum(np.abs(vals_loc[1::2, im, :, :])**2., axis=0)
                parity_loc[im, :, :, 1] =\
                        np.sum(np.abs(vals_loc[::2, im, :, :])**2., axis=0)

        vals.append(parity_loc.tolist())
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
'qv': qv, 'nq': spec0.nq, 'rvals': rvals, 'rinds': spec0.inds,\
'nr': spec0.nr, 'lvals': lvals, 'nell': spec0.nell, 'lmax': lmax,\
'mvals': mvals, 'nm': nm, 'mmax': spec0.mmax},\
    f, protocol=4)
f.close()
