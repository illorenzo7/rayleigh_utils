# Routine to trace Rayleigh G_Avgs data in time
# Created by: Loren Matilsky
# On: 02/26/2020
############################################################################
# This routine computes the trace in time/spherical harmonics l and m
# of the absolute-squared power associated with  v_r, v_theta, and v_phi
# 
# By default, the routine traces over the last  niter = 100 files available,# though the user can specify a different range in sevaral ways:
# -n 10 (last 10 files)
# -range iter1 iter2 (no.s for start/stop data files; iter2 can be "last")
# -centerrange iter0 nfiles (trace about central file iter0 over nfiles)
# -all (all the files in the data directory)
#
# The final datacube output ('vals') will have shape
# (ntimes, nell, ndepths, 3)

# Import relevant modules
import numpy as np
import pickle
import sys, os
sys.path.append(os.environ['rapp'])
sys.path.append(os.environ['raco'])
from rayleigh_diagnostics import Shell_Spectra
from common import *
from get_parameter import get_parameter

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

if nargs == 0:
    index_first, index_last = nfiles - 101, nfiles - 1  
    # By default trace over the last 100 files
else:
    index_first, index_last = get_desired_range(int_file_list, args)

# Read in first Shell_Spectra file for (spectral) grid info
spec0 = Shell_Spectra(radatadir + file_list[index_first], '')

# Create 1d array of (float) l-values 
lmax = spec0.lmax
mmax = spec0.mmax
nell = spec0.nell
nm = spec0.nm
lvals = np.arange(nell, dtype='float')
mvals = np.arange(nm, dtype='float')

# Set the timetrace savename by the directory, what we are saving, 
# and first and last iteration files for the trace
savename = dirname_stripped + '_trace_lmvals_v_' +\
        file_list[index_first] + '_' + file_list[index_last] + '.pkl'
savefile = datadir + savename    

# Trace over the relevant data range
print ('Considering Shell_Spectra files %s through %s for the trace ...'\
       %(file_list[index_first], file_list[index_last]))

iter1, iter2 = int_file_list[index_first], int_file_list[index_last]

vals_vr = []
vals_vt = []
vals_vp = []
times = []
iters = []

for i in range(index_first, index_last + 1):
    print ('Adding Shell_Spectra/%s to the trace ...' %file_list[i])
    if i == index_first:
        spec = spec0
    else:   
        spec = Shell_Spectra(radatadir + file_list[i], '')

    local_ntimes = spec.niter
    lut = spec.lut
    for j in range(local_ntimes):
        vals_vr_loc = spec.vals[:, :, :, lut[1], j]
        vals_vt_loc = spec.vals[:, :, :, lut[2], j]
        vals_vp_loc = spec.vals[:, :, :, lut[3], j]
        vals_vr.append(vals_vr_loc.tolist())
        vals_vt.append(vals_vt_loc.tolist())
        vals_vp.append(vals_vp_loc.tolist())
        times.append(spec.time[j])
        iters.append(spec.iters[j])

# Convert lists into arrays
vals_vr = np.array(vals_vr)
vals_vt = np.array(vals_vt)
vals_vp = np.array(vals_vp)
times = np.array(times)
iters = np.array(iters)

# Get final shape of "vals" array
niter = len(iters)

print ('Traced over %i Shell_Spectra slice(s) ...' %niter)

# Save the trace
print ('Saving file at ' + savefile + ' ...')
f = open(savefile, 'wb')
pickle.dump({'vals_vr': vals_vr, 'vals_vt': vals_vt, 'vals_vp': vals_vp, \
        'times': times, 'iters': iters, \
'iter1': iter1, 'iter2': iter2, 'niter': niter, 'lut': spec0.lut, \
'qv': spec0.qv, 'nq': spec0.nq, 'rvals': spec0.radius, 'rinds': spec0.inds,\
'nr': spec0.nr, 'lvals': lvals, 'mvals': mvals, 'nell': spec0.nell,\
'lmax': lmax},\
    f, protocol=4)
f.close()
