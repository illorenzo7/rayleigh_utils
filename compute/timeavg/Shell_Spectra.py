# Routine to average Rayleigh Shell_Spectra data in time
# For the full power (both l, m) the "power" starts complex, so the 
# square-average is taken before averaging. 
# For the lpower, the square has already been taken in the Rayleigh
# diagnostics.
# Created by: Loren Matilsky
# On: 12/07/2018
##################################################################
# This routine computes the average in time of the values in the Shell_Spectra data 
# for a particular simulation. 

# By default, the routine averages over the last 100 files of the data directory, though
# the user can specify a different range in sevaral ways:
# -n 10 (last 10 files)
# -range iter1 iter2 (names of start and stop data files; iter2 can be "last")
# -centerrange iter0 nfiles (average about central file iter0 over nfiles)

# Import relevant modules
import numpy as np
import pickle
import sys, os
sys.path.append(os.environ['rapp'])
from rayleigh_diagnostics import Shell_Spectra
sys.path.append(os.environ['co'])
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

radatadir = dirname + '/Shell_Spectra/'

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

# Get grid_info for radius stuff (not contained in Shell_Spectra objects
# If no grid_info.npy, prompt user to run grid_info.py
try: 
    dummy, dummy, dummy, dummy, dummy, ri, ro, d = np.load(datadir + 'grid_info.npy')
except:
    print('The compute/timeavg/Shell_Spectra routine needs pre-computed file grid_info.npy. Run compute/grid_info.py and try again. Exiting...')
    sys.exit(1)

# Set the timeavg savename by the directory, what we are saving, and first and last
# iteration files for the average
savename = dirname_stripped + '_Shell_Spectra_' + file_list[index_first] + '_' +\
    file_list[index_last] + '.pkl'
savefile = datadir + savename    

# Initialize empty "vals" array for the time average
spec0 = Shell_Spectra(radatadir + file_list[index_first], '')
fullpower = np.zeros_like(spec0.vals[:, :, :, :, 0], dtype='float') 
lpower = np.zeros_like(spec0.lpower[:, :, :, 0, :])

# Average over the relevant data range, summing everything and then dividing
#   by the number of "slices" added at the end
print ('Considering Shell_Spectra files %s through %s for the average ...'\
       %(file_list[index_first], file_list[index_last]))

count = 0
iter1, iter2 = int_file_list[index_first], int_file_list[index_last]

for i in range(index_first, index_last + 1):
    print ('Adding Shell_Spectra/%s to the average ...' %file_list[i])
    if i == index_first:
        spec = spec0
    else:   
        spec = Shell_Spectra(radatadir + file_list[i], '')

    local_ntimes = spec.niter
    for j in range(local_ntimes):
        fullpower += np.abs(spec.vals[:, :, :, :, j])**2
        lpower += spec.lpower[:, :, :, j, :]
        count += 1

fullpower /= count
lpower /= count
print ('Averaged over %i Shell_Spectra ...' %count)

# Create 2d arrays corresponding to l, m values
lmax = spec0.lmax
mmax = spec0.mmax

lvals = np.arange(lmax + 1, dtype='float')
mvals = np.arange(mmax + 1, dtype='float')

lvals_2d, mvals_2d = np.meshgrid(lvals, mvals, indexing='ij')

# Compute some quantities having to do with radius
rvals = spec0.radius
rvals_depth = (ro - rvals)/d
rvals_height = (rvals - ri)/d

# Save the avarage
print ('Saving file at ' + savefile + ' ...')
f = open(savefile, 'wb')
pickle.dump({'fullpower': fullpower, 'lpower': lpower, 'lut': spec0.lut, 'count': count, 'iter1': iter1, 'iter2': iter2, 'qv': spec0.qv, 'nq': spec0.nq, 'rinds': spec0.inds, 'rvals': rvals, 'rvals_depth': rvals_depth, 'rvals_height': rvals_height, 'nr': spec0.nr, 'ri': ri, 'ro': ro, 'd': d, 'lvals': lvals, 'lvals_2d': lvals_2d, 'nell': spec0.nell, 'lmax': lmax, 'mvals': mvals, 'mvals_2d': mvals_2d, 'nm': spec0.nm, 'mmax': mmax}, f, protocol=4)
f.close()
