import numpy as np
import sys, os
sys.path.append('/altair/loma3853/rayleigh/compute')
from diagnostic_reading import AzAverage, ReferenceState
from common import strip_dirname, get_widest_range_file, get_iters_from_file
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

# Get grid info (if it's not computed already using grid_info.py, this will fail)
rr,tt,cost,sint,rr_depth,ri,ro,d = np.load(datadir + 'grid_info.npy')
nr = len(rr)
nt = len(tt)

# Read in pressure and entropy
pressure_file = get_widest_range_file(datadir, 'p_spherical_mean')
pressure = np.load(datadir + pressure_file)

entropy_file = get_widest_range_file(datadir, 's_spherical_mean')
entropy = np.load(datadir + entropy_file)

# Get reference values
cp = get_parameter(dirname, 'pressure_specific_heat')
poly_n = get_parameter(dirname, 'poly_n')
gamma = 1. + 1./poly_n
ref = ReferenceState(dirname + '/reference', '')
p_ref = (ref.pressure)
t_ref = (ref.temperature)

temperature_ratio = pressure/p_ref*(1. - 1./gamma) + entropy/cp
temperature = temperature_ratio*t_ref

# Set the savename by the directory, what we are saving, and first and last
# iteration files for the average
iter1, iter2 = get_iters_from_file(pressure_file)
savename = dirname_stripped + '_t_spherical_mean_' + str(iter1).zfill(8) + '_' + str(iter2).zfill(8) + '.npy'
savefile = datadir + savename    

print ('Saving spherical temperature mean at ' + savefile + ' ...')
np.save(savefile, temperature)
