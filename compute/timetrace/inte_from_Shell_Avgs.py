# Routine to trace Rayleigh Shell_Avgs data in time
# Created by: Loren Matilsky
# On: 11/08/2019
##################################################################
# This routine computes the trace in time of inte directly by 
# averaging the Shell_Avgs data in radius--currently there is a bug
# in quantity codes 701 (thermal_energy_full) and 707 (thermal_energy_sq)
#
# Assumes user has already created a trace_Shell_Avgs file for the sim

# Import relevant modules
import numpy as np
import pickle
import sys, os
sys.path.append(os.environ['rapp'])
sys.path.append(os.environ['raco'])
from rayleigh_diagnostics import GridInfo, ReferenceState
from reference_tools import equation_coefficients
from common import strip_dirname, get_widest_range_file, get_dict

# Get the name of the run directory
dirname = sys.argv[1]
# Get the stripped name to use in file naming
dirname_stripped = strip_dirname(dirname)

# Find the relevant place to store the data
datadir = dirname + '/data/'

# Get the Shell_Averaged entropy as a function of time
the_file = get_widest_range_file(datadir, 'trace_Shell_Avgs')
print("Getting trace of Shell_Avgs data from %s ..." %the_file)
di = get_dict(datadir + the_file)
iter1 = di['iter1']
iter2 = di['iter2']
vals = di['vals']
lut = di['lut']
entropy = vals[:, :, lut[501]]

# Get density * temperature
nr = di['nr']
try:
    ref = ReferenceState(dirname + '/reference')
    rhot = (ref.density*ref.temperature).reshape((1, nr))
except:
    eq = equation_coefficients()
    eq.read(dirname + '/equation_coefficients')
    rhot = (eq.functions[0]*eq.functions[3]).reshape((1, nr))

# Get the radial integration weights
gi = GridInfo(dirname + '/grid_info')
rw = gi.rweights.reshape((1, nr))

# Calculate the internal energy
inte = np.sum(rhot*entropy*rw, axis=1)

# Set the timetrace savename by the directory, what we are saving, and first and last
# iteration files for the trace
savename = dirname_stripped + '_inte_from_Shell_Avgs_' +\
        str(iter1).zfill(8) + '_' + str(iter2).zfill(8) + '.pkl'
savefile = datadir + savename    

# Save the thermal energy
print ('Saving file at ' + savefile + ' ...')
f = open(savefile, 'wb')
pickle.dump({'inte': inte, 'times': di['times'], 'iters': di['iters'], 'ntimes': di['ntimes'], 'iter1': iter1, 'iter2': iter2}, f, protocol=4)
f.close()
