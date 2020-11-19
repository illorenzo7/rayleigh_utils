# Routine to trace Rayleigh Shell_Avgs data in time
# Created by: Loren Matilsky
# On: 03/31/2020
##################################################################
# This routine computes the trace in time of the values in the Shell_Avgs data 
# for a particular simulation, averaging separately over the RZ and CZ
# (So, essentially G_Avgs data over each domain)
# Also computes the internal energy (rho*T*S), the full one, and one with
# the top value of the entropy subtracted, since the dynamics should
# be insensitive to absolute values of the entropy

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
from rayleigh_diagnostics import Shell_Avgs, GridInfo
from common import get_file_lists, get_desired_range, strip_dirname,\
        get_parameter
from get_eq import get_eq

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

# Get grid information
gi = GridInfo(dirname + '/grid_info', '')
ir_bcz = get_parameter(dirname, 'ncheby')[1] - 1

# Get averaging weights for CZ and RZ separately
rw = gi.rweights
nr = gi.nr
rw_cz = rw[:ir_bcz + 1]
nr_cz = len(rw_cz)
rw_rz = rw[ir_bcz + 1:]
nr_rz = len(rw_rz)
rw_cz /= np.sum(rw_cz)
rw_rz /= np.sum(rw_rz)
rw = rw.reshape((nr, 1))
rw_cz = rw_cz.reshape((nr_cz, 1))
rw_rz = rw_rz.reshape((nr_rz, 1))

# Get rho*T
eq = get_eq(dirname)
rhot = (eq.density*eq.temperature).reshape((1, nr))

# Read in CLAs
args = sys.argv[2:]
nargs = len(args)

if nargs == 0:
    index_first, index_last = 0, nfiles - 1  
    # By default trace over all files
else:
    index_first, index_last = get_desired_range(int_file_list, args)

# Set the timetrace savename by the directory, what we are saving, and first and last
# iteration files for the trace
savename = dirname_stripped + '_trace_2dom_G_Avgs_' + file_list[index_first] + '_' +\
    file_list[index_last] + '.pkl'
savefile = datadir + savename    

# Initialize empty "vals" array for the timetrace
sh0 = Shell_Avgs(radatadir + file_list[index_first], '')

print ('Considering Shell_Avgs files %s through %s for the trace ...'\
       %(file_list[index_first], file_list[index_last]))

count = 0
iter1, iter2 = int_file_list[index_first], int_file_list[index_last]

vals = []
vals_cz = []
vals_rz = []
times = []
iters = []

for i in range(index_first, index_last + 1):
    print ('Adding Shell_Avgs/%s to the trace ...' %file_list[i])
    if i == index_first:
        sh = sh0
    else:   
        sh = Shell_Avgs(radatadir + file_list[i], '')

    #local_ntimes = sh.niter
    local_ntimes = sh.niter
    for j in range(local_ntimes):
        vals_loc = sh.vals[:, 0, :, j]

        # add in internal energy
        inte_loc = rhot*sh.vals[:, 0, sh.lut[501], j]
        # top S subtracted
        inte_loc_subt = rhot*(sh.vals[:, 0, sh.lut[501], j] -\
                sh.vals[0, 0, sh.lut[501], j])
        # bottom S subtracted
        inte_loc_subb = rhot*(sh.vals[:, 0, sh.lut[501], j] -\
                sh.vals[-1, 0, sh.lut[501], j])

        # add in the three energies
        vals_loc = np.hstack((vals_loc, inte_loc.T, inte_loc_subt.T,\
                inte_loc_subb.T))

        # Get the values in the CZ/RZ separately
        vals_cz_loc = vals_loc[:ir_bcz + 1]
        vals_rz_loc = vals_loc[ir_bcz + 1:]

        gav = np.sum(rw*vals_loc, axis=0)
        gav_cz = np.sum(rw_cz*vals_cz_loc, axis=0)
        gav_rz = np.sum(rw_rz*vals_rz_loc, axis=0)

        vals.append(list(gav)) 
        vals_cz.append(list(gav_cz)) 
        vals_rz.append(list(gav_rz)) 
        times.append(sh.time[j])
        iters.append(sh.iters[j])
        count += 1

vals = np.array(vals)
vals_cz = np.array(vals_cz)
vals_rz = np.array(vals_rz)
times = np.array(times)
iters = np.array(iters)

# Also append the lut (making the inte, inte_subt, and inte_subb quantities
# (4000, 4001, 4002)
lut_app = np.array([sh0.nq, sh0.nq + 1, sh0.nq + 2])
lut = np.hstack((sh0.lut, lut_app))

print ('Traced over %i Shell_Avgs slice(s)' %count)

# Save the avarage
print ('Saving file at ' + savefile)
f = open(savefile, 'wb')
pickle.dump({'vals': vals, 'vals_cz': vals_cz, 'vals_rz': vals_rz, 'times': times, 'iters': iters, 'lut': lut, 'ntimes': count, 'iter1': iter1, 'iter2': iter2, 'rr': sh0.radius, 'nr': sh0.nr, 'qv': sh0.qv, 'nq': sh0.nq},\
        f, protocol=4)
f.close()
