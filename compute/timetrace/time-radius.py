# Routine to trace Rayleigh G_Avgs data in time
# Created by: Loren Matilsky
# On: 02/27/2019
############################################################################
# This routine computes the trace in time/radius of quantities in the 
# AZ_Avgs data for a particular simulation. 
# 
# By default, the quantities traced are the nq (=5, 8) fluid variables 
# (vr, vt, vp, s, p, [bp, bt, and bp] (if magnetism present).
# This may be changed via the '-vars' CLA, e.g., -vars '1 2 3 1425 1427'.
#
# By default, the 8 variables are computed at each time (for all radii)
# at nlats = 13 latitudes equally spaced from -90 degrees to 90 degrees,
# in increments of 15 degrees.
# This may be changed via the '-lats' CLA, e.g., -lats '0 5 25'
#
# By default, the routine traces over the last  niter = 100 files available,# though the user can specify a different range in sevaral ways:
# -n 10 (last 10 files)
# -range iter1 iter2 (no.s for start/stop data files; iter2 can be "last")
# -centerrange iter0 nfiles (trace about central file iter0 over nfiles)
#
# The final datacube output ('vals') will have shape
# (niter, nlats, nr, nq)

# Import relevant modules
import numpy as np
import pickle
import sys, os
sys.path.append(os.environ['rasource'] + '/post_processing')
sys.path.append(os.environ['raco'])
from rayleigh_diagnostics import AZ_Avgs
from common import get_file_lists, get_desired_range, strip_dirname
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

radatadir = dirname + '/AZ_Avgs/'

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

# Set other defaults
qvals = [1, 2, 3, 301, 302]
magnetism = get_parameter(dirname, 'magnetism')
if magnetism:
    qvals.append(801)
    qvals.append(802)
    qvals.append(803)
lats = [-85., -75., -60., -45., -30., -15., 0., 15., 30., 45., 60., 75.,\
        85.]
user_specified_vals = False
user_specified_lats = False
tag_already_provided = False

for i in range(nargs):
    arg = args[i]
    if arg == '-vars':
        qvals = []
        qvals_str = args[i+1].split() 
        for qval_str in qvals_str:
            qvals.append(int(qval_str))
        user_specified_vals = True
    elif arg == '-lats':
        lats = []
        lats_str = args[i+1].split()
        for lat_str in lats_str:
            lats.append(float(lat_str))
            user_specified_lats = True
    elif arg == '-tag':
        tag_already_provided = True
        tag = args[i+1] + '_'

# If the user chose different depths/vals, they must put a unique tag
# on the data
if not tag_already_provided:
    tag = ''
    if user_specified_lats or user_specified_vals:
        tag = input("You chose latitudes/quantities different from the\n\
    default; please enter a unique tag for the output data file:\n")
        tag = tag + '_'

# Set the timetrace savename by the directory, what we are saving, 
# and first and last iteration files for the trace (and optional tag)
savename = dirname_stripped + '_time-radius_' + tag +\
        file_list[index_first] + '_' + file_list[index_last] + '.pkl'
savefile = datadir + savename    
print('Your data will be saved in the file %s' %savename)

# Read in first AZ_Avgs file for grid info
az0 = AZ_Avgs(radatadir + file_list[index_first], '')

# Get quantity indices associated with qvals
qinds = az0.lut[qvals]

# Get bunch of grid info
rr = az0.radius
ri, ro = np.min(rr), np.max(rr)
d = ro - ri
rr_depth = (ro - rr)/d
rr_height = (rr - ri)/d
sint = az0.sintheta
cost = az0.costheta
tt = np.arccos(cost)
tt_lat = (np.pi/2 - tt)*180/np.pi
nr = az0.nr
nt = az0.ntheta

# compute some derivative quantities for the grid
tt_2d, rr_2d = np.meshgrid(tt, rr, indexing='ij')
sint_2d = np.sin(tt_2d); cost_2d = np.cos(tt_2d)
xx = rr_2d*sint_2d
zz = rr_2d*cost_2d

# get theta-indices associated with depths
theta_inds = []
for lat in lats:
    theta_inds.append(np.argmin(np.abs(tt_lat - lat)))

print ('Considering AZ_Avgs files %s through %s for the trace ...'\
       %(file_list[index_first], file_list[index_last]))

iter1, iter2 = int_file_list[index_first], int_file_list[index_last]

vals = []
times = []
iters = []

for i in range(index_first, index_last + 1):
    print ('Adding AZ_Avgs/%s to the trace ...' %file_list[i])
    if i == index_first:
        az = az0
    else:   
        az = AZ_Avgs(radatadir + file_list[i], '')

    local_ntimes = az.niter
    for j in range(local_ntimes):
        vals_reduced = az.vals[:, :, qinds, j]
        vals_loc = vals_reduced[theta_inds, :, :]
        vals.append(vals_loc.tolist())
        times.append(az.time[j])
        iters.append(az.iters[j])

# Convert lists into arrays
vals = np.array(vals)
times = np.array(times)
iters = np.array(iters)

# Get final shape of "vals" array
niter = len(iters)
ntheta = az.ntheta
nlats = len(lats)
nq = len(qvals)

print ('Traced over %i AZ_Avgs slice(s) ...' %niter)

# Save the avarage
print ('Saving file at ' + savefile + ' ...')
f = open(savefile, 'wb')
pickle.dump(\
{'vals': vals, 'times': times, 'iters': iters, 'lats': lats,\
'qvals': qvals, 'theta_inds': theta_inds, 'qinds': qinds, 'niter': niter,\
'nlats': nlats, 'nq': nq, 'iter1': iter1, 'iter2': iter2, 'rr': rr,\
'rr_depth': rr_depth, 'rr_height': rr_height, 'nr': nr, 'ri': ri,\
'ro': ro, 'd': d, 'tt': tt, 'tt_lat': tt_lat, 'sint': sint, 'cost': cost,\
'ntheta': ntheta, 'rr_2d': rr_2d, 'tt_2d': tt_2d, 'sint_2d': sint_2d,\
'cost_2d': cost_2d, 'xx': xx, 'zz': zz}, f, protocol=4)
f.close()
