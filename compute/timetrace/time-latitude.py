# Routine to trace Rayleigh G_Avgs data in time
# Created by: Loren Matilsky
# On: 02/27/2019
############################################################################
# This routine computes the trace in time/latitude of quantities in the 
# AZ_Avgs data for a particular simulation. 
# 
# By default, the quantities traced are the nq (=5, or 8) fluid variables 
# (vr, vt, vp, s, p, [bp, bt, and bp] (if magnetism present).
# This may be changed via the '-vars' CLA, e.g., -vars '1 2 3 1425 1427'.
#
# By default, the 8 variables are computed at each time (for all latitudes)
# at ndepths = 9 depths equally spaced throughout the shell (like the shell
# slice levels).
# This may be changed via the '-depths' CLA, e.g., -depths '0 0.4 0.75 0.95'
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
sys.path.append(os.environ['co'])
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
depths = [0.05, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 0.95]
user_specified_vals = False
user_specified_depths = False

for i in range(nargs):
    arg = args[i]
    if arg == '-vars':
        qvals = []
        qvals_str = args[i+1].split() 
        for qval_str in qvals_str:
            qvals.append(int(qval_str))
        user_specified_vals = True
    elif arg == '-depths':
        depths = []
        depths_str = args[i+1].split()
        for depth_str in depths_str:
            depths.append(float(depth_str))
            user_specified_depths = True

# Set the timetrace savename by the directory, what we are saving, 
# and first and last iteration files for the trace
basename = dirname_stripped + '_time-latitude_' +\
        file_list[index_first] + '_' + file_list[index_last]

# If the user chose different depths/vals, they must put a unique tag
# on the data
more = ''
if user_specified_depths or user_specified_vals:
    more = input("You chose latitudes/quantities different from the\n\
default; please enter a unique tag for the output data file:\n")
    more = '_' + more
    print('Your data will be saved in the file %s'\
            %(basename + more + '.pkl'))

savename = basename + more + '.pkl'
savefile = datadir + savename    

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

# get r-indices associated with depths
rinds = []
for depth in depths:
    rinds.append(np.argmin(np.abs(rr_depth - depth)))

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
        vals_loc = vals_reduced[:, rinds, :]
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
ndepths = len(depths)
nq = len(qvals)

print ('Traced over %i AZ_Avgs slice(s) ...' %niter)

# Save the avarage
print ('Saving file at ' + savefile + ' ...')
f = open(savefile, 'wb')
pickle.dump({'vals': vals, 'times': times, 'iters': iters,\
'depths': depths,'qvals': qvals, 'rinds': rinds, 'qinds': qinds,\
'niter': niter,\
'ndepths': ndepths, 'nq': nq, 'iter1': iter1, 'iter2': iter2, 'rr': rr,\
'rr_depth': rr_depth, 'rr_height': rr_height, 'nr': nr, 'ri': ri,\
'ro': ro, 'd': d, 'tt': tt, 'tt_lat': tt_lat, 'sint': sint, 'cost': cost,\
'ntheta': ntheta, 'rr_2d': rr_2d, 'tt_2d': tt_2d, 'sint_2d': sint_2d,\
'cost_2d': cost_2d, 'xx': xx, 'zz': zz}, f, protocol=4)
f.close()
