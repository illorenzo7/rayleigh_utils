# Routine to trace Rayleigh G_Avgs data in time
# Created by: Loren Matilsky
# On: 08/15/2019
############################################################################
# This routine computes the trace in time/longitude of quantities in the 
# Shell_Slices data for a particular simulation. 
# 
# By default, the quantities traced are the nq (=5, or 8) fluid variables 
# (vr, vt, vp, s, p, [bp, bt, and bp] (if magnetism present).
# This may be changed via the '-qvals' CLA, e.g., -qvals '1 2 3 1425 1427'.
#
# By default, the 8 variables are computed at each time, in a latitude strip defined by
# (clat, dlat) = ([central latitude for average], [range of latitudes to average over])
# at ndepths = 9 depths equally spaced throughout the shell (like the shell
# slice levels).
# This may be changed via the '-depths' CLA, e.g., -depths '0 0.4 0.75 0.95'
#
# By default, the routine traces over all Shell_Slices, though user can specify an 
# alternate range, by, e.g.,
# -n 10 (last 10 files)
# -range iter1 iter2 (no.s for start/stop data files; iter2 can be "last")
# -centerrange iter0 nfiles (trace about central file iter0 over nfiles)
#
# The final datacube output ('vals') will have shape
# (nphi, niter, ndepths, nq)

# Import relevant modules
import numpy as np
import pickle
import sys, os
sys.path.append(os.environ['rasource'] + '/post_processing')
sys.path.append(os.environ['co'])
from rayleigh_diagnostics import AZ_Avgs, Shell_Slices
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

# We are interested in latitude strips from the shell slices
radatadir = dirname + '/Shell_Slices/'

# Get all the file names in datadir and their integer counterparts
file_list, int_file_list, nfiles = get_file_lists(radatadir)

# Read in CLAs
args = sys.argv[2:]
nargs = len(args)

# By default, have the iter indices range over full file range
index_first, index_last = 0, nfiles - 1
newrange = False
for arg in args:
    if arg in ['-range', '-centerrange', '-leftrange', '-rightrange', '-n', '-f', '-all', '-iter']:
        newrange = True
if newrange:
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
tag = ''
clat = 10.
dlat = 20. # by default average over the first 20 degrees of the Northern hemisphere
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
    elif arg == '-tag':
        tag = args[i+1]
    elif arg == '-clat':
        clat = float(args[i+1])
    elif arg == '-dlat':
        dlat = float(args[i+1])
        
# Set the timetrace savename by the directory, what we are saving, 
# and first and last iteration files for the trace (and optional tag)
if clat > 0.:
    hemisphere = 'N'
else:
    hemisphere = 'S'
    
savename = dirname_stripped + ('_time-longitude_clat%s%02.0f_dlat%02.0f' %(hemisphere, clat, dlat)) +\
        tag + file_list[index_first] + '_' + file_list[index_last] + '.pkl'
savefile = datadir + savename    
print('Your data will be saved in the file %s' %savename)

# Read in first Shell_Slices file 
a0 = Shell_Slices(radatadir + file_list[index_first], '')

# Also read in time-latitude file (needed for differential rotation)
tl_file = get_widest_range_file(datadir, 'time-latitude')
print("Getting time-latitude trace data from %s ..." %tl_file)
di_tl = get_dict(datadir + tl_file)
vals_tl = di_tl['vals']
qvals_tl = di_tl['qvals']
times_tl = di_tl['times']

# Read in grid info from time-latitude trace data
rr, ri, ro = di_tl['rr'], di_tl['ri'], di_tl['ro']
sint, tt_lat = di_tl['sint'], di_tl['tt_lat']
nt, nr = len(sint), len(rr)
sint_2d = sint.reshape((1, nt))
nphi = 2*nt
lons = np.arange(0., 360., 360/nphi)

# Desired latitude range
ith1 = np.argmin(np.abs(tt_lat - (clat - dlat/2.)))
ith2 = np.argmin(np.abs(tt_lat - (clat + dlat/2.)))
lats_strip = tt_lat[ith1:ith2+1]
# Get quantity indices associated with qvals
qinds = a0.lut[qvals]





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
