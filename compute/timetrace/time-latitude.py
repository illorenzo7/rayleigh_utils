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
sys.path.append(os.environ['rapp'])
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

if nargs == 0:
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
tag_already_provided = False
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
        tag_already_provided = True
        tag = args[i+1] + '_'
    elif arg == '-rzquarter': # 9 depths in RZ and CZ, with RZ depth
        # 0.25 of CZ depth
        print("Taking 9 depths in CZ and RZ each")
        print("assuming depth RZ = (1/4) depth CZ")
        depths = [0.04, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.76,\
               0.81, 0.825, 0.85, 0.875, 0.9, 0.925, 0.95, 0.975, 0.990]
    elif arg == '-rzhalf': # 9 depths in RZ and CZ, with RZ depth
        # 0.5 of CZ depth
        print("Taking 9 depths in CZ and RZ each")
        print("assuming depth RZ = (1/2) depth CZ")
        depths =[0.03333333, 0.08333333, 0.16666667, 0.25, 0.33333333,\
            0.41666667, 0.5, 0.58333333, 0.63333333, 0.68333333,\
            0.70833333, 0.75, 0.79166667, 0.83333333, 0.875,\
            0.91666667, 0.95833333, 0.9833333]
    elif arg == '-rz75': # 9 depths in RZ and CZ, with RZ depth
        # 0.75 of CZ depth
        print("Taking 9 depths in CZ and RZ each")
        print("assuming depth RZ = (3/4) depth CZ")
        depths = 1.0 - np.array([0.02142857, 0.05357143, 0.10714286,\
                0.16071429, 0.21428571, 0.26785714, 0.32142857, 0.375,\
                0.40714286, 0.45714286, 0.5, 0.57142857, 0.64285714,\
                0.71428571, 0.78571429, 0.85714286, 0.92857143, 0.97142857])
        depths = depths.tolist()
    elif arg == '-torques':
        print("tracing over TORQUES")
        qvals = [3, 1801, 1802, 1803, 1804, 1819]
        if magnetism:
            qvals.append(1805)
            qvals.append(1806)
        tag_already_provided = True
        tag = 'torques' + '_'
    elif arg == '-induction':
        print("tracing over INDUCTION QUANTITIES")
        qvals = [1604, 1605, 1609, 1610, 1614, 1615, 1619, 1620, 1623,\
                1624, 1629, 1630, 1601, 1602, 1603, 1606, 1607, 1608,\
                1611, 1612, 1613, 1616, 1617, 1618, 1621, 1622, 1623,\
                1626, 1627, 1628]
        tag_already_provided = True
        tag = 'induction' + '_'

# If the user chose different depths/vals, they must put a unique tag
# on the data
if not tag_already_provided:
    tag = ''
    if user_specified_depths or user_specified_vals:
        tag = input("You chose latitudes/quantities different from the\n\
    default; please enter a unique tag for the output data file:\n")
        tag = tag + '_'

# Set the timetrace savename by the directory, what we are saving, 
# and first and last iteration files for the trace (and optional tag)
savename = dirname_stripped + '_time-latitude_' + tag +\
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
print ('Saving file at ' + savefile)
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
