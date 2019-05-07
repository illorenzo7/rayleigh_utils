# NEEDS UPDATED DOCS AND TESTING
# Routine to figure out how much simulation time a given directory contains
# Created: 12/23/2018 (but really before)
#
# By default the simulation time is calculated from the longest possible time
# interval in the G_Avgs data (i.e., first output file to last output file)
# User can specify another time interval in the standard way, e.g.,
# -n 10 (last 10 G_Avgs files)
# -range iter1 iter2 (files closest to iter1 to iter2 range)
# -centerrange iter0 nfiles (nfiles centered about iter0)
# etc., etc.
import numpy as np
import os, sys
sys.path.append(os.environ['rapp'])
sys.path.append(os.environ['co'])
from rayleigh_diagnostics import G_Avgs
from common import get_file_lists, get_desired_range,\
        strip_dirname
from get_parameter import get_parameter

dirname = sys.argv[1]
dirname_stripped = strip_dirname(dirname)
gavg_dir = dirname + '/G_Avgs/'

# Get all files in gavg_dir and their integer countparts
file_list, int_file_list, nfiles = get_file_lists(gavg_dir)

# Read in CLAs
args = sys.argv[2:]
nargs = len(args)

if (nargs == 0):
    index_first, index_last = 0, nfiles - 1
else:
    index_first, index_last = get_desired_range(int_file_list,\
            args)

f1 = file_list[index_first]
f2 = file_list[index_last]

a1 = G_Avgs(gavg_dir + f1, '')
a2 = G_Avgs(gavg_dir + f2, '')

niter1 = a1.niter
niter2 = a2.niter

simtime_sec = a2.time[niter2 - 1] - a1.time[0]
simtime_days = simtime_sec/86400.
starttime = a1.time[0]/86400.
rotation = True
try:
    Om0 = get_parameter(dirname, 'angular_velocity')
    P_rot = 2*np.pi/Om0
    simtime_rot = simtime_sec/P_rot
except:
    rotation = False

# Compute thermal diffusion time
ktop = get_parameter(dirname, 'kappa_top')
try:
    rmin = get_parameter(dirname, 'rmin')
    rmax = get_parameter(dirname, 'rmax')
except: # two domains stitched together
    domain_bounds = get_parameter(dirname, 'domain_bounds')
    rmin = np.min(domain_bounds)
    rmax = np.max(domain_bounds)
depth = rmax - rmin
tdt = depth**2/ktop
simtime_diff = simtime_sec/tdt


print ('Simulation %s has run from iter1=%s to iter2=%s' \
        %(dirname_stripped, f1, f2))
print ('Started at %.1f days = %.1f yr'\
        %(starttime,starttime/365.24))
print('Simulation %s has run for a total of %.1f days'\
        %(dirname_stripped, simtime_days))
print('Or %.1f years' %(simtime_days/365.24))
print('Kappa_top = %.1e cm^2/sec; shell_depth=%.1e cm' %(ktop, depth))
print('T_diff = %.1e sec = %.1e days = %.1f yr'\
        %(tdt, tdt/86400, tdt/86400/365.24))
print('Simulation %s has run for a total of %.1f thermal diffusions times'\
        %(dirname_stripped, simtime_diff))
if rotation:
    print('Frame rotation is %.1e rad/sec; rot. period is %.1f days'\
            %(Om0, P_rot/86400))
    print('Simulation %s has run for a total of %.1f rotation periods'\
        %(dirname_stripped, simtime_rot))
else:
    print('Simulation %s is nonrotating' %dirname_stripped)

