# Routine to figure out how much simulation time a given directory contains
# -- where it started, where it ended. 
# Created: 12/23/2018 (but really before)
#
# The simulation time is calculated from the longest possible time
# interval in the G_Avgs data (i.e., first output file to last output file)
import numpy as np
import os, sys
sys.path.append(os.environ['rapp'])
sys.path.append(os.environ['co'])
from rayleigh_diagnostics import G_Avgs, AZ_Avgs
from common import get_file_lists, get_desired_range, strip_dirname,\
        get_widest_range_file
from get_parameter import get_parameter

dirname = sys.argv[1]
dirname_stripped = strip_dirname(dirname)

datadir = dirname + '/data/'

try: # first try to get data from pre-computed time trace:
    trace_file = get_widest_range_file(datadir, 'trace_G_Avgs')
    di = get_dict(datadir + trace_file)
    times = di['times']
    iters = di['iters']
    t1, t2 = times[0], times[-1]
    iter1, iter2 = iters[0], iters[-1]
except:
    try:
        gavg_dir = dirname + '/G_Avgs/'
        file_list, int_file_list, nfiles = get_file_lists(gavg_dir)

        f1 = file_list[0]
        f2 = file_list[-1]

        a1 = G_Avgs(gavg_dir + f1, '')
        a2 = G_Avgs(gavg_dir + f2, '')

        t1 = a1.time[0]
        t2 = a2.time[-1]

        iter1 = a1.iters[0]
        iter2 = a2.iters[-1]
    except: # and if that fails, try the AZ_Avgs
        azavg_dir = dirname + '/AZ_Avgs/'
        file_list, int_file_list, nfiles = get_file_lists(azavg_dir)

        f1 = file_list[0]
        f2 = file_list[-1]

        a1 = AZ_Avgs(azavg_dir + f1, '')
        a2 = AZ_Avgs(azavg_dir + f2, '')

        t1 = a1.time[0]
        t2 = a2.time[-1]

        iter1 = a1.iters[0]
        iter2 = a2.iters[-1]

simtime = t2 - t1

rotation = True
try:
    Om0 = get_parameter(dirname, 'angular_velocity')
    P_rot = 2*np.pi/Om0
    simtime /= P_rot
    t1 /= P_rot
    t2 /= P_rot
    unit = 'P_rot'
except:
    simtime /= 86400.
    t1 /= 86400.
    t2 /= 86400.
    unit = 'days'
    rotation = False

print ('iter1: %s, iter2: %s' %(f1, f2))
print ('t1: %.1f, t2: %.1f (%s)' %(t1, t2, unit))
print ('t2 - t1: %.1f (%s)' %(simtime, unit))

if rotation:
    print('P_rot: %.1f days' %(P_rot/86400))
