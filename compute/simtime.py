# Routine to figure out how much simulation time a given directory contains
# -- where it started, where it ended. 
# Created: 12/23/2018 (but really before)
#
# The simulation time is calculated from the longest possible time
# interval in the G_Avgs data (i.e., first output file to last output file)
import numpy as np
import os, sys
sys.path.append(os.environ['rapp'])
sys.path.append(os.environ['raco'])
from rayleigh_diagnostics import G_Avgs, AZ_Avgs, TransportCoeffs
from reference_tools import equation_coefficients
from common import get_file_lists, get_desired_range, strip_dirname,\
        get_widest_range_file, get_dict
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
    print ("Got timing info from data/*trace_G_Avgs* file")
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
        print ("Got timing info from G_Avgs directory")
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
        print ("Got timing info from AZ_Avgs directory")

simtime = t2 - t1

# Get the thermal diffusion time
try:
    t = TransportCoeffs(dirname + '/transport')
    nu_top = t.nu[0]
    radius = t.radius
except:
    eq = equation_coefficients()
    eq.read(dirname + '/equation_coefficients')
    nu = eq.constants[4]*eq.functions[2]
    nu_top = nu[0]
    radius = eq.radius

shell_depth = np.max(radius) - np.min(radius)
tdt = shell_depth**2/nu_top

simtime_tdt = simtime/tdt

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

print ('t2 - t1: %.1f %s' %(simtime, unit))
print ('t2 - t1: %.2f TDTs' %simtime_tdt)
print ('----------------------')
print ('iter1: %s, iter2: %s' %(str(iter1).zfill(8), str(iter2).zfill(8)))
print ('t1: %.1f, t2: %.1f (%s)' %(t1, t2, unit))
print ('----------------------')

if rotation:
    print('P_rot: %.1f days' %(P_rot/86400))
print('TDT: %.1f days' %(tdt/86400))
