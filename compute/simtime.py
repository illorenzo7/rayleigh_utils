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
from rayleigh_diagnostics import G_Avgs, Shell_Slices
from common import get_file_lists, get_desired_range, strip_dirname,\
        get_widest_range_file, get_dict
from get_parameter import get_parameter
from time_scales import compute_Prot, compute_tdt

dirname = sys.argv[1]

# Primary ways to get the 
use_gav = True      # Use the G_Avgs/ directory directly
use_sslice = False  # Use the Shell_Slices/ directory directly
use_gtr = False     # Use the pre-computed G_Avgs time trace
verbose = False

# Get command-line arguments
args = sys.argv[2:]
nargs = len(args)
for i in range(nargs):
    arg = args[i]
    if arg == '-sslice': 
        use_gav = False
        use_sslice = True
        use_gtr = False   
    elif arg == '-gtr':
        use_gav = False
        use_sslice = False
        use_gtr = True
    elif arg == '-v': # verbose
        verbose = True

datadir = dirname + '/data/'

if use_gav:
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
    print ("simtime(): Got timing info from G_Avgs/ directory")

elif use_sslice:
    sslice_dir = dirname + '/Shell_Slices/'
    file_list, int_file_list, nfiles = get_file_lists(sslice_dir)

    f1 = file_list[0]
    f2 = file_list[-1]

    a1 = Shell_Slices(sslice_dir + f1, '')
    a2 = Shell_Slices(sslice_dir + f2, '')

    t1 = a1.time[0]
    t2 = a2.time[-1]

    iter1 = a1.iters[0]
    iter2 = a2.iters[-1]
    print ("simtime(): Got timing info from Shell_Slices/ directory")

elif use_gtr:
    trace_file = get_widest_range_file(datadir, 'trace_G_Avgs')
    di = get_dict(datadir + trace_file)
    times = di['times']
    iters = di['iters']
    t1, t2 = times[0], times[-1]
    iter1, iter2 = iters[0], iters[-1]
    print ("simtime(): Got timing info from data/*trace_G_Avgs* file")

simtime = t2 - t1

# Print run time info in various units (only in Prot if rotation = True)
rotation = get_parameter(dirname, 'rotation')
if rotation:
    Prot = compute_Prot(dirname)
    print ('t2 - t1 = %1.2e Prot' %(simtime/Prot))
# Get the diffusion time
TDT = compute_tdt(dirname)
print ('Delta_t = %.2f TDTs' %(simtime/TDT))
if verbose:
    print ('Delta_t = %1.2e sec' %simtime)
    print ('Delta_t = %.1f days' %(simtime/86400.))
    print ('----------------------')
    print ('iter1 = ', iter1)
    print ('iter2 = ', iter2)
    print ('Delta_iter = %1.2e' %(iter2 - iter1))
    if rotation:
        print('1 Prot = %.1f days' %(Prot/86400.))
    print('1 TDT = %.1f days' %(TDT/86400))
