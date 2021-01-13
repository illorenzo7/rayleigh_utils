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

dirname = sys.argv[1]

# Primary ways to get the 
use_gav = True      # Use the G_Avgs/ directory directly
use_sslice = False  # Use the Shell_Slices/ directory directly
use_gtr = False     # Use the pre-computed G_Avgs time trace
verbose = False
trace_file = None
mag = False # don't use magnetic diffusion time by default
tach = False #  by default don't calculate two diffusion times

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
    elif arg == '-usefile':
        trace_file = args[i+1]
        trace_file = Shell_Avgs_file.split('/')[-1]
    elif arg == '-mag':
        mag = True
    elif arg == '-tach':
        tach = True

datadir = dirname + '/data/'

print ("-------------------------------------------------")
if use_gav:
    gavg_dir = dirname + '/G_Avgs/'
    file_list, int_file_list, nfiles = get_file_lists(gavg_dir)

    the_tuple = get_desired_range(int_file_list, args)
    if the_tuple is None:
        index_first, index_last = 0, -1 # by default get the WHOLE sim time
    else:
        index_first, index_last = the_tuple

    f1 = file_list[index_first]
    f2 = file_list[index_last]

    a1 = G_Avgs(gavg_dir + f1, '')
    a2 = G_Avgs(gavg_dir + f2, '')

    t1 = a1.time[0]
    t2 = a2.time[-1]

    iter1 = a1.iters[0]
    iter2 = a2.iters[-1]
    print ("simtime(): Got timing info from G_Avgs/ directory")
    print (f1, " to ", f2)

elif use_sslice:
    sslice_dir = dirname + '/Shell_Slices/'
    file_list, int_file_list, nfiles = get_file_lists(sslice_dir)

    the_tuple = get_desired_range(int_file_list, args)
    if the_tuple is None:
        index_first, index_last = 0, -1 # by default get the WHOLE sim time
    else:
        index_first, index_last = the_tuple

    f1 = file_list[index_first]
    f2 = file_list[index_last]

    a1 = Shell_Slices(sslice_dir + f1, '')
    a2 = Shell_Slices(sslice_dir + f2, '')

    t1 = a1.time[0]
    t2 = a2.time[-1]

    iter1 = a1.iters[0]
    iter2 = a2.iters[-1]
    print ("simtime(): Got timing info from Shell_Slices/ directory")
    print (f1, " to ", f2)

elif use_gtr:
    if trace_file == None:
        trace_file = get_widest_range_file(datadir, 'trace_G_Avgs')
    di = get_dict(datadir + trace_file)
    times = di['times']
    iters = di['iters']
    t1, t2 = times[0], times[-1]
    iter1, iter2 = iters[0], iters[-1]
    print ("simtime(): Got timing info from data/*trace_G_Avgs* file")
    print ("fname = ", trace_file)

# calculate simulation time (seconds)
simtime = t2 - t1

# Get viscous/thermal diffusion time(s), and magnetic diffusion time
# (if present)
if tach:
    VDTrz, VDTcz, VDT = compute_tdt(dirname, visc=True, tach=True)
    TDTrz, TDTcz, TDT = compute_tdt(dirname, tach=True)
    if mag:
        MDTrz, MDTcz, MDT = compute_tdt(dirname, mag=True, tach=True)
else:
    VDT = compute_tdt(dirname, visc=True)
    TDT = compute_tdt(dirname)
    if mag:
        MDT = compute_tdt(dirname, mag=True)

Pr = TDT/VDT

if mag: # Get the thermal Prandtl number
    # (assumes same scaling for all diffusions)
    Prm = MDT/TDT

# Print Prandtl numbers
fmt = "%7.2f"
print ("-------------------------------------------------")
print ("Pr = %.1f" %Pr)
if mag:
    print ("Pr_m = %.1f (mag = True)" %Prm)
print ("-------------------------------------------------")

# run time in millions of iterations
print ('Delta_t = %7.2f M iter' %((iter2 - iter1)/1.0e6))
print ("-------------------------------------------------")

# Run time, viscous DTs
if tach:
    print (("Delta_t = " + fmt + " RZ VDTs") %(simtime/VDTrz))
    print (("Delta_t = " + fmt + " CZ VDTs") %(simtime/VDTcz))
print (("Delta_t = " + fmt + " VDTs") %(simtime/VDT))
print ("-------------------------------------------------")

# Run time, thermal DTs (only if Prandtl number != 1)
tol = 1.0e-12
if np.abs(Pr - 1.0) > tol:
    if tach:
        print (("Delta_t = " + fmt + " RZ TDTs") %(simtime/TDTrz))
        print (("Delta_t = " + fmt + " CZ TDTs") %(simtime/TDTcz))
    print (("Delta_t = " + fmt + " TDTs") %(simtime/TDT))
    print ("-------------------------------------------------")

# Run time, magnetic DTs (if available)
if mag:
    if np.abs(Prm - 1.0) > tol:
        if tach:
            print (("Delta_t = " + fmt + " RZ MDTs") %(simtime/MDTrz))
            print (("Delta_t = " + fmt + " CZ MDTs") %(simtime/MDTcz))
        print (("Delta_t = " + fmt + " MDTs") %(simtime/MDT))
        print ("-------------------------------------------------")

# Run time, P_rot (if available)
fmt = "%7.1f"
rotation = get_parameter(dirname, 'rotation')
if rotation:
    Prot = compute_Prot(dirname)
    print (("Delta_t = " + fmt + " P_rot") %(simtime/Prot))
    print ("-------------------------------------------------")

if verbose:
    print ('Delta_t = %7.2e sec' %simtime)
    print ('Delta_t = %7.2e days' %(simtime/86400.))
    print ("-------------------------------------------------")

    print ('iter1 = %08i' %iter1)
    print ('iter2 = %08i' %iter2)
    print ("-------------------------------------------------")

    fmt = "%7.2f"
    print (("t1 = " + fmt + " VDTs") %(t1/VDT))
    print (("t2 = " + fmt + " VDTs") %(t2/VDT))
    print ("-------------------------------------------------")

    print (("t1 = " + fmt + " TDTs") %(t1/TDT))
    print (("t2 = " + fmt + " TDTs") %(t2/TDT))
    print ("-------------------------------------------------")
    if mag:
        print (("t1 = " + fmt + " MDTs") %(t1/MDT))
        print (("t2 = " + fmt + " MDTs") %(t2/MDT))
        print ("-------------------------------------------------")
    if rotation:
        fmt = "%7.1f"
        print (("t1 = " + fmt + " P_rot") %(t1/Prot))
        print (("t2 = " + fmt + " P_rot") %(t2/Prot))
        print ("-------------------------------------------------")

# Print the various time scales
if rotation:
    print('1 Prot = %.2f days' %(Prot/86400.))
    print('Omega_0 = 2*pi/Prot = %.1f nHz' %(1.0/Prot*1.0e9))
    print('Omega_0 = 2*pi/Prot = %1.1e s^{-1}' %(2*np.pi/Prot))
    print ("-------------------------------------------------")

if rotation:
    print('1 VDT    = %7.1f P_rot = %7.1f days' %(VDT/Prot, VDT/86400.))
else:
    print('1 VDT    = %7.1f P_rot' %(VDT/Prot))
if tach:
    if rotation:
        print('1 RZ VDT = %7.1f P_rot = %7.1f days' %(VDTrz/Prot,\
                VDTrz/86400.))
    else:
        print('1 RZ VDT = %7.1f P_rot' %(VDTrz/Prot))
    if rotation:
        print('1 CZ VDT = %7.1f P_rot = %7.1f days' %(VDTcz/Prot,\
                VDTcz/86400.))
    else:
        print('1 CZ VDT = %7.1f P_rot' %(VDTcz/Prot))
print ("-------------------------------------------------")

# only print for Pr != 1
if np.abs(Pr - 1.0) > tol:
    if rotation:
        print('1 TDT    = %7.1f P_rot = %7.1f days' %(TDT/Prot, TDT/86400.))
    else:
        print('1 TDT    = %7.1f P_rot' %(TDT/Prot))
    if tach:
        if rotation:
            print('1 RZ TDT = %7.1f P_rot = %7.1f days' %(TDTrz/Prot,\
                    TDTrz/86400.))
        else:
            print('1 RZ TDT = %7.1f P_rot' %(TDTrz/Prot))
        if rotation:
            print('1 CZ TDT = %7.1f P_rot = %7.1f days' %(TDTcz/Prot,\
                    TDTcz/86400.))
        else:
            print('1 CZ TDT = %7.1f P_rot' %(TDTcz/Prot))
    print ("-------------------------------------------------")

if mag:
    # only print for Pr != 1
    if np.abs(Prm - 1.0) > tol:
        if rotation:
            print('1 MDT    = %7.1f P_rot = %7.1f days' %(MDT/Prot,\
                    MDT/86400.))
        else:
            print('1 MDT = %7.1f P_rot = %7.1f' %(MDT/Prot))
        if tach:
            if rotation:
                print('1 RZ MDT = %7.1f P_rot = %7.1f days' %(MDTrz/Prot,\
                        MDTrz/86400.))
            else:
                print('1 RZ MDT = %7.1f P_rot' %(MDTrz/Prot))
            if rotation:
                print('1 CZ MDT = %7.1f P_rot = %7.1f days' %(MDTcz/Prot,\
                        MDTcz/86400.))
            else:
                print('1 CZ MDT = %7.1f P_rot' %(MDTcz/Prot))
        print ("-------------------------------------------------")
