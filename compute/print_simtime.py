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
from common import *
from cla_util import *

# read in args
args = sys.argv 
clas0, clas = read_clas(args)
dirname = clas0.dirname
rotation = clas0.rotation

# get equation coefficients (for t_omega and t_kappa)
eq = get_eq(dirname)

# Primary ways to get the simulation time
kw_default = dotdict()
kw_default.gav = True      # Use the G_Avgs/ directory directly (default)
kw_default.sslice = False  # Use the Shell_Slices/ directory directly
kw_default.gtr = False     # Use the pre-computed G_Avgs time trace
kw_default.verbose = False
kw_default.the_file = None

# change kwargs with clas
kw = update_dict(kw_default, clas)

# now print the simulation time
print (buff_line)
if kw.sslice:
    sslice_dir = dirname + '/Shell_Slices/'
    file_list, int_file_list, nfiles = get_file_lists_all(sslice_dir)

    f1 = file_list[0]
    f2 = file_list[-1]

    a1 = Shell_Slices(sslice_dir + f1, '')
    a2 = Shell_Slices(sslice_dir + f2, '')

    t1 = a1.time[0]
    t2 = a2.time[-1]

    iter1 = a1.iters[0]
    iter2 = a2.iters[-1]
    print ("simtime(): Got timing info from Shell_Slices/ directory")
    print ("iters:", f1, " to ", f2)
    print ("times: %1.4e to %1.4e (simulation units)" %(t1, t2))

elif kw.gtr:
    if kw.the_file == None:
        kw.the_file = get_widest_range_file(dirname + '/data/', 'G_Avgs_trace')
    di = get_dict(kw.the_file)
    times = di['times']
    iters = di['iters']
    t1, t2 = times[0], times[-1]
    iter1, iter2 = iters[0], iters[-1]
    print ("simtime(): Got timing info from data/*trace_G_Avgs* file")
    print ("fname = ", kw.the_file)
    print ("iters:", str(iter1).zfill(8), " to ", str(iter2).zfill(8))
    print ("times: %1.4e to %1.4e (simulation units)" %(t1, t2))

elif kw.gav:
    gav_dir = dirname + '/G_Avgs/'
    file_list, int_file_list, nfiles = get_file_lists_all(gav_dir)

    f1 = file_list[0]
    f2 = file_list[-1]

    a1 = G_Avgs(gav_dir + f1, '')
    a2 = G_Avgs(gav_dir + f2, '')

    t1 = a1.time[0]
    t2 = a2.time[-1]

    iter1 = a1.iters[0]
    iter2 = a2.iters[-1]
    print ("simtime(): Got timing info from G_Avgs/ directory")
    print ("iters:", f1, " to ", f2)
    print ("times: %1.4e to %1.4e (simulation units)" %(t1, t2))

# calculate simulation time (sim units)
simtime = t2 - t1

# run time in millions of iterations
print (buff_line)
print ('Delta_t = %7.3f M iter' %((iter2 - iter1)/1.0e6))

# run time in thermal diffusion times
fmt = "%1.4e"
print (buff_line)
print (("Delta_t = " + fmt + " t_kappa") %(simtime/eq.tkappa))

# Run time, P_rot (if available)
if rotation:
    print (buff_line)
    print (("Delta_t = " + fmt + " t_omega") %(simtime/eq.tomega))

# more stuff in verbose mode
if kw.verbose:
    print (buff_line)
    print ('iter1 = %08i' %iter1)
    print ('iter2 = %08i' %iter2)

    print (("t1 = " + fmt + " t_kappa") %(t1/eq.tkappa))
    print (("t2 = " + fmt + " t_kappa") %(t2/eq.tkappa))

    if rotation:
        print (("t1 = " + fmt + " t_omega") %(t1/eq.tomega))
        print (("t2 = " + fmt + " t_omega") %(t2/eq.tomega))

# Print the various time scales
print (buff_line)
print (("1 t_kappa      = " + fmt + " simulation units") %eq.tkappa)
if rotation:
    print (("1 t_omega = " + fmt + " simulation units") %eq.tomega)
    print(("1 t_kappa      = " + fmt + " t_omega") %(eq.tkappa/eq.tomega))
print (buff_line)
