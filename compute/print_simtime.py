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

# get equation coefficients (for prot and tdt)
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
if kw.gav:
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

elif kw.sslice:
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
    print (f1, " to ", f2)

elif kw.gtr:
    if kw.the_file == None:
        kw.the_file = get_widest_range_file(datadir, 'G_Avgs_trace')
    di = get_dict(the_file)
    times = di['times']
    iters = di['iters']
    t1, t2 = times[0], times[-1]
    iter1, iter2 = iters[0], iters[-1]
    print ("simtime(): Got timing info from data/*trace_G_Avgs* file")
    print ("fname = ", kw.the_file)

# calculate simulation time (seconds)
simtime = t2 - t1

# run time in millions of iterations
print (buff_line)
print ('Delta_t = %7.2f M iter' %((iter2 - iter1)/1.0e6))

# run time in thermal diffusion times
fmt = "%1.3e"
print (buff_line)
print (("Delta_t = " + fmt + " TDTs") %(simtime/eq.tdt))

# Run time, P_rot (if available)
if rotation:
    print (buff_line)
    print (("Delta_t = " + fmt + " rotations") %(simtime/eq.prot))

# more stuff in verbose mode
if kw.verbose:
    print (buff_line)
    print (("Delta_t = " + fmt + " sec") %simtime)
    print (("Delta_t = " + fmt + " days") %(simtime/86400.))

    print (buff_line)
    print ('iter1 = %08i' %iter1)
    print ('iter2 = %08i' %iter2)

    print (("t1 = " + fmt + " TDTs") %(t1/eq.tdt))
    print (("t2 = " + fmt + " TDTs") %(t2/eq.tdt))

    if rotation:
        print (("t1 = " + fmt + " rotations") %(t1/eq.prot))
        print (("t2 = " + fmt + " rotations") %(t2/eq.prot))

# Print the various time scales
print (buff_line)
print (("1 TDT      = " + fmt + " sec = " + fmt + " days" ) %(eq.tdt, eq.tdt/86400.))
if rotation:
    print (("1 rotation = " + fmt + " sec = " + fmt + " days" ) %(eq.prot, eq.prot/86400.))
    print(("1 TDT      = " + fmt + " rotations") %(eq.tdt/eq.prot))
print (buff_line)
