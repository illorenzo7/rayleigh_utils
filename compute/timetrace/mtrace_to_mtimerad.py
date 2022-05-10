##################################################################
# Routine to trace AZ_Avgs in time/latitude space (pick different radii)
# Author: Loren Matilsky
# Created: 02/27/2019
# Parallelized: 12/12/2020
##################################################################
# This routine computes the trace in time/latitude (by default) or
# time/radius of quantities in the 
# AZ_Avgs data for a particular simulation. 
##################################################################

# improt modules
import numpy as np
# data type and reading function
import sys, os
sys.path.append(os.environ['rapp'])
sys.path.append(os.environ['raco'])
from common import *
from cla_util import *
import pickle

# get CLAS
args = sys.argv
clas0, clas = read_clas(args)
dirname = clas0['dirname']
magnetism = clas0['magnetism']
kwargs_default = dict({'latvals': default_latvals, 'qvals': None, 'groupname': 'b', 'mval': 1, 'mmax': None})
kwargs = update_dict(kwargs_default, clas)

if kwargs.qvals is None: # it's a quantity group
    groupname = kwargs.groupname
    qgroup = get_quantity_group(groupname, magnetism)
    qvals = qgroup['qvals']
else:
    qvals = kwargs.qvals
    groupname = input("choose a groupname to save your data: ")

nq = len(qvals)

mval = kwargs['mval']
mmax = kwargs.mmax # this just tells us which directory to get the data from

# get grid information
di_grid = get_grid_info(dirname)
rr = di_grid['rr']
tt_lat = di_grid['tt_lat']

# get indices associated with desired sample vals
samplevals = kwargs['latvals']
sampleaxis = tt_lat

isamplevals = []
for sampleval in samplevals:
    isamplevals.append(np.argmin(np.abs(sampleaxis - sampleval)))

isamplevals = np.array(isamplevals)

# recompute the actual sample values we get
samplevals = sampleaxis[isamplevals]
nsamplevals = len(samplevals)
print ('converting mtrace data to mtimerad')
print ('groupname =', groupname)
print ("qvals = " + arr_to_str(qvals, "%i"))
print ("sampling lats = " + arr_to_str(samplevals, '%.1f'))

# get mtrace directory
if mmax is None:
    datadir_mtrace = clas0['datadir'] + 'mtrace/'
else:
    datadir_mtrace = clas0['datadir'] + ('mtrace_mmax%03i/' %mmax)

# get the levels associated with Shell_Slices
radlevs = get_slice_levels(dirname)
irvals = radlevs.inds
rr_avail = radlevs.radius
nr_avail = len(rr_avail)
print ("irvals avail =", irvals)
print ("rvals  avail =", rr_avail/rsun)

# get the mtrace data
firstfile = True
for irval in range(nr_avail):
    # vals (initialized later) should have shape
    # (ntimes, nlats, nr_avail, len(qvals))
    for iqval in  range(nq):
        qval = qvals[iqval]
        # names of datafile
        dataname_mtrace = ('mtrace_qval%04i_irval%02i' %(qval, irval)) +\
                clas0['tag']
        # get mtrace data (local)
        the_file = get_widest_range_file(datadir_mtrace, dataname_mtrace)
        print (buff_line)
        print ('reading ' + the_file)
        print ('rval =  %.3f' %(rr_avail[irval]/rsun))
        di = get_dict(the_file)

        if firstfile:
            # only need to read in certain things
            # (length of time axis, first and last iters, etc.) once
            times = di['times']
            iters = di['iters']
            ntimes = len(times)
            vals = np.zeros((ntimes, nsamplevals, nr_avail, nq), 'complex')
            # get the first and last iters for the savename
            iter1, iter2 = get_iters_from_file(the_file)

            firstfile = False
       
        # add in relevant vals (only at the lats we wish to sample)
        vals_loc = di['vals'][:, mval, isamplevals]
        vals[:, :, irval, iqval] = vals_loc

# create data directory if it doesn't already exist
basename = 'mtimerad'
datadir = clas0['datadir'] + basename + '_mval%03i/' %mval
if not os.path.isdir(datadir):
    os.makedirs(datadir)

basename += '_' + groupname
savename = basename + clas0['tag'] + '-' +\
        str(iter1).zfill(8) + '_' + str(iter2).zfill(8) + '.pkl'
savefile = datadir + savename

# save the data
f = open(savefile, 'wb')
di_sav = dict({'vals': vals, 'times': times, 'iters': iters, 'qvals': qvals, 'latvals': samplevals, 'mval': mval, 'rr_avail': rr_avail})
#pickle.dump(di_sav, f, protocol=4)
f.close()

print ('data saved at ')
print (make_bold(savefile))
