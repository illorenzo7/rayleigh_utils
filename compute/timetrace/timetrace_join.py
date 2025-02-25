# Routine to join Rayleigh time traces
# Created by: Loren Matilsky
# On: 07/18/2019
# Revised: 04/08/2021
############################################################################
# This routine joins multiple traces (hopefully at 
# contiguous intervals, e.g.,
# python timetrace_join.py [...]_iter1_iter2.pkl  [...]_iter2_iter3.pkl 
#  [...]_iter3_iter4.pkl -->  
#           [...]_iter1_iter4.pkl  in directory [wdir]/data

# Import relevant modules
import numpy as np
import pickle
import sys, os
sys.path.append(os.environ['rapp'])
sys.path.append(os.environ['raco'])
from common import *
from grid_util import compute_grid_info, compute_theta_grid

# Find the relevant place to store the data, and create the directory if it
# doesn't already exist
files = []
dirname = sys.argv[1]

delete_old_files = False # delete the partial files by specifying --del
args = sys.argv[2:]
nargs = len(args)
for i in range(nargs):
    arg = args[i]
    if arg == '--del':
        delete_old_files = True
    if arg[-4:] == '.pkl':
        files.append(arg)
nfiles = len(files)

# get the relative data directory
datadir_rel = get_datadir_from_file(files[0])

# the dataname
dataname = get_dataname_from_file(files[0])

# full datadir
datadir = dirname + '/' + datadir_rel
# create it if non-existent
if not os.path.isdir(datadir):
    os.makedirs(datadir)

# for timelat or mtimelat, see if we need to interpolate
if 'timelat' in dataname:
    # for now this will only work if I don't change the radial resolution
    interp = False # only interpolate if at least one axis
    # is a different size than the rest
    nts = np.zeros(nfiles, dtype='int')
    di_list = []
    for i in range(nfiles):
        di_list.append(get_dict(files[i]))
        nt_loc = np.shape(di_list[i]['vals'])[1]
        nts[i] = nt_loc
    nt_min, nt_max = np.min(nts), np.max(nts)
    if nt_min < nt_max:
        print('interpolating all horizontal grids onto nt=%i' %nt_min)
        interp = True
        tt_min, tt_min_tw = compute_theta_grid(nt_min)
else:
    interp = False

# get the first array to add on the others
print(make_bold('starting joined %s with' %dataname))
print(files[0])

# start (maybe interpolating) the first array

di0 = di_list[0]
if interp:
    if nts[0] > nt_min:
        print('interpolating nt = %i onto the coarser grid nt=%i' %(nts[0], nt_min))
        tt_loc, tw_loc = compute_theta_grid(nts[0])
        di0['vals'] = vals = interp_nd(di0['vals'], tt_loc, tt_min, axis=1) 
else:
    vals = di0['vals'] # start with arrays from first dictionary and then 
        # append corresponding arrays from all the data files

# some other values to join in this case
if 'nquad' in dataname:
    vals_full = di0['vals_full'] 

# basic time data, to join
times = di0['times']
iters = di0['iters']
iter1, dummy = get_iters_from_file(files[0])

# now create the final dictionary with all the others to be appended
di_all = dict(di0)

# Now append vals, times, and iters with data from joining data files
for i in range(nfiles - 1):
    print(make_bold('appending'))
    print(files[i+1])

    # join these two dictionaries
    di1 = di_list[i]
    di2 = di_list[i+1]

    di1_iter2 = di1['iters'][-1]
    di2_iters = di2['iters']
    niters2 = len(np.where(di2_iters > di1_iter2)[0]) # this is the number
            # of NON-OVERLAPPING values in di2_iters
            
    # Now join the dictionary to append
    vals_append = di2['vals'][-niters2:]

    # see if we need to interpolate
    if interp:
        nt_loc = nts[i+1]
        if nt_loc > nt_min:
            print('interpolating nt = %i onto the coarser grid nt=%i' %(nt_loc, nt_min))
            tt_loc, tw_loc = compute_theta_grid(nt_loc)
            vals_append = interp_nd(vals_append, tt_loc, tt_min, axis=1)

    di_all['vals'] = np.vstack((di_all['vals'], vals_append))
    if 'nquad' in dataname:
        vals_full = np.vstack((vals_full, di2['vals_full'][-niters2:]))
    times = np.hstack((times, di2['times'][-niters2:]))
    iters = np.hstack((iters, di2['iters'][-niters2:]))
    if i == nfiles - 2:
        dummy, iter2 = get_iters_from_file(files[i+1])

if 'nquad' in dataname:
    di_all['vals_full'] = vals_full

di_all['times'] = times
di_all['iters'] = iters
di_all['iter1'] = iter1
di_all['iter2'] = iter2

savename = dataname + '-' + str(iter1).zfill(8) +\
        '_' + str(iter2).zfill(8) + '.pkl'
savefile = datadir + savename
f = open(savefile, 'wb')
pickle.dump(di_all, f, protocol=4)
f.close()
print ("Saved joined trace in")
print (make_bold(savefile))

# only do this after proper save
if delete_old_files:
    print (make_bold("deleting"))
    for i in range(nfiles):
        fname = files[i]
        print (fname)
        os.remove(fname)
