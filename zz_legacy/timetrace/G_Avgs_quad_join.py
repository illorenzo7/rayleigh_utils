# Routine to join Rayleigh time traces
# Created by: Loren Matilsky
# On: 07/18/2019
# Revised: 11/15/2023
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

print(make_bold('starting joined %s with' %dataname))
print(files[0])
di0 = get_dict(files[0])

# figure out the overlap of the values that are output in each dictionary
qv0 = di0['qv']
qv = np.copy(qv0)
for i in range(nfiles - 1):
    qv = np.intersect1d(qv, get_dict(files[i + 1])['qv'])
nq = len(qv)

# also need to update the lookup table
lut = np.zeros_like(di0['lut']) + 4000
lut[qv] = np.arange(nq)

# rearrange first axis in vals to correspond to qv
vals = di0['vals'] # start with arrays from first dictionary and then 
    # append corresponding arrays from all the data files
q_inds0 = np.zeros(nq, dtype=int)
for iq in range(nq):
    q_inds0[iq] = np.argmin(np.abs(qv0 - qv[iq]))
vals = vals[:, q_inds0, ...]

if 'nquad' in dataname:
    vals_full = di0['vals_full'] 
    vals_full = vals_full[:, q_inds0, ...]

times = di0['times']
iters = di0['iters']
iter1, dummy = get_iters_from_file(files[0])

di_all = dict(di0)

# Now append vals, times, and iters with data from joining data files
for i in range(nfiles - 1):
    print(make_bold('appending'))
    print(files[i+1])
    di1 = get_dict(files[i])
    di2 = get_dict(files[i + 1])

    di1_iter2 = di1['iters'][-1]
    di2_iters = di2['iters']
    niters2 = len(np.where(di2_iters > di1_iter2)[0]) # this is the number
            # of NON-OVERLAPPING values in di2_iters
    vals_loc = di2['vals'][-niters2:]

    # build array of (sorted) qv indices the dictionary to average
    q_inds2 = np.zeros(nq, dtype='int')
    for iq in range(nq):
        q_inds2[iq] = np.argmin(np.abs(di2['qv'] - qv[iq]))
    vals_loc = vals_loc[:,q_inds2,...]

    # Now join the dictionary to append
    vals = np.vstack((vals, vals_loc))

    if 'nquad' in dataname:
        vals_full_loc = di2['vals_full'][-niters2:]
        vals_full_loc = vals_full_loc[:, q_inds2,...]
        vals_full = np.vstack((vals_full, vals_full_loc))

    times = np.hstack((times, di2['times'][-niters2:]))
    iters = np.hstack((iters, di2['iters'][-niters2:]))
    if i == nfiles - 2:
        dummy, iter2 = get_iters_from_file(files[i+1])

di_all['vals'] = vals
if 'nquad' in dataname:
    di_all['vals_full'] = vals_full

di_all['times'] = times
di_all['lut'] = lut
di_all['qv'] = qv
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
