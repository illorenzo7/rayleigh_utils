# Routine to trace Rayleigh G_Avgs data in time
# Created by: Loren Matilsky
# On: 07/18/2019
############################################################################
# This routine joins multiple traces in time/latitude (hopefully at 
# contiguous intervals, i.e., 
# python time-latitude_join.py [...]_iter1_iter2.pkl  [...]_iter2_iter3.pkl 
#  [...]_iter3_iter4.pkl -->  
#           [...]_iter1_iter4.pkl  in data directory

# Import relevant modules
import numpy as np
import pickle
import sys, os
sys.path.append(os.environ['rapp'])
sys.path.append(os.environ['raco'])
from rayleigh_diagnostics import AZ_Avgs
from common import *

# Find the relevant place to store the data, and create the directory if it
# doesn't already exist

# Make opportunity for command-line args...
args = sys.argv[1:]
n_total_args = len(args)
n_for_cla = 0
for arg in args:
    if arg == '-tag':
        n_for_cla += 2

nfiles = n_total_args - 1 - n_for_cla
print(nfiles)

files = sys.argv[1:nfiles + 1]
dirname = sys.argv[nfiles + 1]
datadir = dirname + '/data/'
dirname_stripped = strip_dirname(dirname)

tag = ''
for i in range(n_total_args):
    arg = args[i]
    if arg == '-tag':
        tag = args[i+1] + '_'

# Read in all the dictionaries to be conjoined
di_list = []
for i in range(nfiles):
    di = get_dict(files[i])
    di_list.append(di)

# Calculate cutoffs for each of the nfiles - 1 "anterior" arrays
niter = 0
cutoffs = np.zeros(nfiles - 1, dtype='int')
for i in range(nfiles - 1):
    di1_iters = di_list[i]['iters']
    di2_iter1 = di_list[i+1]['iters'][0]
    niters1 = len(np.where(di1_iters < di2_iter1)[0])
    niter += niters1
    cutoffs[i] = niters1
niter += di_list[-1]['niter']
# niter is the total number of times

# Initialize the new arrays for vals, times, iters
dummy, nl, nm, nrvals, nq = np.shape(di_list[0]['vals'])
vals = np.zeros((niter, nl, nm, nrvals, nq), dtype='complex')
times = np.zeros(niter)
iters = np.zeros(niter)

length = 0 # Current (just past?) length of the time axis of nonzero values
    # in the various arrays
for i in range(nfiles):
    if i != nfiles - 1:
        vals[length:length + cutoffs[i]] = di_list[i]['vals'][:cutoffs[i]]
        times[length:length + cutoffs[i]] = di_list[i]['times'][:cutoffs[i]]
        iters[length:length + cutoffs[i]] = di_list[i]['iters'][:cutoffs[i]]
        length += cutoffs[i]
    else:
        vals[length:] = di_list[i]['vals']
        times[length:] = di_list[i]['times']
        iters[length:] = di_list[i]['iters']

# Initialize joined dictionary, then change it
di_all = dict(di_list[0])

di_all['vals'] = vals
di_all['times'] = times
di_all['iters'] = iters
di_all['niter'] = niter

iter1, iter2 = di_list[0]['iter1'], di_list[nfiles - 1]['iter2']
di_all['iter1'] = iter1
di_all['iter2'] = iter2

savename = dirname_stripped + '_trace_lmvals_' + tag + str(iter1).zfill(8) +\
        '_' + str(iter2).zfill(8) + '.pkl'
savefile = datadir + savename
print('saving ', savefile)
f = open(savefile, 'wb')
pickle.dump(di_all, f, protocol=4)
f.close()
