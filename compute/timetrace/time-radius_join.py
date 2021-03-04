# Routine to join several time-radius data sets together
# Created by: Loren Matilsky
# On: 10/21/2020
############################################################################
# This routine joins multiple traces in time/radius (hopefully at 
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
files = []
dirname = sys.argv[1]
datadir = dirname + '/data/' # data subdirectory of output directory
dirname_stripped = strip_dirname(dirname)

# CLAs
tag = ''
delete_old_files = True # delete the partial files by default
args = sys.argv[2:]
nargs = len(args)
for i in range(nargs):
    arg = args[i]
    if arg == '-tag':
        tag = args[i+1] + '_'
    if arg == '-nodel':
        delete_old_files = False
    if arg[-4:] == '.pkl':
        files.append(arg)
nfiles = len(files)

# Find the relevant place to store the data, and create the directory if it
# doesn't already exist

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
dummy, nlats, nr, nq = np.shape(di_list[0]['vals'])
vals = np.zeros((niter, nlats, nr, nq))
times = np.zeros(niter)
iters = np.zeros(niter)

length = 0 # Current (just past?) length of the time axis of nonzero values
    # in the various arrays
for i in range(nfiles):
    print(make_bold('appending'))
    print(files[i])
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

savename = dirname_stripped + '_time-radius_' + tag + str(iter1).zfill(8) +\
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
