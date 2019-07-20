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
sys.path.append(os.environ['rasource'] + '/post_processing')
sys.path.append(os.environ['co'])
from rayleigh_diagnostics import AZ_Avgs
from common import get_file_lists, get_desired_range, strip_dirname,\
        strip_filename, get_dict
from get_parameter import get_parameter

# Find the relevant place to store the data, and create the directory if it
# doesn't already exist

files = sys.argv[1:-1]
dirname = sys.argv[-1]
datadir = dirname + '/data/'
dirname_stripped = strip_dirname(dirname)
nfiles = len(files)
di0 = get_dict(files[0])
vals = di0['vals']
times = di0['times']
iters = di0['iters']
niter = di0['niter']
iter1 = di0['iter1']
iter2 = di0['iter2']

di_all = dict(di0)

for i in range(nfiles - 1):
    di1 = get_dict(files[i])
    di2 = get_dict(files[i + 1])

    di1_iters = di1['iters']
    di2_iter1 = di2['iters'][0]
    niters1 = len(np.where(di1_iters < di2_iter1)[0])

    vals = np.vstack((vals[:niters1], di2['vals']))
    times = np.hstack((times[:niters1], di2['times']))
    iters = np.hstack((iters[:niters1], di2['iters']))
    niter += (di2['niter'] - (len(di1_iters) - niters1))
    if i == nfiles - 2:
        iter2 = di2['iter2']

di_all['vals'] = vals
di_all['times'] = times
di_all['iters'] = iters
di_all['niter'] = niter
di_all['iter1'] = iter1
di_all['iter2'] = iter2


savename = dirname_stripped + '_time-latitude_' + str(iter1).zfill(8) +\
        '_' + str(iter2).zfill(8) + '.pkl'
savefile = datadir + savename
f = open(savefile, 'wb')
pickle.dump(di_all, f, protocol=4)
