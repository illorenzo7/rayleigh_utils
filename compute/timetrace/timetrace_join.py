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

# Find the relevant place to store the data, and create the directory if it
# doesn't already exist
files = []
dirname = sys.argv[1]
datadir = dirname + '/data/' # data subdirectory of output directory

delete_old_files = True # delete the partial files by default
args = sys.argv[2:]
nargs = len(args)
for i in range(nargs):
    arg = args[i]
    if arg == '--nodel':
        delete_old_files = False
    if arg[-4:] == '.pkl':
        files.append(arg)
nfiles = len(files)
dataname = get_dataname_from_file(files[0])
print(make_bold('starting joined %s with' %dataname))
print(files[0])
di0 = get_dict(files[0])
vals = di0['vals'] # start with arrays from first dictionary and then 
    # append corresponding arrays from all the data files
if dataname[:17] == 'G_Avgs_trace_2dom':
    vals_cz = di0['vals_cz'] 
    vals_rz = di0['vals_rz'] 

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

    # Now join the dictionary to append
    vals = np.vstack((vals, di2['vals'][-niters2:]))
    if dataname[:17] == 'G_Avgs_trace_2dom':
        vals_cz = np.vstack((vals_cz, di2['vals_cz'][-niters2:]))
        vals_rz = np.vstack((vals_rz, di2['vals_rz'][-niters2:]))
    times = np.hstack((times, di2['times'][-niters2:]))
    iters = np.hstack((iters, di2['iters'][-niters2:]))
    if i == nfiles - 2:
        dummy, iter2 = get_iters_from_file(files[i+1])

di_all['vals'] = vals
if dataname[:17] == 'G_Avgs_trace_2dom':
    di_all['vals_cz'] = vals_cz
    di_all['vals_rz'] = vals_rz

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
