##################################################################
# Routine to average Rayleigh data in time (generic)
# Author: Loren Matilsky
# Created: 09/13/2019
##################################################################
# This routine joines several time averages of Rayleigh data
# averaging in proportion each one's interval (iter2 - iter1)
# Usage: 
# python timeavg_join.py dirname dir1/data/[dataname]_iter1_iter2.pkl  
# dir2/data/[dataname]_iter2_iter3.pkl dir3/data/[dataname]_iter3_iter4.pkl
#  -->  
# dirname/data/[dataname]_iter1_iter4.pkl
##################################################################

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

# Make opportunity for command-line args
files = []
dirname = sys.argv[1]
datadir = dirname + '/data/'

delete_old_files = True # delete the partial files by default
args = sys.argv[2:]
nargs = len(args)
tag = ''
for i in range(nargs):
    arg = args[i]
    if arg == '--nodel':
        delete_old_files = False
    if arg == '--tag':
        tag = '_' + args[i+1]
    if arg[-4:] == '.pkl':
        files.append(arg)

nfiles = len(files)
dataname = get_dataname_from_file(files[0])
print ('averaging %i %s files' %(nfiles, dataname))
# Read in all the dictionaries to be conjoined
di_list = []
for i in range(nfiles):
    di = get_dict(files[i])
    di_list.append(di)

# calculate the averaging weights
interval_tot = 0.
for i in range(nfiles):
    iter1, iter2 = get_iters_from_file(files[i])
    interval_tot += (iter2 - iter1)
weights = np.zeros(nfiles)
for i in range(nfiles):
    iter1, iter2 = get_iters_from_file(files[i])
    interval_loc = iter2 - iter1
    weights[i] = interval_loc/interval_tot
    if i == 0:
        iter1_glob = iter1
    if i == nfiles - 1:
        iter2_glob = iter2

# average the files together
di0 = di_list[0]

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
vals = np.zeros_like(di0['vals'])
q_inds0 = np.zeros(nq, dtype=int)
for iq in range(nq):
    q_inds0[iq] = np.argmin(np.abs(qv0 - qv[iq]))
vals = vals[..., q_inds0]

if dataname == 'Shell_Spectra':
    lpower = np.zeros_like(di0['lpower'])
    lpower = lpower[..., q_inds0, :]

for i in range(nfiles):
    print('adding %s' %files[i])
    di_loc = di_list[i]
    weight_loc = weights[i]
    print ('weight = %.3f' %weight_loc)
    vals_loc = di_loc['vals']
    # build array of (sorted) qv indices the dictionary to average
    q_inds2 = np.zeros(nq, dtype='int')
    for iq in range(nq):
        q_inds2[iq] = np.argmin(np.abs(di_loc['qv'] - qv[iq]))
    vals_loc = vals_loc[..., q_inds2]

    vals += vals_loc*weight_loc
    if dataname == 'Shell_Spectra':
        lpower_loc = di_loc['lpower']
        lpower_loc = lpower_loc[..., q_inds2, :]
        lpower += lpower_loc*weight_loc

# Initialize joined dictionary, then change it
di = dict(di_list[0])
di['vals'] = vals
di['qv'] = qv
di['lut'] = lut
if dataname == 'Shell_Spectra':
    di['lpower'] = lpower

savename = dataname + tag + '-' +  str(iter1_glob).zfill(8) + '_' +\
        str(iter2_glob).zfill(8) + '.pkl'
savefile = datadir + savename
f = open(savefile, 'wb')
pickle.dump(di, f, protocol=4)
f.close()
print ("Saved joined average in")
print (make_bold(savefile))

# only do this after proper save
if delete_old_files:
    print (make_bold("deleting"))
    for i in range(nfiles):
        fname = files[i]
        print (fname)
        os.remove(fname)
