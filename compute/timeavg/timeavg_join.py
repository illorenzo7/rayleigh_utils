# Created by: Loren Matilsky
# On: 09/13/2019
# Routine to average Rayleigh time-averaged data subsets.
# Weights each average in proportion to "count" (time-averaging interval).
# Averages only the intersection of the different qv's updating 
# "nq" and "lut" appropriately in the output dictionary.
# Usage: python timeavg_join.py dir1/data/[...]_iter1_iter2.pkl  
# dir2/data/[...]_iter2_iter3.pkl dir3/data/[...]_iter3_iter4.pkl
# dirname/ -->  
# dirname/data/[...]_iter1_iter4.pkl
# Assumes default data type is AZ_Avgs. To specify another datatype, use
# -dtype key, where key = 
# specav: Shell_Spectra
# shav: Shell_Avgs
# gav: G_Avgs
# merav: Meridional_Slices

# Import relevant modules
import numpy as np
import pickle
import sys, os
sys.path.append(os.environ['rapp'])
sys.path.append(os.environ['raco'])
from rayleigh_diagnostics import AZ_Avgs
from common import *
from get_parameter import get_parameter

# Find the relevant place to store the data, and create the directory if it
# doesn't already exist

# Make opportunity for command-line args...
args = sys.argv[1:]
n_total_args = len(args)
n_for_cla = 0
for arg in args:
    if arg == '-tag':
        n_for_cla += 2
    elif arg == '-dtype':
        n_for_cla += 2

nfiles = n_total_args - 1 - n_for_cla

files = sys.argv[1:nfiles + 1]
dirname = sys.argv[nfiles + 1]
print('dirname = ', dirname)
datadir = dirname + '/data/'
dirname_stripped = strip_dirname(dirname)

tag = ''
dtype = 'azav'
# Other choices are gav, shav, specav, merav, enstr
for i in range(n_total_args):
    arg = args[i]
    if arg == '-tag':
        tag = args[i+1] + '_'
    elif arg == '-dtype':
        dtype = args[i+1]

# Read in all the dictionaries to be conjoined
di_list = []
for i in range(nfiles):
    di = get_dict(files[i])
    di_list.append(di)
di0 = di_list[0]

# Get intersection of qv
qv = np.copy(di_list[0]['qv'])
for i in range(nfiles - 1):
    qv = np.intersect1d(qv, di_list[i + 1]['qv'])
nq = len(qv)

# Create new data array of zeros, overwriting nq with value from intersection
if dtype == 'azav':
    nt, nr, nq_old = np.shape(di0['vals'])
    vals = np.zeros((nt, nr, nq))
elif dtype == 'gav':
    nq_old, = np.shape(di0['vals'])
    vals = np.zeros((nq,))
elif dtype == 'shav':
    nr, nq_old = np.shape(di0['vals'])
    vals = np.zeros((nr, nq))
elif dtype == 'merav':
    nphi, nt, nr, nq_old = np.shape(di0['vals'])
    vals = np.zeros((nphi, nt, nr, nq))
elif dtype == 'specav':
    nl, nm, nr, nq_old = np.shape(di0['fullpower'])
    fullpower = np.zeros((nl, nm, nr, nq))
    lpower = np.zeros((nl, nr, nq, 3))
elif dtype == 'enstr':
    nphi, nt, nr = np.shape(di0['vals'])
    vals = np.zeros((nphi, nt, nr))

lut = np.zeros_like(di_list[0]['lut'])
lut[qv] = np.arange(nq)

count = 0
for i in range(nfiles):
    print('Averaging data from %s ...' %files[i])
    this_di = di_list[i]
    this_count = this_di['count']
    print('count = %i ...' %this_count)
    this_lut = this_di['lut']
    if dtype == 'azav':
        these_vals = this_di['vals']
        vals += this_count*these_vals[:, :, this_lut[qv]]
    elif dtype == 'gav':
        these_vals = this_di['vals']
        vals += this_count*these_vals[this_lut[qv]]
    elif dtype == 'shav':
        these_vals = this_di['vals']
        vals += this_count*these_vals[:, this_lut[qv]]
    elif dtype == 'merav':
        these_vals = this_di['vals']
        vals += this_count*these_vals[:, :, :, this_lut[qv]]
    elif dtype == 'specav':
        this_fullpower = this_di['fullpower']
        this_lpower = this_di['lpower']
        fullpower += this_count*this_fullpower[:, :, :, this_lut[qv]]
        lpower += this_count*this_lpower[:, :, this_lut[qv], :]
    elif dtype == 'enstr':
        these_vals = this_di['vals']
        vals += this_count*these_vals

    count += this_count
    print ('total count = %i ...' %count)
    print ('-------------------------------')

if dtype == 'specav':
    fullpower /= count
    lpower /= count
else:
    vals /= count

# Initialize joined dictionary, then change it
di = dict(di_list[0])

if dtype == 'specav':
    di['fullpower'] = fullpower
    di['lpower'] = lpower
else:
    di['vals'] = vals

iter1, iter2 = di_list[0]['iter1'], di_list[nfiles - 1]['iter2']
di['iter1'] = iter1
di['iter2'] = iter2
di['count'] = count

di['qv'] = qv
di['nq'] = nq
di['lut'] = lut

basename = '_AZ_Avgs_'
if dtype == 'gav':
    basename = '_G_Avgs_'
elif dtype == 'shav':
    basename = '_Shell_Avgs_'
elif dtype == 'specav':
    basename = '_Shell_Spectra_'
elif dtype == 'merav':
    basename = '_Meridional_Slices_'
elif dtype == 'enstr':
    basename = '_enstrophy_from_mer_'

savename = dirname_stripped + basename + tag + str(iter1).zfill(8) +\
        '_' + str(iter2).zfill(8) + '.pkl'
savefile = datadir + savename
print ("Saving joined average data in %s ..." %savefile)
f = open(savefile, 'wb')
pickle.dump(di, f, protocol=4)
f.close()
