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

# Find the relevant place to store the data, and create the directory if it
# doesn't already exist

# Make opportunity for command-line args
files = []
dirname = sys.argv[1]
datadir = dirname + '/data/'
dirname_stripped = strip_dirname(dirname)

tag = ''
# Other choices are gav, shav, specav, merav, enstr
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
print ('nfiles =', nfiles)
# Read in all the dictionaries to be conjoined
di_list = []
for i in range(nfiles):
    di = get_dict(files[i])
    di_list.append(di)
di0 = di_list[0]

vals = np.zeros_like(di0['vals'])

# calculate the total count for weights
count = 0
for i in range(nfiles):
    count += di_list[i]['count']

for i in range(nfiles):
    print('adding %s' %files[i])
    this_di = di_list[i]
    this_count = this_di['count']
    print('weight = %.3f' %(this_count/count))
    these_vals = this_di['vals']
    vals += this_count*these_vals
vals /= count

# Initialize joined dictionary, then change it
di = dict(di_list[0])
di['vals'] = vals

iter1, iter2 = di_list[0]['iter1'], di_list[nfiles - 1]['iter2']
di['iter1'] = iter1
di['iter2'] = iter2
di['count'] = count

basename = '_mag_torque_from_mer_'

savename = dirname_stripped + basename + tag + str(iter1).zfill(8) +\
        '_' + str(iter2).zfill(8) + '.pkl'
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
