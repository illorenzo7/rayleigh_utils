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
kwargs_default = dict({'latvals': default_latvals, 'qvals': None, 'groupname': 'b', 'mval': 1})
kwargs = update_dict(kwargs_default, clas)

if kwargs.qvals is None: # it's a quantity group
    groupname = kwargs.groupname
    qgroup = get_quantity_group(groupname, magnetism)
    qvals = qgroup['qvals']
else:
    qvals = kwargs.qvals
    groupname = input("choose a groupname to save your data: ")

mval = kwargs['mval']

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

basename = 'mertimerad'
    else:
        basename = 'mertimelat'

    # create data directory if it doesn't already exist
    datadir = clas0['datadir'] + basename + '_mval%03i/' %mval
    if not os.path.isdir(datadir):
        os.makedirs(datadir)

    basename += '_' + groupname + rtag
    savename = basename + clas0['tag'] + '-' +\
            file_list[0] + '_' + file_list[-1] + '.pkl'
    savefile = datadir + savename

    # save the data
    f = open(savefile, 'wb')
    di_sav = dict({'vals': vals, 'times': times, 'iters': iters, 'qvals': qvals, 'samplevals': samplevals, 'mval': mval})
    pickle.dump(di_sav, f, protocol=4)
    f.close()
    t2 = time.time()
    print (format_time(t2 - t1))
    print(make_bold(fill_str('total time', lent, char)), end='')
    print (make_bold(format_time(t2 - t1_glob)))
    print ('data saved at ')
    print (make_bold(savefile))
