# Author: Loren Matilsky
# Created: 02/10/2021
# script to compute the average output frequency for a given run directory
# (first argument, as relative path to current working directory)

from common import *
from cla_util import *

# get the CLAs
args = sys.argv
clas0, clas = read_clas(args)
dirname = clas0['dirname']
rotation = clas0['rotation']

# overwrite defaults
kwargs_default = dict({'radtype': 'azav'})
kw = update_dict(kwargs_default, clas)
radtype = kw.radtype

# get the data directory
datadir = dirname + '/data/'
radatadir = dirname + '/' + di_radtypes[radtype].dataname

print ("In %s:" %radatadir)

# Get the baseline time unit
eq = get_eq(dirname)
if rotation:
    time_unit = eq.trot
    time_label = 'rotations'
else:
    time_unit = eq.tdt
    time_label = 'TDTs'

the_file = get_widest_range_file(datadir, 'G_Avgs_trace')
if the_file is None:
    print ("need a valid time trace file for this to work,")
    print ("e.g., G_Avgs_trace")
    print ("Now exiting")
    sys.exit()

di = get_dict(the_file) 
times = di['times']
iters = di['iters']
print ("got time/iter translation from")
print (the_file)

dummy, int_file_list, nfiles = get_file_lists_all(radatadir)
times_data = np.zeros(nfiles)
for i in range(nfiles):
    iter_loc = int_file_list[i]
    j = np.argmin(np.abs(iters - iter_loc))
    times_data[i] = times[j]/time_unit
times_diff = np.diff(times_data)
iters_diff = np.diff(int_file_list)
print (fill_str('avg cadence') +\
        make_bold(('%1.3e '  + time_label) %np.mean(times_diff)))
print (fill_str('stddev cadence') +\
        ('%1.3e '  + time_label) %np.std(times_diff))
print (fill_str('avg output frequency') +\
        ('%1.3e iters') %np.mean(iters_diff))
