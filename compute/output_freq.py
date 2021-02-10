# Author: Loren Matilsky
# Created: 02/10/2021
# script to compute the average output frequency for a given run directory
# (first argument, as relative path to current working directory)

from common import *

dirname = sys.argv[1]
radatadir = dirname + '/' + sys.argv[2]
datadir = dirname + '/data/'

print ("In %s:" %radatadir)

# Get the baseline time unit
rotation = get_parameter(dirname, 'rotation')
if rotation:
    time_unit = compute_Prot(dirname)
    time_label = 'rotations'
else:
    time_unit = compute_tdt(dirname)
    time_label = 'TDTs'
print(time_label)
print(datadir)
the_file = get_widest_range_file(datadir, 'trace_G_Avgs')
if the_file == '':
    the_file = get_widest_range_file(datadir, 'trace_2dom_G_Avgs')
print(the_file)
di = get_dict(datadir + the_file) 
times = di['times']
iters = di['iters']
print ("got time/iter translation from %s" %the_file)

dummy, int_file_list, nfiles = get_file_lists(radatadir)
times_data = np.zeros(nfiles)
for i in range(nfiles):
    iter_loc = int_file_list[i]
    j = np.argmin(np.abs(iters - iter_loc))
    times_data[i] = times[j]/time_unit
times_diff = np.diff(times_data)
iters_diff = np.diff(int_file_list)
lent = 30
char = ' '
print (fill_str('avg cadence', lent, char) +\
        make_bold(('%.3f '  + time_label) %np.mean(times_diff)))
print (fill_str('stddev cadence', lent, char) +\
        ('%.3f '  + time_label) %np.std(times_diff))
print (fill_str('avg output frequency', lent, char) +\
        ('%.1f iters') %np.mean(iters_diff))

