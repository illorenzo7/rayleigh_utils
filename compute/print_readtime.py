##################################################################
# Routine to average Rayleigh data in time (generic)
# Author: Loren Matilsky
# Created: 04/08/2021
##################################################################
# This routine computes how long it takes to read Rayleigh data
# default data type is azav. to specify another 
# (specav, gav, shav, ssav, merav, eqav) use
# --radtype [radataname]
##################################################################

# import bunch of stuff
import time
import numpy as np
import sys, os
sys.path.append(os.environ['raco'])
sys.path.append(os.environ['rapp'])
from rayleigh_diagnostics import AZ_Avgs, Shell_Avgs, G_Avgs,\
        Shell_Spectra, Shell_Slices, Meridional_Slices, Equatorial_Slices
import rayleigh_diagnostics_alt as rdalt

from common import *
from cla_util import *

# args
args = sys.argv
clas0, clas = read_clas(args)
dirname = clas0['dirname']

# set default values for qval and irval
kwargs_default = dict({'radtype': 'azav', 'irvals': None, 'qvals': None})

# overwrite defaults
kw = update_dict(kwargs_default, clas)
radtype = kw.radtype

kw.irvals = make_array(kw.irvals)
kw.qvals = make_array(kw.qvals)

lent = 50
char = '.'
t1_glob = time.time()
t1 = t1_glob + 0.0

# get reading function and dataname from the di_radtypes container
reading_func = di_radtypes[radtype].reading_func
dataname = di_radtypes[radtype].dataname

# Get the Rayleigh data directory
radatadir = dirname + '/' + dataname + '/'

# Get desired file names in radatadir and their integer counterparts
file_list, int_file_list, nfiles = get_file_lists(radatadir, clas)

print (buff_line)
print ('Considering %i %s files: %s through %s'\
    %(nfiles, dataname, file_list[0], file_list[-1]))
if radtype == 'specalt':
    print ("irvals = ", kw.irvals)
    print ("qvals = ", kw.qvals)
print (buff_line)

# now time the reading
total_size = 0
for i in range(nfiles):
    t1 = time.time()
    fname = radatadir + str(file_list[i]).zfill(8)
    print (fill_str("reading " + fname, lent, char), end='')
    if radtype == 'specalt':
        a = reading_func(fname, '', irvals=kw.irvals, qvals=kw.qvals)
    else:
        a = reading_func(fname, '')
    t2 = time.time()
    the_size = sys.getsizeof(a.vals)
    total_size += the_size
    print (format_size(the_size) + 3*' ', end='')
    print (format_time(t2 - t1), end='')
    io_rate = the_size/(t2 - t1)
    print (5*' ' + format_size(io_rate) + '/s')

t2 = time.time()
print (fill_str("total time", lent, char) + format_time(t2 - t1_glob))
print (fill_str("avg. time", lent, char) + format_time((t2 - t1_glob)/nfiles))
print (fill_str("avg. I/O speed", lent, char) + format_size(total_size/(t2 - t1_glob)) + '/s')
