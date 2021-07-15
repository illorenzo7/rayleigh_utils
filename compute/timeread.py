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
kwargs_default = dict({'radtype': 'azav', 'irvals': np.array([0]), 'qvals': np.array([1])})

# overwrite defaults
kw = update_dict(kwargs_default, clas)
radtype = kw.radtype

lent = 50
char = '.'
t1_glob = time.time()
t1 = t1_glob + 0.0

if radtype == 'azav':
    reading_func = AZ_Avgs
    dataname = 'AZ_Avgs'
if radtype == 'shav':
    reading_func = Shell_Avgs
    dataname = 'Shell_Avgs'
if radtype == 'gav':
    reading_func = G_Avgs
    dataname = 'G_Avgs'
if radtype == 'spec':
    reading_func = Shell_Spectra
    dataname = 'Shell_Spectra'
if radtype == 'specalt':
    reading_func = rdalt.Shell_Spectra
    dataname = 'Shell_Spectra'
if radtype == 'sslice':
    reading_func = Shell_Slices
    dataname = 'Shell_Slices'
if radtype == 'merslice':
    reading_func = Meridional_Slices
    dataname = 'Meridional_Slices'
if radtype == 'eqslice':
    reading_func = Equatorial_Slices
    dataname = 'Equatorial_Slices'

# Get the Rayleigh data directory
radatadir = dirname + '/' + dataname + '/'

# Get desired file names in radatadir and their integer counterparts
file_list, int_file_list, nfiles = get_file_lists(radatadir, args)

print (buff_line)
print ('Considering %i %s files: %s through %s'\
    %(nfiles, dataname, file_list[0], file_list[-1]))
if radtype == 'specalt':
    print ("irvals = ", kw.irvals)
    print ("qvals = ", kw.qvals)
print (buff_line)

# now time the reading
total_read = 0
for i in range(nfiles):
    t1 = time.time()
    fname = radatadir + str(file_list[i]).zfill(8)
    print (fill_str("reading " + fname, lent, char), end='')
    if radtype == 'specalt':
        a = reading_func(fname, '', irvals=kw.irvals, qvals=kw.qvals)
    else:
        a = reading_func(fname, '')
    t2 = time.time()
    size_in_M = sys.getsizeof(a.vals)/1024**2
    total_read += size_in_M
    print ('%.1f M' %size_in_M + 3*' ', end='')
    print (format_time(t2 - t1), end='')
    io_rate = size_in_M/(t2 - t1)
    print (5*' ' + '%.1f M/s' %io_rate)

t2 = time.time()
print (fill_str("total time", lent, char) + format_time(t2 - t1_glob))
print (fill_str("avg. time", lent, char) + format_time((t2 - t1_glob)/nfiles))
print (fill_str("avg. I/O speed", lent, char) + "%.1f M/s" %(total_read/(t2 - t1_glob)))
