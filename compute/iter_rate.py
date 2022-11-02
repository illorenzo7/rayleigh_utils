# Author: Loren Matilsky
# Date created: 05/14/2020
# Prints iteration rate for directory [dirname] (argument #1),
# namely iters_per_second (in real time) and
# seconds_per_iter (in sim. time)
# By default uses the last of the numbered logfiles in the directory
# or else the logfile specified by --fname
import numpy as np
import os, sys
sys.path.append(os.environ['raco'])
from cla_util import *
from common import *

# read in args
clas0, clas = read_clas(sys.argv)
dirname = clas0.dirname
rotation = clas0.rotation

if clas.fname is None: # default
    lognames = []
    for name in os.listdir(dirname):
        if 'logfile' in name:
            lognames.append(name)
    lognames = np.sort(lognames)
    if len(lognames) > 0:
        fname = lognames[-1]
    else:
        print ('no logfiles found')
        sys.exit()
else:
    fname = clas.fname

fullname = dirname + '/' + fname
print (buff_line)
print ("In %s:" %fullname)
di = read_log(fullname)
iters = di['iters'][1:]
iters_per_sec = di['iters_per_sec'][1:] # first one is zero for some reason
delta_t = di['delta_t'][1:]

# treat the iterations just after output separately
iters_io = di['iters_io']
iters_per_sec_io = di['iters_per_sec_io'] 
delta_t_io = di['delta_t_io']

# now adjust this by desired range
it1, it2 = 0, len(iters) # default look at full range
args = sys.argv
nargs = len(args)
for i in range(len(sys.argv)):
    arg = sys.argv[i]
    if arg == '--f':
        it2 = int(float(args[i+1]))
    if arg == '--n':
        it1 = it2 - int(float(args[i+1]))
    if arg == '--range':
        it1 = np.argmin(np.abs(iters - float(args[i+1])))
        it2 = np.argmin(np.abs(iters - float(args[i+2])))
    if arg == '--centerrange':
        itc = np.argmin(np.abs(iters - float(args[i+1])))
        length = int(float(args[i+1]))
        it1 = itc - length//2
        it2 = it1 + length

# safety check
if it1 < 0:
    it1 = 0
if it2 > len(iters):
    it2 = len(iters)

iter1, iter2 = iters[it1], iters[it2 - 1]
it1_io = np.argmin(np.abs(iters_io - iter1))
it2_io = np.argmin(np.abs(iters_io - iter2))

# shorten the arrays possibly (by default, no)
iters = iters[it1:it2]
iters_per_sec = iters_per_sec[it1:it2]
delta_t = delta_t[it1:it2]

iters_io = iters_io[it1_io:it2_io]
iters_per_sec_io = iters_per_sec_io[it1_io:it2_io]
delta_t_io = delta_t_io[it1_io:it2_io]

print ("range = " + str(iter1).zfill(8) + ' to ' + str(iter2).zfill(8))
print ("ncpu: ", di['ncpu'])

dt_av = np.mean(delta_t)
niter = len(iters)
secs_per_iter = 1/iters_per_sec
secs_per_iter_io = 1/iters_per_sec_io
runtime = np.sum(secs_per_iter)
iotime = np.sum(secs_per_iter_io)

iters_per_sec_av = niter/runtime
print (make_bold("avg. iters/sec: %.2f" %iters_per_sec_av))
print ("mean of rates:  %.2f" %np.mean(iters_per_sec))
print (make_bold("avg. timestep:  %1.2e sec" %dt_av))
print ("run time     :  " + format_time(runtime))
print ("est. I/O time:  " + format_time(iotime))
print ("frac run     :    %.3f" %((runtime - iotime)/runtime))
print ("frac I/O     :    %.3f" %(iotime/runtime))
print ("niter        :  ", niter)

if clas.verbose:
    # Get fancy now ...
    print ("===============================")
    # Get the baseline time unit
    time_unit, time_label, rotation, simple_label = get_time_unit(dirname)
    unit_name = simple_label

    simtime_per_hour = dt_av*(iters_per_sec_av*3600.)/\
            time_unit
    print ("Simulation rate = %.1f %s/hour" %(simtime_per_hour, unit_name))
    print ("Simulation rate = %.1f %s/day" %(simtime_per_hour*12.,\
            unit_name))
    print ("Simulation rate = %.1f %s/(5 days)" %(simtime_per_hour*12.*5.,\
            unit_name))
    print ("===============================")
    print ("Simulation rate = %1.2e iters/hour" %(iters_per_sec_av*3600.))
    print ("Simulation rate = %1.2e iters/day"\
            %(iters_per_sec_av*3600.*24.))
    print ("Simulation rate = %1.2e iters/(5 days)"\
            %(iters_per_sec_av*3600.*24.*5.))
    print ("===============================")
    print ("Min. time step = %1.2e s" %np.min(di['delta_t']))
    print ("Max. time step = %1.2e s" %np.max(di['delta_t']))
