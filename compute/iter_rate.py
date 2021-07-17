# Author: Loren Matilsky
# Date created: 05/14/2020
# Prints iteration rate for directory [dirname] (argument #1),
# namely iters_per_second (in real time) and
# seconds_per_iter (in sim. time)
# By default uses the last of the numbered logfiles in the directory
# or else the logfile specified by -fname
import numpy as np
import os, sys
sys.path.append(os.environ['rapp'])
sys.path.append(os.environ['raco'])
from common import *

# Get the run directory on which to perform the analysis
dirname = sys.argv[1]
args = sys.argv[2:]
nargs = len(args)

# See if user wants to use a different file name than last logfile
fname = None
verbose = False
for i in range(nargs):
    arg = args[i]
    if arg == '--fname':
        fname = args[i+1]
    elif arg == '--v':
        verbose = True

lognames = []
lognumbers = []
for name in os.listdir(dirname):
    if 'logfile' in name:
        # Make sure there is number following logfile
        lastpart = name[7:]
        try:
            lognumbers.append(int(lastpart))
            lognames.append(name)
        except:
            print ("Can't discern a number for %s" %name)
            print ("Not considering %s by default" %name)
            print ("To use it, specify -fname %s" %name)

# convert to arrays
lognames = np.array(lognames)
lognumbers = np.array(lognumbers)

if fname is None:
    imax = np.argmax(lognumbers)
    fname = lognames[imax]

print ("In %s:" %fname)
di = read_log(dirname + '/' + fname)
print ("ncpu: ", di['ncpu'])
iters = di['iters']
niter = iters[-1] - iters[0]
iters_per_sec = di['iters_per_sec'][1:] # first one is zero for some reason
secs_per_iter = 1/iters_per_sec
runtime = np.sum(secs_per_iter)
iters_per_sec_av = niter/runtime
dt_av = np.mean(di['delta_t'])
print ("avg. iters/sec: %.2f" %iters_per_sec_av)
print ("mean of rates:  %.2f" %np.mean(iters_per_sec))
print ("avg. timestep:  %1.2e sec" %dt_av)
print ("run time     :  " + format_time(runtime))
print ("niter        :  ", niter)
print ("len(iters)   :  ", len(iters_per_sec))

if verbose:
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
