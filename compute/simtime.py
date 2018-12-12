import numpy as np
import os
import sys
sys.path.append('~/00_Rayleigh-clone/post-processing/')
from rayleigh_diagnostics import G_Avgs

dirname = sys.argv[1]
dirname_stripped = (dirname.split('/'))[-1]
gavg = dirname + '/G_Avgs/'

files = os.listdir(gavg)
files.sort()
nfiles = len(files)

args = sys.argv[2:]
nargs = len(args)
for i in range(nargs):
    arg = args[i]
    if arg == '-n':
        start_iter = int(args[i+1])

f1 = files[0]
f2 = files[nfiles - 1]

a1 = G_Avgs(gavg + f1, '')
a2 = G_Avgs(gavg + f2, '')

simtime_sec = a2.time[0] - a1.time[0]
simtime_days = simtime_sec/86400.
starttime = a1.time[0]/86400.

print ('Simulation %s started at %.1f days' %(dirname_stripped, starttime))
print('Simulation %s has run for a total of %.1f days'\
        %(dirname_stripped, simtime_days))
print('Or %.1f years' %(simtime_days/365))
