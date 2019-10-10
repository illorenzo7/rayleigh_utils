# Author: Loren Matilsky
# Date created: well before 05/06/2019
# Prints various timescales associated with a given Rayleigh simulation
# directory, including how long it was run, the thermal diffusion time, 
# rotation period, etc.

import sys
from translate_times import translate_times
from common import allthrees_start
import numpy as np

# Get the run directory on which to perform the analysis
dirname = sys.argv[1]
args = sys.argv[2:]
nargs = len(args)

# See if user wants to work in frame where allthrees_start has been subtracted
# from time
sub = False
for i in range(nargs):
    arg = args[i]
    if arg == '-sub':
        sub = True

for i in range(nargs):
    arg = args[i]
    if arg == '-sec': 
        time = float(args[i+1])
        di = translate_times(time, dirname, translate_from='sec')
        val_sec = di['val_sec']; val_iter = di['val_iter']
        val_day = di['val_day']; val_prot = di['val_prot']
    elif arg == '-iter': 
        time = int(args[i+1])
        di = translate_times(time, dirname, translate_from='iter')
        val_sec = di['val_sec']; val_iter = di['val_iter']
        val_day = di['val_day']; val_prot = di['val_prot']
    elif arg == '-day':
        time = float(args[i+1])
        di = translate_times(time, dirname, translate_from='day')
        val_sec = di['val_sec']; val_iter = di['val_iter']
        val_day = di['val_day']; val_prot = di['val_prot']
    elif arg == '-prot':
        time = float(args[i+1])
        if sub:
            time += allthrees_start
        di = translate_times(time, dirname, translate_from='prot')
        val_sec = di['val_sec']; val_iter = di['val_iter']
        val_day = di['val_day']; val_prot = di['val_prot']

if sub:
    val_sec -= (allthrees_start*2*np.pi/8.61e-6)
    val_day -= (allthrees_start*2*np.pi/8.61e-6/86400)
    val_prot -= allthrees_start

print ('%1.6e sec   = %08i iter   = %.2f day    = %.3f Prot'\
        %(val_sec, val_iter, val_day, val_prot))
