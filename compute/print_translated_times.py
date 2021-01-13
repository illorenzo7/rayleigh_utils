# Author: Loren Matilsky
# Date created: well before 05/06/2019
# Prints various timescales associated with a given Rayleigh simulation
# directory, including how long it was run, the thermal diffusion time, 
# rotation period, etc.

import sys
import numpy as np

# Get the run directory on which to perform the analysis
dirname = sys.argv[1]
args = sys.argv[2:]
nargs = len(args)

# See if user wants to work in frame where a time (in Prot or TDT)
# has been subtracted from time
sub = False
for i in range(nargs):
    arg = args[i]
    if arg == '-sub':
        sub = True
        start_time = float(args[i+1])

for i in range(nargs):
    arg = args[i]
    if arg == '-sec': 
        time = float(args[i+1])
        di = translate_times(time, dirname, translate_from='sec')
        val_sec = di['val_sec']; val_iter = di['val_iter']
        val_day = di['val_day']; val_unit = di['val_unit']
        time_unit = di['time_unit']; time_label = di['time_label']
    elif arg == '-iter': 
        time = int(args[i+1])
        di = translate_times(time, dirname, translate_from='iter')
        val_sec = di['val_sec']; val_iter = di['val_iter']
        val_day = di['val_day']; val_unit = di['val_unit']
        time_unit = di['time_unit']; time_label = di['time_label']
    elif arg == '-day':
        time = float(args[i+1])
        di = translate_times(time, dirname, translate_from='day')
        val_sec = di['val_sec']; val_iter = di['val_iter']
        val_day = di['val_day']; val_unit = di['val_unit']
        time_unit = di['time_unit']; time_label = di['time_label']
    elif arg in ['-prot', '-tdt', '-unit']:
        time = float(args[i+1])
        if sub:
            time += start_time
        di = translate_times(time, dirname, translate_from='unit')
        val_sec = di['val_sec']; val_iter = di['val_iter']
        val_day = di['val_day']; val_unit = di['val_unit']
        time_unit = di['time_unit']; time_label = di['time_label']

if sub:
    val_sec -= (start_time*time_unit)
    val_day -= (start_time*time_unit/86400.)
    val_unit -= start_time

print ('%1.6e sec   = %08i iter   = %.2f day    = %.3f %s'\
        %(val_sec, val_iter, val_day, val_unit, time_label))
