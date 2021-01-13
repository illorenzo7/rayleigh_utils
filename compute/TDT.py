# Author: Loren Matilsky
# Created: 01/17/2020
# This script computes the volume-averaged Ekman number for a 
# Rayleigh run in directory [dirname], using the (constant) shell depth
# as the length scale. 
# Gets diffusion profiles from transport or equation_coefficients
# Reads grid_info for the radial weights
# Displays the computed Ekman number at the terminal

import numpy as np
import sys, os
sys.path.append(os.environ['raco'])

# Get directory name
dirname = sys.argv[1]

# Specify if there's a tachocline or we want mag. diffusion time
tach = False
mag = False
args = sys.argv[2:]
nargs = len(args)
in_sec = False
for i in range(nargs):
    arg = args[i]
    if arg == '-tach':
        tach = True
    elif arg == '-mag':
        mag = True
    elif arg == '-sec':
        in_sec = True

# Get the diffusion time
if tach:
    TDT, TDT_CZ, TDT_RZ = compute_tdt(dirname, mag=mag, tach=tach)
else:
    TDT = compute_tdt(dirname, mag=mag, tach=tach)

# Get the baseline time unit
rotation = get_parameter(dirname, 'rotation')
if rotation:
    if in_sec:
        time_unit = 1.
        time_string = ' sec'
    else:
        time_unit = compute_Prot(dirname)
        time_string = ' rotations'
else:
    if in_sec:
        time_unit = 1.
        time_string = ' sec'
    else:
        time_unit = 86400.
        time_string = ' days'

if in_sec:
    format_string = "%1.1e"
else:
    format_string = "%.1f"

if tach:
    print(("TDT across layer: " + format_string %(TDT/time_unit)) + time_string)
    print(("TDT across CZ: " + format_string %(TDT_CZ/time_unit)) + time_string)
    print(("TDT across RZ: " + format_string %(TDT_RZ/time_unit)) + time_string)
else:
    print(("TDT across layer: " + format_string %(TDT/time_unit)) + time_string)
