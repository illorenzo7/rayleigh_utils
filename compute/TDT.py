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
from time_scales import compute_tdt, compute_Prot
from get_parameter import get_parameter

# Get directory name
dirname = sys.argv[1]

# Specify if there's a tachocline or we want mag. diffusion time
tach = False
mag = False
args = sys.argv[2:]
nargs = len(args)
for i in range(nargs):
    arg = args[i]
    if arg == '-tach':
        tach = True
    elif arg == '-mag':
        mag = True

# Get the diffusion time
if tach:
    TDT, TDT_CZ, TDT_RZ = compute_tdt(dirname, mag=mag, tach=tach)
else:
    TDT = compute_tdt(dirname, mag=mag, tach=tach)

# Get the baseline time unit
rotation = get_parameter(dirname, 'rotation')
if rotation:
    time_unit = compute_Prot(dirname)
    time_string = ' rotations'
else:
    time_unit = 86400.
    time_string = ' days'

if tach:
    print(("TDT across layer: %.1f" %(TDT/time_unit)) + time_string)
    print(("TDT across CZ: %.1f" %(TDT_CZ/time_unit)) + time_string)
    print(("TDT across RZ: %.1f" %(TDT_RZ/time_unit)) + time_string)
else:
    print(("TDT across layer: %.1f" %(TDT/time_unit)) + time_string)
