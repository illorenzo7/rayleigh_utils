# Author: Loren Matilsky
# Created: 08/19/2019
# This script computes the volume-averaged Reynolds number for a 
# Rayleigh run in directory [dirname], using the (constant) shell depth
# as the length scale. 
# Gets diffusion profiles from transport or equation_coefficients
# Reads grid_info for the radial weights
# Displays the computed Reynolds number at the terminal

import numpy as np
import sys, os
sys.path.append(os.environ['rapp'])
sys.path.append(os.environ['raco'])
from rayleigh_diagnostics import Shell_Avgs, GridInfo
from common import *

# Get directory name
dirname = sys.argv[1]

# Read in the Shell_Avgs data
datadir = dirname + '/data/'
Shell_Avgs_file = get_widest_range_file(datadir, 'Shell_Avgs')
print ('Getting velocity amplitudes from ' + datadir +\
        Shell_Avgs_file + ' ...')
di = get_dict(datadir + Shell_Avgs_file)
vals = di['vals']
lut = di['lut']
rr = di['rr']
H = np.max(rr) - np.min(rr)

# Read in grid info for radial weights and reference velocity
gi = GridInfo(dirname + '/grid_info')
rw = gi.rweights

# Read in transport coefficients for nu-profile
eq = get_eq(dirname)
nu = eq.nu

# Find the rms convective velocity, averaged over spheres
vsq_r, vsq_t, vsq_p = vals[:, lut[422]], vals[:, lut[423]],\
    vals[:, lut[424]], 
vsq = vsq_r + vsq_t + vsq_p
v_rms = np.sqrt(vsq)

# Compute volume-averaged Reynolds number
# using the radial integration weights
Re_vs_r = v_rms*H/nu
Re = np.sum(rw*Re_vs_r)

# And print it
print("The volume-averaged Reynolds number (length scale = shell depth) is %1.3e"\
        %Re)
