# Author: Loren Matilsky
# Created: 10/20/2020
# This script computes the volume-averaged Reynolds number for a 
# Rayleigh run in directory [dirname], using vorticity (v/om) i
# to get the length scale. 
# Takes data from Shell_Avgs file and the rotation rate from the 
# main input file
# Reads grid_info for the radial weights
# Displays the computed Rossby number at the terminal and writes it
# as an empty file in the directory [dirname]

import numpy as np
import sys, os
sys.path.append(os.environ['rapp'])
sys.path.append(os.environ['raco'])
from rayleigh_diagnostics import Shell_Avgs, GridInfo
from common import get_widest_range_file, get_dict
from get_length_scales import get_length_scales

# Get directory name
dirname = sys.argv[1]

# Read in the Shell_Avgs data
datadir = dirname + '/data/'
Shell_Avgs_file = get_widest_range_file(datadir, 'Shell_Avgs')
print ('Getting velocity amplitudes from ' + datadir +\
        Shell_Avgs_file)
di = get_dict(datadir + Shell_Avgs_file)
vals = di['vals']
lut = di['lut']
rr = di['rr']
ri = np.min(rr)
ro = np.max(rr)

di_len = get_length_scales(dirname)
L_om = di_len['L_om']

# Read in grid info for radial weights and reference velocity
gi = GridInfo(dirname + '/grid_info')
rw = gi.rweights

# Read in transport coefficients for nu-profile
eq = get_eq(dirname)
nu = eq.nu

# Find the rms convective velocity
vsq_r, vsq_t, vsq_p = vals[:, lut[422]], vals[:, lut[423]],\
    vals[:, lut[424]], 
vsq = vsq_r + vsq_t + vsq_p
v_rms = np.sqrt(vsq)

# Compute the radially dependent Reynolds number
Re_vs_r = v_rms*L_om/nu

# Compute volume-average of the Reynolds number
# using the radial integration weights
Re = np.sum(Re_vs_r*rw)

# Print the vorticity length scale (max/min/global avg) over the shell depth
L_om_gav = np.sum(L_om*rw)
print ("L_om_min/(r_o - r_i) = %1.3e" %(np.min(L_om)/(ro - ri)))
print ("L_om_max/(r_o - r_i) = %1.3e" %(np.max(L_om)/(ro - ri)))
print ("L_om_gav/(r_o - r_i) = %1.3e" %(L_om_gav/(ro - ri)))

# Print the volume averaged Reynolds number
print("The vorticity Reynolds number (length scale = v/om) is %1.3e"\
        %Re)
