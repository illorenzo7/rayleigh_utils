# Author: Loren Matilsky
# Created: 05/18/2019
# This script computes the volume-averaged Rossby number for a 
# Rayleigh run in directory [dirname], using the (constant) shell depth
# as the length scale. 
# Takes data from Shell_Avgs file and the rotation rate from the 
# main input file
# Reads grid_info for the radial weights
# Displays the computed Rossby number at the terminal and writes it
# as an empty file in the directory [dirname]

import numpy as np
import sys, os
sys.path.append(os.environ['rapp'])
sys.path.append(os.environ['co'])
from get_parameter import get_parameter
from rayleigh_diagnostics import Shell_Avgs, GridInfo
from common import get_widest_range_file, get_dict

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
Om = get_parameter(dirname, 'angular_velocity')

# Find the rms convective velocity
vsq_r, vsq_t, vsq_p = vals[:, lut[422]], vals[:, lut[423]],\
    vals[:, lut[424]], 
vsq = vsq_r + vsq_t + vsq_p

# Compute volume-average of the velocity squared
# using the radial integration weights
vsq_av = np.sum(vsq*rw)

# Compute the Rossby number
Ro = np.sqrt(vsq_av)/(2.0*Om*H)

# And print it
print("The velocity Rossby number (length scale = shell depth) is %1.3e"\
        %Ro)

# Also write it as an empty file in [dirname]
# Remove all other "Ro_vel1_is_" that might be present
# (in case this definition of Ra is replaced
names = os.listdir(dirname)
for name in names:
    if "Ro_vel1_is_" in name:
        os.remove(dirname + '/' + name)
fname = dirname + ("/00_Ro_vel1_is_%1.3e" %Ro)
f = open(fname, "w")
f.close()
