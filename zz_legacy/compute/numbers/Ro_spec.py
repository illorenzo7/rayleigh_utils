# Author: Loren Matilsky
# Created: 09/13/2019
# This script computes the volume-averaged Rossby number for a 
# Rayleigh run in directory [dirname], using a radially dependent spectral
# rms l-value for the length scale. 
# Gets diffusion profiles from transport or equation_coefficients
# Displays the computed Reynolds number at the terminal

import numpy as np
import sys, os
sys.path.append(os.environ['rapp'])
sys.path.append(os.environ['raco'])
from rayleigh_diagnostics import Shell_Avgs
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

# Get angular velocity
Om0 = get_parameter(dirname, 'angular_velocity')

# Find the rms convective velocity, averaged over spheres
vsq_r, vsq_t, vsq_p = vals[:, lut[422]], vals[:, lut[423]],\
    vals[:, lut[424]], 
vsq = vsq_r + vsq_t + vsq_p
v_rms = np.sqrt(vsq)

# Get length-scales via the spectra
spec_file = get_widest_range_file(datadir, 'Shell_Spectra')
print ('Reading Shell_Spectra data from ' + datadir + spec_file + ' ...')
di_spec = get_dict(datadir + spec_file)
lpower = di_spec['lpower']
lvals = di_spec['lvals']
rinds = di_spec['rinds']
nell = di_spec['nell']
lut = di_spec['lut']
rvals_spec = di_spec['rvals']
vrsq_power = lpower[:, :, lut[1], 2] # get the convective power
vtsq_power = lpower[:, :, lut[2], 2] 
vpsq_power = lpower[:, :, lut[3], 2] 
vsq_power = vrsq_power + vtsq_power + vpsq_power
l_rms = np.sum(vsq_power*lvals.reshape((nell, 1)), axis=0)/\
        np.sum(vsq_power, axis=0)
twopir = 2*np.pi*rvals_spec
H_spec = twopir/l_rms

# Compute "volume-averaged" Reynolds number (averaging only over the discrete
# radial locations for the Shell_Spectra

Ro = np.mean(v_rms[rinds]/(2*H_spec*Om0))

# And print it
print("The volume-averaged Rossby number (length scale from Shell_Spectra) is %1.3e" %Ro)
