# Author: Loren Matilsky
# Created: 09/13/2019
# This script computes the "volume-averaged" magnetic Reynolds number for a 
# Rayleigh run in directory [dirname], using a radially dependent spectral
# rms l-value for the length scale (based on the magnetic field). 
# Gets diffusion profiles from transport or equation_coefficients
# Displays the computed Reynolds number at the terminal

import numpy as np
import sys, os
sys.path.append(os.environ['rapp'])
sys.path.append(os.environ['raco'])
from get_parameter import get_parameter
from rayleigh_diagnostics import Shell_Avgs, TransportCoeffs
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

# Read in transport coefficients for nu-profile
t = TransportCoeffs(dirname + '/transport')
eta = t.eta

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
brsq_power = lpower[:, :, lut[801], 2] # get the convective B-field power
btsq_power = lpower[:, :, lut[802], 2] 
bpsq_power = lpower[:, :, lut[803], 2] 
bsq_power = brsq_power + btsq_power + bpsq_power
l_rms = np.sum(bsq_power*lvals.reshape((nell, 1)), axis=0)/\
        np.sum(bsq_power, axis=0)
twopir = 2*np.pi*rvals_spec
H_spec = twopir/l_rms

# Compute "volume-averaged" Reynolds number (averaging only over the discrete
# radial locations for the Shell_Spectra

Rm = np.mean(v_rms[rinds]*H_spec/eta[rinds])

# And print it
print("The magnetic Reynolds number (length scale from Shell_Spectra) is %1.3e" %Rm)
