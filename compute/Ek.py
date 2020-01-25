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
sys.path.append(os.environ['rapp'])
sys.path.append(os.environ['raco'])
from get_parameter import get_parameter
from rayleigh_diagnostics import Shell_Avgs, GridInfo,\
        TransportCoeffs
from common import get_widest_range_file, get_dict
from reference_tools import equation_coefficients

# Get directory name
dirname = sys.argv[1]

# Read in grid info for radial weights and reference velocity
gi = GridInfo(dirname + '/grid_info')
rw = gi.rweights
rr = gi.radius
H = np.max(rr) - np.min(rr)

# Read in transport coefficients for nu-profile
try:
    t = TransportCoeffs(dirname + '/transport')
    nu = t.nu
    print ("Got nu(r) from 'transport' file")
except:
    eq = equation_coefficients()
    eq.read(dirname + '/equation_coefficients')
    nu = eq.constants[4]*eq.functions[2]
    print ("Got nu(r) from 'equation_coefficients' file")

# Get angular velocity 
try:
    Om0 = get_parameter(dirname, 'angular_velocity')
    print ("Got Omega_0 from 'main_input' file")
except:
    eq = equation_coefficients()
    eq.read(dirname + '/equation_coefficients')
    Om0 = eq.constants[0]/2.
    print ("Got Omega_0 from 'equation_coefficients' file")

# Compute volume-averaged Ekman number
# using the radial integration weights
Ek_vs_r = nu/(2*H**2*Om0) # note there is technically a 1/sin(theta) in the 
# Definition of Ek, but this averages to 1 over latitude
Ek = np.sum(rw*Ek_vs_r)

# And print it
print("The volume-averaged Ekman number (length scale = shell depth) is %1.3e"\
        %Ek)
