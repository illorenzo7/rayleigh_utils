# Author: Loren Matilsky
# Created: 10/21/2019
#
# Purpose: modify a binary file (default name custom_reference_binary) 
# to contain heating profile that is confined to CZ, transitioning to
# no heating of the RZ -- by default at the stable/unstable layer
# transition
# Must be run AFTER reference state (which includes the 
# density/temperature) is generated

# Parameters: output_dir (first argument), 
# command line options:
# -rt
# transition radius for heating
#
# -delta
# Transition width delta (as a fraction of rsun) default 0.005
#
# -lum
# Luminosity to be driven through layer, default Lsun

import numpy as np
import sys, os
from scipy.integrate import simps

sys.path.append(os.environ['rapp'])
sys.path.append(os.environ['raco'])

import basic_constants as bc

from reference_tools import equation_coefficients
from common import rsun, lsun

# Set default constants
rt = 5.0e10 # by default transition a bit below RZ-CZ transition
delta = 0.005*rsun
lum = lsun

# Get directory to save binary files for reference state and heating
dirname = sys.argv[1]
fname = 'custom_reference_binary'

args = sys.argv[2:]
nargs = len(args)
for i in range(nargs):
    arg = args[i]
    if arg == '-rt':
        rt = float(args[i+1])
    elif arg == '-delta':
        delta = float(args[i+1])*rsun
    elif arg == '-fname':
        fname = args[i+1]
    elif arg == '-lum':
        lum = float(args[i+1])
        
# Open and read the hopefully already existing reference file!
eq = equation_coefficients()
the_file = dirname + '/' + fname
eq.read(the_file)
rr = eq.radius
rho = eq.functions[0]
T = eq.functions[3]
smooth = 0.5*(1.0 + np.tanh((rr - rt)/delta)) # "detects" CZ
profile = (rho*T - rho[0]*T[0])*smooth

int_profile = -4*np.pi*simps(profile*rr**2.0, rr) # remember rr is reversed

radial_shape = profile/int_profile

print("---------------------------------")
print("Computed radial shape of heating for RZ-CZ, joined with tanh")
print("rt: %1.3e cm" %rt) 
print("delta/rsun: %.3f"  %(delta/rsun))
print("---------------------------------")

# Now write to file using the equation_coefficients framework
print("Setting f_6 and c_10")

eq.set_function(radial_shape, 6)
eq.set_constant(lum, 10) # Total luminosity normalizes heating function

# Will need to figure out how to deal with c_1 (supposed to be 2 x angular velocity, i.e., the Coriolis coefficient. Hopefully we don't need c_1 in the
# custom reference framework and will just specify angular_velocity
# If this doesn't work, will need to use override_constants framework

print("Writing the heating to %s" %the_file)
print("---------------------------------")
eq.write(the_file)
