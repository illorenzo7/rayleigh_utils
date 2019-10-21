# Author: Loren Matilsky
# Created: 10/21/2019
#
# Purpose: modify a binary file (default name custom_reference_binary) 
# to contain diffusive profiles that descend like a power of density 
# and drop at the RZ-CZ transition
# Must be run AFTER reference state (which includes the density) is 
# generated

# Parameters: output_dir (first argument), 
# command line options:
# -rt
# transition radius for diffusions
#
# -delta
# Transition width delta (as a fraction of rsun) default 0.005
#
# -power
# default -0.5
#
# -drop
# Diffusion drop (default 1000) to plunge diffusivities in radiative zone
#
# -nutop
# default 3 x 10^12 c.g.s
#
# -kappatop
# default 3 x 10^12 c.g.s
#
# -etatop
# default 3 x 10^12 c.g.s

import numpy as np
import sys, os


sys.path.append(os.environ['rapp'])
sys.path.append(os.environ['raco'])

import basic_constants as bc

from reference_tools import equation_coefficients
from common import rsun

# Set default constants
rt = 4.87e10 # by default transition a bit below RZ-CZ transition
delta = 0.005*rsun
power = -0.5
drop = 1000.
nutop = 3.0e12
kappatop = 3.0e12
etatop = 3.0e12

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
    elif arg == '-power':
        power = float(args[i+1])
    elif arg == '-drop':
        drop = float(args[i+1])
    elif arg == '-nutop':
        nutop = float(args[i+1])
    elif arg == '-kappatop':
        kappatop = float(args[i+1])
    elif arg == '-etatop':
        etatop = float(args[i+1])
    elif arg == '-fname':
        fname = args[i+1]
        
# Open and read the hopefully already existing reference file!
eq = equation_coefficients()
the_file = dirname + '/' + fname
eq.read(the_file)
rr = eq.radius
rho = eq.functions[0]
dlnrho = eq.functions[7]
rhotop = rho[0]
monotone = (rho/rhotop)**power
dmonotone = power*monotone*dlnrho

smooth1 = 0.5*(1.0 + np.tanh((rr - rt)/delta)) # "detects" CZ
dsmooth1 = -0.5/delta*1.0/(np.cosh((rr - rt)/delta))**2.0
smooth2 = 0.5*(1.0 - np.tanh((rr - rt)/delta)) # "detects" RZ
dsmooth2 = 0.5/delta*1.0/(np.cosh((rr - rt)/delta))**2.0

# Calculate where diffusion drops to...
radial_shape = smooth1*monotone + smooth2*(1.0/drop)
dradial_shape = dsmooth1*monotone + smooth1*dmonotone + dsmooth2*(1.0/drop)
dlnradial_shape = dradial_shape/radial_shape

print("Computed radial shape for RZ-CZ diffusions, joined with tanh")

# Now write to file using the equation_coefficients framework
# nu, kappa, eta, all get the same radial shapes (unless we want radially
# dependent Prandtl numbers, which would be stupid)
print("Setting f_3, f_5, f_7, f_11, f_12, and f_13")
eq.set_function(radial_shape, 3)
eq.set_function(radial_shape, 5)
eq.set_function(radial_shape, 7)

eq.set_function(dlnradial_shape, 11)
eq.set_function(dlnradial_shape, 12)
eq.set_function(dlnradial_shape, 13)

print("Setting c_5, c_6, c_7, c_8, and c_9")
eq.set_constant(nutop, 5) # multiplies viscous force
eq.set_constant(kappatop, 6) # multiplies thermal diffusion term
eq.set_constant(etatop, 7) # multiplies eta in induction equation
eq.set_constant(nutop, 8) # multiplies viscous heating term
eq.set_constant(etatop/4.0/np.pi, 9) # multiplies magnetic diffusion term

# Will need to figure out how to deal with c_1 (supposed to be 2 x angular velocity, i.e., the Coriolis coefficient. Hopefully we don't need c_1 in the
# custom reference framework and will just specify angular_velocity
# If this doesn't work, will need to use override_constants framework

# c_10 will be set in the "generate heating" scripts

print("Writing the atmosphere to %s" %the_file)
eq.write(the_file)
