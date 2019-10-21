# Author: Loren Matilsky
# Created: 05/01/2019
# Modified to work with new cref framework: 10/20/2019
#
# Purpose: generate a binary file (for Rayleigh to read) that contains
# a reference state, consisting of a stable region (stiffness k) underlying 
# a neutrally stable CZ

# Parameters: output_dir (first argument), 
# Inner, outer, and transition radii (default ri = 3.4139791e10, 
# rm = 5e10, ro = 6.5860209e10...for CZ corresponding to 3 density
# scale heights )
# nr1 (number of Chebyshevs to use for CZ, -nr1) default 96
# nr2 (number of Chebyshevs to use for RZ, -nr2) default 64
# Density at transition boundary (-rhom) default 0.18053428
# Temperature at transition boundary (-tm) default 2.111256e6
# Stiffness k (-k) default 2
# Transition width delta (-delta) (as a fraction of rm) default 0.005
# gamma (-gam) default 5/3
# pressure specific heat cp (-cp) default 3.5e8
# central mass (-M) default M_sun

import numpy as np
import sys, os
from arbitrary_atmosphere import arbitrary_atmosphere

import basic_constants as bc

sys.path.append(os.environ['rapp'])
from reference_tools import equation_coefficients

# Set default constants
ri = 4.176e10  # Set RZ width about 0.5x CZ width
rm = bc.rm
ro = bc.ro
cp = bc.cp

Tm = bc.Tm
pm = bc.pm
rhom = bc.rhom
gam = bc.gamma
k = 2.0
delta = 0.005*ro

# Get directory to save binary files for reference state and heating
dirname = sys.argv[1]

args = sys.argv[2:]
nargs = len(args)
for i in range(nargs):
    arg = args[i]
    if arg == '-ri':
        ri = float(args[i+1])
    elif arg == '-rm':
        rm = float(args[i+1])
    elif arg == '-ro':
        ro = float(args[i+1])
    elif arg == '-nr':
        nr = int(args[i+1])
    elif arg == '-rhom':
        rhom = float(args[i+1])        
    elif arg == '-tm':
        tm = float(args[i+1])   
    elif arg == '-k':
        k = float(args[i+1])  
    elif arg == '-delta':
        delta = float(args[i+1])*ro
    elif arg == '-gam':
        gam = float(args[i+1]) 
    elif arg == '-cp':
        cp = float(args[i+1]) 
    elif arg == '-M':
        M = float(args[i+1])  
        
# First, compute reference state on super-fine grid to interpolate onto later    
nr = 5000
rr = np.linspace(ro, ri, nr)

d2sdr2 = -k*cp/rm*(1./2./delta)*(1./np.cosh((rr - rm)/delta))**2.
dsdr = k*cp/rm*0.5*(1.0 - np.tanh((rr - rm)/delta))
s = k*cp*0.5*((rr/rm - 1.0) -\
        (delta/rm)*np.log(np.cosh((rr - rm)/delta)))

g = bc.G*bc.M/rr**2
dgdr = -2.0*g/rr

T, rho, p, dlnT, dlnrho, dlnp, d2lnrho =\
    arbitrary_atmosphere(rr, s, dsdr, d2sdr2, g,\
                         dgdr, rm, Tm, pm, cp, gam)

print("Computed atmosphere for RZ-CZ, ds/dr joined with tanh")

# Now write to file using the equation_coefficients framework
eq = equation_coefficients(rr)

# Set only the thermodynamic functions/constants in this routine
# In other routines, we can set the heating and transport coefficients

print("Setting f_1, f_2, f_4, f_8, f_9, f_10, and f_14")
eq.set_function(rho, 1)
buoy = rho*g/cp
eq.set_function(buoy, 2)
eq.set_function(T, 4)
eq.set_function(dlnrho, 8)
eq.set_function(d2lnrho, 9)
eq.set_function(dlnT, 10)
eq.set_function(dsdr, 14)

print("Setting c_2, c_3, c_4, c_5, c_6, c_7, c_8, and c_9")
eq.set_constant(1.0, 2) # multiplies buoyancy
eq.set_constant(1.0, 3) # multiplies pressure grad.
eq.set_constant(1.0/4.0/np.pi, 4) # multiplies Lorentz force
eq.set_constant(1.0, 5) # multiplies viscous force
eq.set_constant(1.0, 6) # multiplies thermal diffusion term
eq.set_constant(1.0, 7) # multiplies eta in induction equation
eq.set_constant(1.0, 8) # multiplies viscous heating term
eq.set_constant(1.0/4.0/np.pi, 9) # multiplies magnetic diffusion term

# Will need to figure out how to deal with c_1 (supposed to be 2 x angular velocity, i.e., the Coriolis coefficient. Hopefully we don't need c_1 in the
# custom reference framework and will just specify angular_velocity
# If this doesn't work, will need to use override_constants framework

# c_10 will be set in the "generate heating" scripts

# The "generate transport" scripts will only set the profiles, not
# the constants c_5, c_6, c_7, c_8, c_9

the_file = dirname + '/custom_reference_binary'

print("Writing the atmosphere to %s" %the_file)
eq.write(the_file)
