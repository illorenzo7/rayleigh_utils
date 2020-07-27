# Author: Loren Matilsky
# Created: 05/01/2019
# Modified to work with new cref framework: 10/20/2019
#
# Purpose: generate a binary file (for Rayleigh to read) that contains
# a reference state, consisting of a stable region (stiffness k) underlying 
# a neutrally stable CZ

# Parameters: output_dir (first argument), 
# Command-line options:
#
# -ri, -rm, -ro
# Inner, outer, and transition radii (default ri = 3.4139791e10, 
# rm = 5e10, ro = 6.5860209e10...for CZ corresponding to 3 density
# scale heights )
#
# -rhom
# Density at transition boundary default 0.18053428
#
# -tm
# Temperature at transition boundary default 2.111256e6
#
# -k
# Stiffness k, default 2
# 
# -delta
# Transition width delta as a fraction of rsun, default 0.010
#
# -gam
# specific heat ratio gamma, default 5/3
#
# -cp
# pressure specific heat cp, default 3.5e8
#
# -M
# central mass M, default M_sun
#
# -mag
# Whether magnetism is True or False, default False (hydro)

import numpy as np
import sys, os
from arbitrary_atmosphere import arbitrary_atmosphere

sys.path.append(os.environ['rapp'])
sys.path.append(os.environ['raco'])
sys.path.append(os.environ['raco'] + '/tachocline')

from reference_tools import equation_coefficients
from common import rsun
import basic_constants as bc

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
delta = 0.010*rsun
mag = False

# Get directory to save binary files for reference
dirname = sys.argv[1]
fname = 'custom_reference_binary'

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
    elif arg == '-rhom':
        rhom = float(args[i+1])        
    elif arg == '-tm':
        tm = float(args[i+1])   
    elif arg == '-k':
        k = float(args[i+1])  
    elif arg == '-delta':
        delta = float(args[i+1])*rsun
    elif arg == '-gam':
        gam = float(args[i+1]) 
    elif arg == '-cp':
        cp = float(args[i+1]) 
    elif arg == '-M':
        M = float(args[i+1])  
    elif arg == '-fname':
        fname = args[i+1]
    elif arg == '-mag':
        mag = True

# First, compute reference state on super-fine grid to interpolate onto later    
nr = 5000
rr = np.linspace(ro, ri, nr)

# Define an entropy profile that is +1 for r < rm - delta, 0 for r > rm,
# and continuously differentiable (a quartic) in between

s = np.zeros(nr)
dsdr = np.zeros(nr)
d2sdr2 = np.zeros(nr)

for i in range(nr):
    rloc = rr[i]
    if rloc <= rm - delta:
        s[i] = (8.0/15.0)*(k*cp)*(delta/rm) + k*cp*(rloc/rm - 1.0)
        dsdr[i] = k*cp/rm
        d2sdr2[i] = 0.0
    elif rloc > rm - delta and rloc < rm:
        x = (rloc - rm)/delta
        s[i] = (k*cp)*(delta/rm)*((2.0/3.0)*x**3.0 - (1.0/5.0)*x**5.0)
        dsdr[i] = (k*cp/rm)*(1.0 - (1.0 - x**2.0)**2.0)
        d2sdr2[i] = (4.0/delta)*(k*cp/rm)*(1.0 - x**2.0)*x
    else:
        s[i] = 0.0
        dsdr[i] = 0.0
        d2sdr2[i] = 0.0

# Make gravity due to a point mass at the origin
g = bc.G*bc.M/rr**2
dgdr = -2.0*g/rr

T, rho, p, dlnT, dlnrho, dlnp, d2lnrho =\
    arbitrary_atmosphere(rr, s, dsdr, d2sdr2, g,\
                         dgdr, rm, Tm, pm, cp, gam)

print("---------------------------------")
print("Computed atmosphere for RZ-CZ, ds/dr joined with a quartic")
print("ri: %1.3e cm" %ri) 
print("rm: %1.3e cm" %rm) 
print("ro: %1.3e cm" %ro) 
print("delta/rsun: %.3f"  %(delta/rsun))
print("k [(dsdr in RZ)/cp*rm]: %1.1e"  %k)
print("---------------------------------")

# Now write to file using the equation_coefficients framework
eq = equation_coefficients(rr)

# Set only the thermodynamic functions/constants in this routine
# In other routines, we can set the heating and transport coefficients
# Only set c_4 = 1/(4*pi) if mag = True

print("Setting f_1, f_2, f_4, f_8, f_9, f_10, and f_14")
eq.set_function(rho, 1)
buoy = rho*g/cp
eq.set_function(buoy, 2)
eq.set_function(T, 4)
eq.set_function(dlnrho, 8)
eq.set_function(d2lnrho, 9)
eq.set_function(dlnT, 10)
eq.set_function(dsdr, 14)

print("Setting c_2, c_3, c_7, and c_8")
eq.set_constant(1.0, 2) # multiplies buoyancy
eq.set_constant(1.0, 3) # multiplies pressure grad.
eq.set_constant(1.0, 8) # multiplies viscous heating

if mag:
    print("magnetism = True, so setting c_4 = 1/(4*pi), c_7 = c_9 = 1")
    eq.set_constant(1.0/4.0/np.pi, 4) # multiplies Lorentz force
    eq.set_constant(1.0, 7) # multiplies eta in induction-diffusion term
    eq.set_constant(1.0, 9) # multiplies Joule heating
else:
    print("magnetism = False, so setting c_4, c_7, c_9 = 0.0")
    eq.set_constant(0.0, 4) # multiplies Lorentz force
    eq.set_constant(0.0, 7) # multiplies eta in induction-diffusion term
    eq.set_constant(0.0, 9) # multiplies Ohmic heating

# c_10 will be set in the "generate_heating" scripts

# The "generate_transport" scripts will set the transport
# "radial shapes", and the constants c_5, c_6, c_7

the_file = dirname + '/' + fname

print("Writing the atmosphere to %s" %the_file)
print("---------------------------------")
eq.write(the_file)
