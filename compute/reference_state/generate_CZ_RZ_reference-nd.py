# Author: Loren Matilsky
# Created: 04/24/2024
#
# Purpose: generate a binary file (for Rayleigh to read) that contains
# a reference state, consisting of a stable region
# i.e., a radiative zone (RZ) adjacent to a 
# a neutrally stable CZ

# assumes g(r) = const x 1/r^2

# Parameters: output_dir (first argument), 

# Command-line options:
#
# --alpha : ratio of RZ width to CZ width
# --beta : ratio of bottom of CZ to top of CZ 
# --delta : transition width between CZ and RZ via quartic matching of dS/dr
# --jup : if "jup" is specified, RZ lies above CZ
# 
# --gamma : ratio of specific heats
# --nrho : number of density scale heights across CZ
# --amp : amplitude of 1/cp (dS/dr) -- yes it's another independent parameter

# --fname
# File to save reference state in (default custom_reference_binary)
#
# --nr
# Default number of radial (evenly spaced) grid points. 
# Default 10,000 (very fine)

import numpy as np
import sys, os
from arbitrary_atmosphere import arbitrary_atmosphere_nd

sys.path.append(os.environ['rapp'])
sys.path.append(os.environ['raco'])

from reference_tools import equation_coefficients
from common import *
from cla_util import *

# Get CLAs
args = sys.argv
clas0, clas = read_clas_raw(args)
dirname = clas0['dirname']

# Set default kwargs
kw_default = dotdict(dict({'alpha': 1., 'beta': 0.759, 'gamma': 1.67, 'delta': 0.219, 'nrho': 3., 'fname': 'custom_reference_binary', 'nr': 10000, 'jup': False, 'amp': 0.453}))
# creates profiles from tachocline cases, Matilsky et al. (2022, 2024)

# overwrite defaults
kw = update_dict(kw_default, clas)

# check for bad keys
find_bad_keys(kw_default, clas, clas0['routinename'], justwarn=True)

# compute geometry of grid
rbcz = kw.beta/(1.-kw.beta)
rtcz = 1./(1.-kw.beta)
if kw.jup: # RZ above CZ
    rt = rbrz = rtcz
    rtrz = rbrz + kw.alpha
    rmin, rmax = rbcz, rtrz
else: # CZ above RZ
    rt = rtrz = rbcz
    rbrz = rtrz - kw.alpha
    rmin, rmax = rbrz, rtcz

# compute reference state on super-fine grid to interpolate onto later    
r = np.linspace(rmax, rmin, kw.nr) # keep radius in decreasing order for consistency with Rayleigh convention

# Define an entropy profile that is +1 for r < rt - delta, 0 for r > rt,
# and continuously differentiable (a quartic) in between

dsdr = np.zeros(kw.nr)

for i in range(kw.nr):
    rloc = rr[i]
    if kw.jup: # RZ is above
        if rloc <= kw.rt:
            dsdr[i] = 0.
        elif rloc < kw.rt + kw.delta and rloc > kw.rt:
            x = (rloc - kw.rt)/kw.delta
            dsdr[i] = 1.0 - (1.0 - x**2.0)**2.0
        else:
            dsdr[i] = 1.0
    else: # CZ is above
        if rloc <= kw.rt - kw.delta:
            dsdr[i] = 1.
        elif rloc > kw.rt - kw.delta and rloc < kw.rt:
            x = (rloc - kw.rt)/kw.delta
            dsdr[i] = 1.0 - (1.0 - x**2.0)**2.0
        else:
            dsdr[i] = 0.0
dsdr *= kw.amp # scale by the non-dimensional amplitude

# compute the atmosphere
rho, tmp, dlnrho, d2lnrho, dlnt, g =\
        arbitrary_atmosphere_nd(r, dsdr, rbcz, rtcz, kw.gamma, kw.nrho)

print(buff_line)
print("Computed atmosphere for RZ-CZ, ds/dr joined with a quartic")
if kw.jup:
    print ("geometry : Jovian (RZ atop CZ)")
else:
    print ("geometry : solar (CZ atop RZ)")
print("nr         : %i" %kw.nr) 
print("alpha      : %1.3f" %kw.alpha)
print("beta       : %1.3f" %kw.beta)
if kw.jup:
    print("   (rbcz, rtcz=rbrz, rtrz): (%1.3f, %1.3f, %1.3f)"\
            %(rbcz,rtcz,rtrz))
else:
    print("   (rbrz, rtrz=rbcz, rtrz): (%1.3f, %1.3f, %1.3f)"\
            %(rbrz,rtrz,rtcz))
print("delta      : %1.3f" %kw.delta)
print("gamma      : %1.3f" %kw.gamma)
print("   n=1/(gamma-1)      : %1.3f" %(1./(kw.gamma-1.)))
print("Nrho       : %1.3f" %kw.nrho)
print("amp        : %1.3f" %kw.amp)
print(buff_line)

# Now write to file using the equation_coefficients framework
eq = equation_coefficients(r)

# Set only the thermodynamic functions/constants in this routine
# In other routines, we can set the heating and transport coefficients
# Only set c_4 = 1/(4*pi) if mag = True

print("Setting f_1, f_2, f_4, f_8, f_9, f_10, and f_14")
eq.set_function(rho, 1)
buoy = rho*g/kw.cp
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

if kw.mag:
    print("magnetism = True, so setting c_4 = c_9 =  1/(4*pi), c_7 = 1")
    eq.set_constant(1.0/4.0/np.pi, 4) # multiplies Lorentz force
    eq.set_constant(1.0, 7) # multiplies eta in induction-diffusion term
    eq.set_constant(1.0/4.0/np.pi, 9) # multiplies Joule heating
else:
    print("magnetism = False, so setting c_4, c_7, c_9 = 0.0")
    eq.set_constant(0.0, 4) # multiplies Lorentz force
    eq.set_constant(0.0, 7) # multiplies eta in induction-diffusion term
    eq.set_constant(0.0, 9) # multiplies Ohmic heating

# c_10 will be set in the "generate_heating" scripts

# The "generate_transport" scripts will set the transport
# "radial shapes", and the constants c_5, c_6, c_7

the_file = dirname + '/' + kw.fname

print("Writing the atmosphere to %s" %the_file)
print("---------------------------------")
eq.write(the_file)
