# Author: Loren Matilsky
# Created: 05/01/2019
# Modified to work with new cref framework: 10/20/2019
# New CLA interface: 08/30/2021
#
# Purpose: generate a binary file (for Rayleigh to read) that contains
# a reference state, consisting of a stable region (stiffness k) underlying 
# a neutrally stable CZ
# assumes g(r) = -G * mstar / r^2

# Parameters: output_dir (first argument), 
# Command-line options:
#
# --rmin, --rmax, --rt
# Inner, outer, and transition radii (default rmin = 0.1*rstar
# (i.e. very deep shell; don't need to use all of it)
# rt = 5e10, rmax = 0.99*rsun (as close to surface as my simple
# hydro model allows)

# --mstar
# central mass, default mstar = msun
# 
# --rstar 
# stellar radius, default rstar = rsun
#
# --rhot
# Density at transition boundary, default rhobcz ~ 0.18
#
# --tempt
# Temperature at transition boundary default tempbcz ~ 2.1e6
#
# --k
# Stiffness k, default 2.0 (define in units of dsdr in RZ / (cp/rstar) for consistency
# 
# --delta
# Transition width delta as a fraction of rstar, default 0.05
#
# --gam
# specific heat ratio gamma, default thermo_gamma = 5.0/3.0
#
# --cp
# pressure specific heat cp, default c_P = 3.5e8
#
# --mag
# Whether magnetism is True or False, default False (hydro)

import numpy as np
import sys, os
from arbitrary_atmosphere import arbitrary_atmosphere

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
kw_default = dotdict(dict({'rmin': 0.1, 'rmax': 0.99, 'rt': rbcz_nond, 'mstar': msun, 'rstar': rsun, 'rhot': rhobcz, 'tempt': tempbcz, 'k': 2.0, 'delta': 0.05, 'gam': thermo_gamma, 'cp': c_P, 'mag': False, 'fname': 'custom_reference_binary', 'nr': 10000}))

# overwrite defaults
kw = update_dict(kw_default, clas)

# check for bad keys
find_bad_keys(kw_default, clas, clas0['routinename'], justwarn=True)

# compute gas constant R
kw.R = kw.cp*(kw.gam-1)/kw.gam
prst = kw.rhot*kw.R*kw.tempt

# Make everything physical length
kw.delta *= kw.rstar
kw.rmin *= kw.rstar
kw.rt *= kw.rstar
kw.rmax *= kw.rstar

# compute reference state on super-fine grid to interpolate onto later    
rr = np.linspace(kw.rmax, kw.rmin, kw.nr) # keep radius in decreasing order for consistency with Rayleigh convention

# Define an entropy profile that is +1 for r < rt - delta, 0 for r > rt,
# and continuously differentiable (a quartic) in between

s = np.zeros(kw.nr)
dsdr = np.zeros(kw.nr)
d2sdr2 = np.zeros(kw.nr)

for i in range(kw.nr):
    rloc = rr[i]
    if rloc <= kw.rt - kw.delta:
        s[i] = (8.0/15.0)*(kw.k*kw.cp)*(kw.delta/kw.rstar) +\
                kw.k*kw.cp*(rloc/kw.rstar - kw.rt/kw.rstar)
        dsdr[i] = kw.k*kw.cp/kw.rstar
        d2sdr2[i] = 0.0
    elif rloc > kw.rt - kw.delta and rloc < kw.rt:
        x = (rloc - kw.rt)/kw.delta
        s[i] = (kw.k*kw.cp)*(kw.delta/kw.rstar)*((2.0/3.0)*x**3.0 - (1.0/5.0)*x**5.0)
        dsdr[i] = (kw.k*kw.cp/kw.rstar)*(1.0 - (1.0 - x**2.0)**2.0)
        d2sdr2[i] = (4.0/kw.delta)*(kw.k*kw.cp/kw.rstar)*(1.0 - x**2.0)*x
    else:
        s[i] = 0.0
        dsdr[i] = 0.0
        d2sdr2[i] = 0.0

# Make gravity due to a point mass at the origin
g = G*kw.mstar/rr**2.0
dgdr = -2.0*g/rr

T, rho, p, dlnT, dlnrho, dlnp, d2lnrho =\
    arbitrary_atmosphere(rr, s, dsdr, d2sdr2, g,\
                         dgdr, kw.rt, kw.tempt, prst, kw.cp, kw.gam)

print(buff_line)
print("Computed atmosphere for RZ-CZ, ds/dr joined with a quartic")
print("nr         : %i" %kw.nr) 
print("rstar      : %1.8e cm" %kw.rstar) 
print("mstar      : %1.8e  g" %kw.mstar) 
print("rmin/rstar : %.8f    " %(kw.rmin/kw.rstar))
print("rt/rstar   : %.8f    " %(kw.rt/kw.rstar))
print("rmax/rstar : %.8f    " %(kw.rmax/kw.rstar))
print("delta/rstar: %.8f"  %(kw.delta/kw.rstar))
print("stiffness k: %1.8e"  %kw.k)
print("stiffness  : %1.8e"  %(kw.k*kw.cp/kw.rstar))
print("stiffness = (dsdr in RZ)")
print("k = (stiffness)/(cp/rstar)")
print(buff_line)

# Now write to file using the equation_coefficients framework
eq = equation_coefficients(rr)

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
