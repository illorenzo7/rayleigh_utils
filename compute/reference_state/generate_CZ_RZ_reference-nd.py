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
# --amp : amplitude of 1/cp (dS/dr) 
    # -- yes it's another independent parameter, Loren realized on 04/24

# --fname
# File to save reference state in (default customfile)
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
# use double prec for gamma
# 3 sig figs otherwise
kw_default = dotdict(dict({'alpha': 1., 'beta': 0.759, 'gamma': 1.6666666666666667, 'delta': 0.219, 'nrho': 3., 'fname': 'customfile', 'nr': 10000, 'jup': False, 'amp': 0.456}))
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

# note that whatever we get in the end should be rounded to match
# what is manually entered in the main_input file
# for convenience, round radii to two decimal places
rt = round(rt, 2)
rmin = round(rmin, 2)
rmax = round(rmax, 2)
rbcz = round(rbcz, 2)
rtcz = round(rtcz, 2)
rbrz = round(rbrz, 2)
rtrz = round(rtrz, 2)

# recompute alpha and beta in light of rounding
kw.beta = rbcz/rtcz
kw.alpha = (rtrz - rbrz)/(rtcz - rbcz)

# compute reference state on super-fine grid to interpolate onto later    
r = np.linspace(rmax, rmin, kw.nr) # keep radius in decreasing order for consistency with Rayleigh convention

# Define an entropy profile that is +1 for r < rt - delta, 0 for r > rt,
# and continuously differentiable (a quartic) in between

dsdr = np.zeros(kw.nr)

def psifunc(r, r0, delta):
    arr = np.zeros_like(r)
    for i in range(len(r)):
        rloc = r[i]
        if rloc <= r0:
            arr[i] = 0.
        elif rloc < r0 + delta and rloc > r0:
            x = (rloc - r0)/delta
            arr[i] = 1.0 - (1.0 - x**2.0)**2.0
        else:
            arr[i] = 1.0
    return arr

if kw.jup: # RZ is above
    dsdr = kw.amp*psifunc(r, rtcz, kw.delta)
else: # CZ is above
    dsdr = kw.amp*(1. - psifunc(r, rbcz-kw.delta, kw.delta))

# compute the atmosphere (this depends on dsdr, not nsq)
rho, tmp, dlnrho, d2lnrho, dlnt, g =\
        arbitrary_atmosphere_nd(r, dsdr, rbcz, rtcz, kw.gamma, kw.nrho)

# compute the normalized buoyancy frequency
nsq = g*dsdr
nsq_norm = definite_integral(nsq*r**2, r, rbrz, rtrz)
nsq_norm /= 1./3.*(rtrz**3. - rbrz**3.)
nsq /= nsq_norm

print(buff_line)
print("Computed atmosphere for RZ-CZ, ds/dr joined with a quartic")
if kw.jup:
    print ("geometry : Jovian (RZ atop CZ)")
else:
    print ("geometry : solar (CZ atop RZ)")
print("nr         : %i" %kw.nr) 
print("alpha      : %1.5f" %kw.alpha)
print("beta       : %1.5f" %kw.beta)
if kw.jup:
    print("   (rbcz, rtcz=rbrz, rtrz): (%1.2f, %1.2f, %1.2f)"\
            %(rbcz,rtcz,rtrz))
else:
    print("   (rbrz, rtrz=rbcz, rtcz): (%1.2f, %1.2f, %1.2f)"\
            %(rbrz,rtrz,rtcz))
print("delta      : %1.5f" %kw.delta)
print("gamma      : %1.5f" %kw.gamma)
print("   n=1/(gamma-1)      : %1.5f" %(1./(kw.gamma-1.)))
print("Nrho       : %1.5f" %kw.nrho)
print("amp        : %1.5f" %kw.amp)
print(buff_line)

# Now write to file using the equation_coefficients framework
eq = equation_coefficients(r)

# Set only the thermodynamic functions/constants in this routine
# In other routines, we can set the heating and transport coefficients
# Only set c_4 = 1/(4*pi) if mag = True

print("Setting f_1, f_2, f_4, f_8, f_9, f_10, and f_14")
eq.set_function(rho, 1)
eq.set_function(rho*g, 2)
eq.set_function(tmp, 4)
eq.set_function(dlnrho, 8)
eq.set_function(d2lnrho, 9)
eq.set_function(dlnt, 10)
eq.set_function(nsq/g, 14)

print("Setting all constants to 1.0")
# set all the constants to 1. (no harm I can see for now)
for i in range(eq.nconst):
    eq.set_constant(1.0, i)
the_file = dirname + '/' + kw.fname

print("Writing the atmosphere to %s" %the_file)
eq.write(the_file)

# write metadata to separate file
metafile = the_file + '_meta.txt'
f = open(dirname + '/' + metafile, 'w')

f.write(buff_line + '\n')
f.write("Created custom reference state using the\n")
f.write("generate_CZ_RZ_reference-nd routine.\n")
f.write("The CZ/RZ system has the folowing attributes:\n")

if kw.jup:
     f.write("geometry : Jovian (RZ atop CZ)\n")
else:
     f.write("geometry : solar (CZ atop RZ)\n")
f.write("nr         : %i\n" %kw.nr) 
f.write("alpha      : %1.5f\n" %kw.alpha)
f.write("beta       : %1.5f\n" %kw.beta)
if kw.jup:
    f.write("   (rbcz, rtcz=rbrz, rtrz): (%1.2f, %1.2f, %1.2f)\n"\
            %(rbcz,rtcz,rtrz))
else:
    f.write("   (rbrz, rtrz=rbcz, rtcz): (%1.2f, %1.2f, %1.2f)\n"\
            %(rbrz,rtrz,rtcz))
# somehow these precisions are confusing me...anyway forget it for now
f.write("delta      : %1.5f\n" %kw.delta)
f.write("gamma      : %1.16f\n" %kw.gamma)
f.write("   n=1/(gamma-1)      : %1.5f\n" %(1./(kw.gamma-1.)))
f.write("Nrho       : %1.5f\n" %kw.nrho)
f.write("amp        : %1.5f\n" %kw.amp)
f.write(buff_line + '\n')
f.close()
print("Writing the metadata   to %s" %metafile)
print(buff_line)
