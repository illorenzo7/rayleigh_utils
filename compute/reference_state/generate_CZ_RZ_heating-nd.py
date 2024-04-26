#/ Author: Loren Matilsky
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
#
# --rt
# transition radius for heating
#
# --rstar (default rsun)
#
# --delta
# Transition width delta default 0.050*rsun
#
# --lum
# Luminosity to be driven through layer, default Lsun
# 
# --rmax
# radius at which to normalize integral (rmin doesn't matter because
# the tanh will be ~zero here; revisit for very thin RZs)
#
# --fname
# default custom_reference_binary. Place to read current reference state and# write heating to

import numpy as np
import sys, os
from scipy.integrate import simps

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
kw_default = dotdict(dict({'rmax': sun.r_nrho3, 'rt': sun.r_bcz, 'rstar': rsun, 'delta': 0.05*sun.r, 'lum': sun.l, 'fname': 'custom_reference_binary'}))

# overwrite defaults
kw = update_dict(kw_default, clas)

# check for bad keys
find_bad_keys(kw_default, clas, clas0['routinename'], justwarn=True)

# Open and read the hopefully already existing reference file!
eq = equation_coefficients()
the_file = dirname + '/' + kw.fname
eq.read(the_file)
rr = eq.radius
nr = eq.nr
rho = eq.functions[0]
tmp = eq.functions[3]
smooth = 0.5*(1.0 + np.tanh((rr - kw.rt)/kw.delta)) # "detects" CZ
profile = (rho*tmp - rho[0]*tmp[0])*smooth

irmax = np.argmin(np.abs(rr - kw.rmax))
int_profile = -4*np.pi*simps((profile*rr**2.0)[irmax:], rr[irmax:])
# remember rr is reversed

radial_shape = profile/int_profile

print(buff_line)
print("Computed heating for RZ-CZ, joined with tanh")
print("nr                  : %i" %nr) 
print("rstar               : %1.8e cm" %kw.rstar) 
print("rt/rstar            : %1.8e" %(kw.rt/kw.rstar) )
print("delta/rstar         : %.8f"  %(kw.delta/kw.rstar))
print("rmax/rstar          : %.8f"  %(kw.rmax/kw.rstar))
print("c_10 (lum)          : %1.8e erg/s"  %kw.lum)
print(buff_line)

# Now write to file using the equation_coefficients framework
print("Setting f_6 and c_10")

eq.set_function(radial_shape, 6)
eq.set_constant(kw.lum, 10) # Total luminosity normalizes heating function

print("Writing the heating to %s" %the_file)
print("---------------------------------")
eq.write(the_file)
