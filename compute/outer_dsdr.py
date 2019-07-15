# Author: Loren Matilsky
# Created: well before 05/06/2019
# Computes the necessary outer entropy gradient to carry out a given
# luminosity (default solar) out of a spherical shell with given 
# reference state an dtransport coefficients profile
# Currently reference state can either be a polytrope or from a custom
# file

import numpy as np
import sys
from polytrope import compute_polytrope
from common import lsun, rhom, rm, ro
from read_reference import read_reference

# Set default to compute necessary dsdr to carry out a solar luminosity,
# for ktop = 3e12, cp = 3.5x10^8, ri, ro = 5, 6.5860209 x 10^10, nrho = 3,
# adiabatic polytrope (n = 1.5), rho_i = 0.18053428 (units all c.g.s.)

lum = lsun
ktop = 3.0e12
poly_n = 1.5
nrho = 3.0
polytropic_reference = True
file_reference = False

# Now change the defaults via command-line arguments
args = sys.argv[1:]
nargs = len(args)
for i in range(nargs):
    arg = args[i]
    if arg == '-lum':
        lum = float(args[i+1])
    elif arg == '-ktop':
        ktop = float(args[i+1])
    elif arg == '-polyn':
        poly_n = float(args[i+1])
    elif arg == '-nrho':
        nrho = float(args[i+1])
    elif arg == '-ref':
        polytropic_reference = False
        file_reference = True
        ref_file = args[i+1]
    elif arg == '-rhom':
        rhom = float(args[i+1])
    elif arg == '-ri':
        ri = float(args[i+1])
    elif arg == '-ro':
        ro = float(args[i+1])

nr = 128 # this is arbitrary -- the top thermodynamic values will be
         # independent of the number of interior grid points

if polytropic_reference:
    di = compute_polytrope(rm, ro, nrho, nr, poly_n, rhom)
    rho = di['density']
    T = di['temperature']
elif file_reference:
    nr, rr, rho, dlnrho, d2lnrho, p, T, dlnT, dsdr, s, g =\
            read_reference(ref_file)
    ro = np.max(rr)

flux_top = lum/(4*np.pi*ro**2)
desired_dsdr = -flux_top/rho[0]/T[0]/ktop

print('For lum=%1.3e, ro=%1.7e, ktop=%1.3e' %(lum, ro, ktop))
if polytropic_reference:
    print('and polytropic reference: nrho=%1.1f, poly_n=%1.1f, rho_i=%1.7e' %(nrho, poly_n, rho_i))
elif file_reference:
    print('and custom reference ' + ref_file)
print ('Set outer_dsdr (dtdr_top) to %1.8e' %desired_dsdr)
