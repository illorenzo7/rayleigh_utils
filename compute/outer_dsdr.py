# Author: Loren Matilsky
# Created: well before 05/06/2019
# Computes the necessary outer entropy gradient to carry out a given
# luminosity (default solar) out of a spherical shell with given 
# reference state an dtransport coefficients profile
# Currently reference state can either be a polytrope or from a custom
# file

import numpy as np
import sys, os
from polytrope import compute_polytrope
from common import lsun, rhom, rm, ro

sys.path.append(os.environ['rapp'])
from rayleigh_diagnostics import ReferenceState
from reference_tools import equation_coefficients
from get_parameter import get_parameter

# Set default to compute necessary dsdr to carry out a solar luminosity,
# for kappa_top = 3e12, cp = 3.5x10^8, ri, ro = 5, 6.5860209 x 10^10,
# nrho = 3,
# adiabatic polytrope (n = 1.5), rho_i = 0.18053428 (units all c.g.s.)

lum = lsun
poly_n = 1.5
poly_nrho = 3.0
polytropic_reference = True

# Get directory name
dirname = sys.argv[1]

# Now change the defaults via command-line arguments
args = sys.argv[2:]
nargs = len(args)
for i in range(nargs):
    arg = args[i]
    if arg == '-lum':
        lum = float(args[i+1])
    elif arg == '-n':
        poly_n = float(args[i+1])
    elif arg == '-nrho':
        poly_nrho = float(args[i+1])
    elif arg == '-file':
        polytropic_reference = False
    elif arg == '-rhom':
        rhom = float(args[i+1])
    elif arg == '-rm':
        rm = float(args[i+1])
    elif arg == '-ro':
        ro = float(args[i+1])

nr = 128 # this is arbitrary -- the top thermodynamic values will be
         # independent of the number of interior grid points
         # We just need arrays (nr > 1) for compute_polytrope to work
         # properly

if polytropic_reference:
    di = compute_polytrope(rm, ro, poly_nrho, nr, poly_n, rhom)
    rho = di['density']
    T = di['temperature']
else: # get rho, T from equation_coefficients file, or reference/transport pair
    try:
        eq = equation_coefficients()
        eq.read(dirname + '/equation_coefficients')
        rho = eq.functions[0]
        T = eq.functions[3]
    except:
        ref = ReferenceState(dirname + '/reference')
        rho = ref.density
        T = ref.temperature

# get kappa_top from main_input
kappa_top = get_parameter(dirname, 'kappa_top')

# Now compute desired entropy gradient
flux_top = lum/(4*np.pi*ro**2)
desired_dsdr = -flux_top/rho[0]/T[0]/kappa_top

print('For lum=%1.3e, ro=%1.7e, kappa_top=%1.3e' %(lum, ro, kappa_top))
if polytropic_reference:
    print('and polytropic reference: nrho=%1.1f, poly_n=%1.1f, rhom=%1.7e' %(poly_nrho, poly_n, rhom))
else:
    print('equation_coefficients/reference file ')
print ('Set outer_dsdr (dtdr_top) to %1.8e' %desired_dsdr)
