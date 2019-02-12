# Computes the necessary outer entropy gradient to carry out a solar 
# luminosity in a fixed-flux outer boundary condition
# This script is designed with the Sun in mind, thus assuming:
# M_sun = 1.98891e33 (in c.g.s.)
# L_sun = 3.846e33
# cP = 3.5e8
# n = 1.5 (polytropic index)
# rho_i = 0.18053428 @ ri = 5e10

import numpy as np
import sys
from polytrope import compute_polytrope

# User supplies Nrho, ro, and kappa_o (in that order)
Nrho = float(sys.argv[1])
ro = float(sys.argv[2])
ktop = float(sys.argv[3])

# Pre-set values appropriate for the sun
ri = 5e10
nr = 128 # this is arbitrary -- the top thermodynamic values will be
         # independent of the number of interior grid points
poly_n = 1.5
rho_i = 0.18053428
L_sun = 3.846e33

di = compute_polytrope(ri, ro, Nrho, nr, poly_n, rho_i)
rho = di['density']
T = di['temperature']

flux_top = L_sun/(4*np.pi*ro**2)
dsdr = -flux_top/rho[0]/T[0]/ktop

print ('Set drdr_top to %1.8e' %dsdr)
