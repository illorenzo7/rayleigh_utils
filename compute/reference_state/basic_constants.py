# Created: 04/27/2019
# Compute T_i and P_i for our "standard" 3-scale height Rayleigh simulations
# Use these in subsequent calculations of reference states

import numpy as np

ri = 4.176e10 # inner boundary of RZ; makes depth of RZ about 0.5 depth of
        # CZ
rm = 5.0e10 # tachocine boundary
ro = 6.5860209e10
Nrho = 3.
rhom = 0.18053428
beta = rm/ro
cp = 3.5e8 # This value is wrong for solar metal abundances, but whatevs...
gamma = 5.0/3.0
n0 = 1.5
cv = cp/gamma
R = (gamma - 1.)*cp/gamma
M = 1.98891e33
G = 6.67e-8
rsun = 6.957e10

# Compute the temperature at the inner boundary (ri = 5.0e10) of a 3-scale-
# height adiabatic CZ, using eq. (26)
exp = np.exp(Nrho/n0)
Tm = (1. - beta)*exp/(exp - 1.)*G*M/((n0 + 1.)*R*rm)
# Compute the pressure at the inner boundary
pm = rhom*R*Tm
