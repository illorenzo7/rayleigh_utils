# Created: 04/27/2019
# Compute T_i and P_i for our "standard" 3-scale height Rayleigh simulations
# Use these in subsequent calculations of reference states

import numpy as np

ri = 5.0e10
ro = 6.5860209e10
Nrho = 3.
rho_i = 0.18053428
beta = ri/ro
cp = 3.5e8 # This value is wrong for solar metal abundances, but whatevs...
gamma = 5.0/3.0
n0 = 1.5
cv = cp/gamma
R = (gamma - 1.0)*cp/gamma
M = 1.98891e33
G = 6.67e-8
rsun = 6.957e10

# Compute the temperature at the inner boundary (ri = 5.0e10) of a 3-scale-
# height adiabatic CZ, using eq. (26)
exp = np.exp(Nrho/n0)
T_i = (1.0 - beta)*exp/(exp - 1.0) * G*M/((n0 + 1.0)*R*ri)
# Compute the pressure at the inner boundary
p_i = rho_i*R*T_i