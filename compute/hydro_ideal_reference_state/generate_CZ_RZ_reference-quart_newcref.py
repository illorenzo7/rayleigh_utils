# Author: Loren Matilsky
# Created: 05/01/2019
# Purpose: generate a binary file (for Rayleigh to read) that contains
# a reference state, consisting of a stable region (stiffness k) underlying 
# a neutrally stable CZ

# Parameters: output_dir (first argument), 
# Inner, outer, and transition radii (default ri = 3.4139791e10, 
# rm = 5e10, ro = 6.5860209e10...for CZ corresponding to 3 density
# scale heights )
# nr1 (number of Chebyshevs to use for CZ, -nr1) default 96
# nr2 (number of Chebyshevs to use for RZ, -nr2) default 64
# Density at transition boundary (-rhom) default 0.18053428
# Temperature at transition boundary (-tm) default 2.111256e6
# Stiffness k (-k) default 2
# Transition width delta (-delta) (as a fraction of rm) default 0.005
# gamma (-gam) default 5/3
# pressure specific heat cp (-cp) default 3.5e8
# central mass (-M) default M_sun

import numpy as np
import matplotlib.pyplot as plt
import sys, os
from arbitrary_atmosphere import arbitrary_atmosphere
sys.path.append(os.environ['co'])
sys.path.append(os.environ['rapp'])
from write_reference import write_reference
import reference_tools as rt

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
delta = 0.015*ro
nr = 5000 # make the grid super-fine by default

# For now I will also have to initialize transport coefficients here :-(
diff_delta = 0.01
diff_r = 4.87e10

nu_top = 3.0e12
nu_power = -0.5
nu_min = 3.0e9

kappa_top = 3.0e12
kappa_power = -0.5
kappa_min = 3.0e9

# and the heating :-( :-(
heating_delta = 0.005
heating_r = 5.0e10

# Get directory to save binary files for reference state and heating
dirname = sys.argv[1]

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
    elif arg == '-nr':
        nr = int(args[i+1])
    elif arg == '-rhom':
        rhom = float(args[i+1])        
    elif arg == '-tm':
        tm = float(args[i+1])   
    elif arg == '-k':
        k = float(args[i+1])  
    elif arg == '-delta':
        delta = float(args[i+1])*ro
    elif arg == '-gam':
        gam = float(args[i+1]) 
    elif arg == '-cp':
        cp = float(args[i+1]) 
    elif arg == '-M':
        M = float(args[i+1])  
    elif arg == '-nr':
        nr = int(args[i+1])
    elif arg == '-diff_delta':
        diff_delta = float(args[i+1])
    elif arg == '-diff_r':
        diff_r = float(args[i+1])
    elif arg == '-heating_delta':
        heating_delta = float(args[i+1])
    elif arg == '-heating_r':
        heating_r = float(args[i+1])
    elif arg == '-nu_top':
        nu_top = float(args[i+1])
    elif arg == '-nu_power':
        nu_power = float(args[i+1])
    elif arg == '-nu_min':
        nu_min = float(args[i+1])
    elif arg == '-kappa_top':
        kappa_top = float(args[i+1])
    elif arg == '-kappa_power':
        kappa_power = float(args[i+1])
    elif arg == '-kappa_min':
        kappa_min = float(args[i+1])

        
# First, compute reference state on evenly spaced grid, possibly letting
# Rayleigh interpolate later    
rr = np.linspace(ri, ro, nr)

# Define an entropy profile that is +1 for r < rm, 0 for r > rm, and 
# continuously differentiable (a quartic) in between

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

g = bc.G*bc.M/rr**2
dgdr = -2.0*g/rr

T, rho, p, dlnT, dlnrho, dlnp, d2lnrho =\
    arbitrary_atmosphere(rr, s, dsdr, d2sdr2, g,\
                         dgdr, rm, Tm, pm, cp, gam)

# Now that density is set, set transport coefficients
rho_top = rho[-1]

transition_width = diff_delta*ro

f1 = 0.5*(1.0 + np.tanh((rr - diff_r)/transition_width))
f2 = 0.5*(1.0 - np.tanh((rr - diff_r)/transition_width))
df1 = 0.5/(np.cosh((rr - diff_r)/transition_width))**2.0/transition_width
df2 = -df1

nu = nu_top*(rho/rho_top)**nu_power*f1 + nu_min*f2
kappa = kappa_top*(rho/rho_top)**kappa_power*f1 + kappa_min*f2

drho = rho*dlnrho
dnu = nu_top*nu_power/rho_top*(rho/rho_top)**(nu_power - 1.0)*drho*f1 +\
        nu_top*(rho/rho_top)**nu_power*df1 + nu_min*df2
dlnu = dnu/nu

dkappa = kappa_top*kappa_power/rho_top*(rho/rho_top)**(kappa_power - 1.0)*\
        drho*f1 + kappa_top*(rho/rho_top)**kappa_power*df1 + kappa_min*df2
dlnkappa = dkappa/kappa

eq = rt.equation_coefficients(rr) # Class to set equation coefficients
    # and write them into binary form
thefile = dirname + '/new_cref_binary'

buoy = rho*g/cp

# Set the heating function
transition_width = heating_delta*ro
f1_heating = 1.0 + np.tanh((rr - heating_r)/transition_width)
heating = rho*T*f1_heating # proportional to pressure
integral = 4*np.pi*np.trapz(heating*rr**2, x=rr)
lsun = 3.846e33
heating = heating/integral*lsun # this is now f_6; c_10 should = 1

# Make a "zero" array for functions we don't want to set
zero = np.zeros(nr)

eq.set_function(rho, 1)
eq.set_function(buoy, 2)
eq.set_function(nu, 3)
eq.set_function(T, 4)
eq.set_function(kappa, 5)
eq.set_function(heating, 6)
eq.set_function(zero, 7) # eta(r) = 0 
eq.set_function(dlnrho, 8)
eq.set_function(d2lnrho, 9)
eq.set_function(dlnT, 10)
# Derivatives of transport coefficients go here...are they needed?
eq.set_function(dlnu, 11)
eq.set_function(dlnkappa, 12)
eq.set_function(zero, 13)
eq.set_function(dsdr, 14)

eq.set_constant(1.0, 1)
eq.set_constant(1.0, 2)
eq.set_constant(1.0, 3)
eq.set_constant(0.0, 4) # multiplies Lorentz force
eq.set_constant(1.0, 5)
eq.set_constant(1.0, 6)
eq.set_constant(0.0, 7) # multiplies eta in ind. eq.
eq.set_constant(1.0, 8)
eq.set_constant(1.0, 9)
eq.set_constant(1.0, 10)

eq.write(thefile)