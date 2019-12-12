# Author: Loren Matilsky
# Created: 05/18/2019
# This script computes the flux Rayleigh number for a Rayleigh run in
# directory [dirname], following the definition (eq. 14) in Featherstone
# & Hindman (2016) (ApJ paper, not letter)
# Takes the heat flux from Shell_Avgs file and the rest from the reference
# state file and grid_info for the radial weights
# Displays the computed Rayleigh number at the terminal. 

import numpy as np
import sys, os
sys.path.append(os.environ['rapp'])
sys.path.append(os.environ['raco'])
from get_parameter import get_parameter
from rayleigh_diagnostics import Shell_Avgs, ReferenceState, GridInfo,\
        TransportCoeffs
from reference_tools import equation_coefficients

# Get directory name
dirname = sys.argv[1]

# Read in reference state, grid info for radial weights,
# and transport coefficients at the top of the domain
gi = GridInfo(dirname + '/grid_info')
rw = gi.rweights



# Read in one Shell_Avgs file to get heat flux,
# the last one just because
files = os.listdir(dirname + '/Shell_Avgs')
files.sort()
sh = Shell_Avgs(dirname + '/Shell_Avgs/' + files[-1], '')
# Get the shell depth:
H = np.max(sh.radius) - np.min(sh.radius)


# Read in necessary radial profiles for the flux Rayleigh number
# "vsr" for "profile vs radius"
F_rad = sh.vals[:, 0, sh.lut[1433], 0]
# Get the luminosity from F_rad
lum = F_rad[-1]*4*np.pi*np.min(sh.radius)**2
print ("Luminosity calculated from F_rad is %1.3e" %lum)
F_vsr = lum/4/np.pi/sh.radius**2 - F_rad

# Get reference info from reference/transport or equation_coefficients
cp = 3.5e8
try:
    ref = ReferenceState(dirname + '/reference')
    t = TransportCoeffs(dirname + '/transport')
    nu_vsr = t.nu
    kappa_vsr = t.kappa
    g_vsr = ref.gravity
    rho_vsr = ref.density
    T_vsr = ref.temperature
    print ("Got reference info from 'reference'/'transport'")
except:
    eq = equation_coefficients()
    eq.read(dirname + '/equation_coefficients')
    rho_vsr = eq.functions[0]
    T_vsr = eq.functions[3]
    g_vsr = eq.functions[1]*cp/rho_vsr
    nu_vsr = eq.functions[2]
    kappa_vsr = eq.functions[4]
    print ("Got reference info from 'equation_coefficients'")

# Compute volume-averages of the radial profiles using the radial 
# integration weights
F = np.sum(rw*F_vsr)
g = np.sum(rw*g_vsr)
rho = np.sum(rw*rho_vsr)
T = np.sum(rw*T_vsr)
nu = np.sum(rw*nu_vsr)
kappa = np.sum(rw*kappa_vsr)

# Compute the flux Rayleigh number!
Ra = g*F*H**4/(cp*rho*T*nu*kappa**2)

# And print it
print("The flux Rayleigh number is %1.3e" %Ra)

# Also write it as an empty file in [dirname]
# Remove all other "RaF_is_" that might be present
# (in case this definition of Ra is replaced
names = os.listdir(dirname)
for name in names:
    if "RaF_is_" in name:
        os.remove(dirname + '/' + name)
fname = dirname + ("/00_RaF_is_%1.3e" %Ra)
f = open(fname, "w")
f.close()
