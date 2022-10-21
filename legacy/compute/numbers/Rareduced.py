# Author: Loren Matilsky
# Created: 01/17/2020
# This script computes the volume-averaged Ekman number for a 
# Rayleigh run in directory [dirname], using the (constant) shell depth
# as the length scale. 
# Gets diffusion profiles from transport or equation_coefficients
# Reads grid_info for the radial weights
# Displays the computed Ekman number at the terminal

import numpy as np
import sys, os
sys.path.append(os.environ['rapp'])
sys.path.append(os.environ['raco'])
from rayleigh_diagnostics import Shell_Avgs, GridInfo
from common import *

# Get directory name
dirname = sys.argv[1]

# Read in grid info for radial weights and reference velocity
gi = GridInfo(dirname + '/grid_info')
rw = gi.rweights
rr = gi.radius
H = np.max(rr) - np.min(rr)

# Read in transport coefficients for nu-profile
eq = get_eq(dirname)
nu = eq.nu

# Get angular velocity 
Om0 = 2.*np.pi/compute_Prot(dirname)

# Compute volume-averaged Ekman number
# using the radial integration weights
Ek_vs_r = nu/(2*H**2*Om0) # note there is technically a 1/sin(theta) in the 
# Definition of Ek, but this averages to 1 over latitude
Ek = np.sum(rw*Ek_vs_r)

try:
    # By default, read in one Shell_Avgs file to get heat flux,
    # the last one just because
    files = os.listdir(dirname + '/Shell_Avgs')
    files.sort()
    sh = Shell_Avgs(dirname + '/Shell_Avgs/' + files[-1], '')
    F_rad = sh.vals[:, 0, sh.lut[1433], 0]
    rr = sh.radius
    print ("Got F_rad and rr from Shell_Avgs/" + files[-1])
except:
    # otherwise, get data from time-averaged data product
    datadir = dirname + '/data/'
    the_file = get_widest_range_file(datadir, 'Shell_Avgs')
    di = get_dict(datadir + the_file)
    F_rad = di['vals'][:, di['lut'][1433]]
    rr = di['rr']
    print ("Got F_rad and rr from ", the_file)

# Get the shell depth:
H = np.max(rr) - np.min(rr)

# Read in necessary radial profiles for the flux Rayleigh number
# "vsr" for "profile vs radius"

# Get the luminosity from F_rad
lum = F_rad[-1]*4*np.pi*np.min(rr)**2
print ("Luminosity calculated from F_rad is %1.3e" %lum)
F_vsr = lum/4/np.pi/rr**2 - F_rad

# Get reference info from reference/transport or equation_coefficients
cp = 3.5e8
nu_vsr = eq.nu
kappa_vsr = eq.kappa
g_vsr = eq.gravity
rho_vsr = eq.density
T_vsr = eq.temperature

# Compute volume-averages of the radial profiles using the radial 
# integration weights
F = np.sum(rw*F_vsr)
g = np.sum(rw*g_vsr)
rho = np.sum(rw*rho_vsr)
T = np.sum(rw*T_vsr)
nu = np.sum(rw*nu_vsr)
kappa = np.sum(rw*kappa_vsr)

# Compute the flux Rayleigh number!
Ra = g*F*H**4./(cp*rho*T*nu*kappa**2.)

# Compute the reduced Rayleigh number
R = Ra*Ek**(4./3.)

# And print it
print ("Ra_F = %1.3e" %Ra)
print ("Ek = %1.3e" %Ek)
print("The reduced Rayleigh number Ra_F*Ek**(4/3) is %1.3e" %R)
