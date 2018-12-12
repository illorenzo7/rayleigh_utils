###########################
# Author: Loren Matilsky  #
# Revised on: 05/03/2018  #
#############################################################################
# Computes the mass of a spherical shell of a Rayleigh simulation.          #
#############################################################################
# Notes: assumes differentials are already computed via "differentials.py." #
#############################################################################
import numpy as np
import os, sys
from diagnostic_reading import ReferenceState

# Get name of desired run directory
dirname = sys.argv[1]
dirname_stripped = dirname.split('/')[-1]

# Create a data directory if it doesn't exist already
datadir = dirname + '/data/'
if (not os.path.isdir(datadir)):
    os.makedirs(datadir)

# Get basic grid info
rr,tt,cost,sint,rr_depth,ri,ro,d = np.load(datadir + 'grid_info.npy')
nr = len(rr); nt = len(tt)

# Get the differentials dt and dr to integrate over r
dphi, dt, dr = np.load(datadir + 'differentials.npy')

# Get reference state density
ref = ReferenceState(dirname + '/reference', '')
rho = ref.density

# Integrate to get the total shell mass in grams
Mshell_g = np.sum(4.*np.pi*rr**2*rho*dr)

# Get the shell mass in solar masses
Mshell_msun  = Mshell_g/1.989e33

print (dirname_stripped + ' has a shell mass of %e g (or %e M_sun)'\
       %(Mshell_g, Mshell_msun))
np.save(datadir + dirname_stripped + '_Mshell.npy', (Mshell_g, Mshell_msun))