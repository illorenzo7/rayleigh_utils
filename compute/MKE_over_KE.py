# Author: Loren Matilsky
# Created: 09/13/2019
# This script computes the ratio of volume-integrated magnetic energy to 
# volume integrated kinetic energy (also time-averaged) using G_Avgs data
# Displays the computed ratio at the terminal

import numpy as np
import sys, os
sys.path.append(os.environ['rapp'])
sys.path.append(os.environ['raco'])
from common import get_widest_range_file, get_dict

# Get directory name
dirname = sys.argv[1]

# Read in the Shell_Avgs data
datadir = dirname + '/data/'
G_Avgs_file = get_widest_range_file(datadir, 'G_Avgs')
print ('Getting volume-integrated energy info from ' + datadir +\
        G_Avgs_file + ' ...')
di = get_dict(datadir + G_Avgs_file)
vals = di['vals']
lut = di['lut']

# Get the volume-integrated energies
MKE_dens = vals[lut[405]] # mean kin. energy
KE_dens = vals[lut[401]] # total kin. energy

ratio = MKE_dens/KE_dens

# Print the ratio
# This had better be 1 - CKE/KE, so this is a good test to see if all is 
# well with the Universe
print("Ratio of MKE to KE: %1.3e" %ratio)
