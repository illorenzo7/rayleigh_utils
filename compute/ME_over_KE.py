# Author: Loren Matilsky
# Created: 09/13/2019
# This script computes the ratio of volume-integrated magnetic energy to 
# volume integrated kinetic energy (also time-averaged) using G_Avgs data
# Displays the computed ratio at the terminal

import numpy as np
import sys, os
sys.path.append(os.environ['rapp'])
sys.path.append(os.environ['raco'])
from get_parameter import get_parameter
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

# Get the magnetic energy 
ME_dens = vals[lut[1101]]
KE_dens = vals[lut[401]]

ratio = ME_dens/KE_dens

# Print the ratio
print("Ratio of ME to KE: %1.3e" %ratio)
