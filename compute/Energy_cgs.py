# Author: Loren Matilsky
# Created: 09/13/2019
# This script computes the ratio of volume-integrated magnetic energy to 
# volume integrated kinetic energy (also time-averaged) using G_Avgs data
# Displays the computed ratio at the terminal

import numpy as np
import sys, os
sys.path.append(os.environ['rapp'])
sys.path.append(os.environ['raco'])
from common import *

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
magnetism = get_parameter(dirname, 'magnetism')

# Get the kinetic energies
KE_tot = vals[lut[401]]
KE_DR = vals[lut[408]]
KE_MC = vals[lut[406]] + vals[lut[407]]
KE_C = vals[lut[409]]

# Print the energies to 3 s.f.
print('-------------------------')
print('percent of total (KE and ME separate) at right')
print('-------------------------')
print('KE_tot\t\t%1.2e' %KE_tot)
print('KE_DR\t\t%1.2e\t%.03f' %(KE_DR, 100*KE_DR/KE_tot))
print('KE_MC\t\t%1.2e\t%.03f' %(KE_MC, 100*KE_MC/KE_tot))
print('KE_C\t\t%1.2e\t%.03f' %(KE_C, 100*KE_C/KE_tot))

# Get the magnetic energy 
if magnetism:
    ME_tot = vals[lut[1101]]
    ME_phi = vals[lut[1108]]
    ME_m = vals[lut[1106]] + vals[lut[1107]]
    ME_C = vals[lut[1109]]

    # Print the energies to 3 s.f.
    print('-------------------------')
    print('ME_tot\t\t%1.2e' %ME_tot)
    print('ME_phi\t\t%1.2e\t%.03f' %(ME_phi, 100*ME_phi/ME_tot))
    print('ME_m\t\t%1.2e\t%.03f' %(ME_m, 100*ME_m/ME_tot))
    print('ME_C\t\t%1.2e\t%.03f' %(ME_C, 100*ME_C/ME_tot))
