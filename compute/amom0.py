# This file has been re-purposed on 07/15/2019 to compute the grid
# information for a Rayleigh run using the domain bounds, 
# number of radial points in each domain, number of theta points,
# and whether use_extrema is True or False
# Computes the Chebyshev (radial) weights in accordance to the code
# version 0.9.1 as it is NOW, although probably these weights are
# incorrect when use_extrema = False

import numpy as np
import sys, os
sys.path.append(os.environ['rapp'])
sys.path.append(os.environ['co'])
from get_parameter import get_parameter
from rayleigh_diagnostics import GridInfo, ReferenceState
from common import strip_dirname
dirname = sys.argv[1]
dirname_stripped = strip_dirname(dirname)

# Get relevant info from main_input file and grid_info/reference
om0 = get_parameter(dirname, 'angular_velocity')
ref = ReferenceState(dirname + '/reference')
gi = GridInfo(dirname + '/grid_info')
rr = gi.radius
ri, ro = np.min(rr), np.max(rr)
rw = gi.rweights
dr = rw*(1./3.)*(ro**3 - ri**3)/rr**2.

rho = ref.density
integral = np.sum(rho*rr**4.*dr)
amom0 = 8.*np.pi*om0/3.*integral

volume = 4./3.*np.pi*(ro**3. - ri**3.)

print ("The shell angular momentum density associated with " +\
        dirname_stripped + " is:")
print ("%1.3e g cm^(-1) s^(-1)" %(amom0/volume))
