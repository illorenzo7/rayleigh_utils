# Routine to compute the shell volume for a simulation (in c.c.)
# Created: 11/08/2019 
# Can specify rmin/rmax to get partial volume; otherwise, returns the
# volume of the whole shell

import numpy as np
import os, sys
sys.path.append(os.environ['raco'])
sys.path.append(os.environ['rapp'])
from get_parameter import get_parameter
from rayleigh_diagnostics import GridInfo

dirname = sys.argv[1]

gi = GridInfo(dirname + '/grid_info', '')

# By default integrate over the whole shell
rr = gi.radius
nr = len(rr)
ir_min = 0
ir_max = nr - 1

# Read in CLAs
args = sys.argv[2:]
nargs = len(args)

for i in range(nargs):
    arg = args[i]
    if arg == '-rmin':
        desired_rmin = float(args[i+1])
        ir_max = np.argmin(np.abs(rr - desired_rmin)) # Note that since the radius
        # arrays are reversed, rmin sets ir_max ...
    elif arg == '-rmax':
        desired_rmax = float(args[i+1])
        ir_min = np.argmin(np.abs(rr - desired_rmax))

# Get the actual rmin/rmax values associated with the index
rmin = rr[ir_max] # remember reversed radius array ...
rmax = rr[ir_min]

# Compute shell volume
shell_volume = 4./3.*np.pi*(rmax**3. - rmin**3.)

print ("Volume between rmin = %1.3e cm and rmax = %1.3e cm:" %(rmin, rmax))
print ("%1.3e cm^3" %shell_volume)
