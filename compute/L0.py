# NEEDS UPDATED DOCS AND TESTING
# Routine to figure out how much simulation time a given directory contains
# Created: 12/23/2018 (but really before)
#
# By default the simulation time is calculated from the longest possible time
# interval in the G_Avgs data (i.e., first output file to last output file)
# User can specify another time interval in the standard way, e.g.,
# -n 10 (last 10 G_Avgs files)
# -range iter1 iter2 (files closest to iter1 to iter2 range)
# -centerrange iter0 nfiles (nfiles centered about iter0)
# etc., etc.
import numpy as np
import os, sys
sys.path.append(os.environ['raco'])
sys.path.append(os.environ['rapp'])
from get_parameter import get_parameter
from rayleigh_diagnostics import GridInfo, ReferenceState

dirname = sys.argv[1]

# get reference state and grid info
ref = ReferenceState(dirname + '/reference', '')
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

# Get the frame rotation rate
Om0 = get_parameter(dirname, 'angular_velocity')

# Get the radial integration weights and density
rw = gi.rweights
rho = ref.density

# Compute L0
scaling = 1./3.*(rmax**3. - rmin**3.)
total_amom = 8.*np.pi*Om0/3.*np.sum((rho*rr**2.*rw)[ir_min:ir_max+1])*scaling
shell_volume = 4./3.*np.pi*(rmax**3. - rmin**3.)
L0 = total_amom/shell_volume

print ("Average amom_z density between ")
print ("rmin = %1.3e cm and rmax = %1.3e cm:" %(rmin, rmax))
print ("%1.3e g/cm/s" %L0)
