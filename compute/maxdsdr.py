# Routine to figure out maxabs dsdr of reference state
# Created: 09/02/2020
#
import numpy as np
import os, sys
sys.path.append(os.environ['raco'])
sys.path.append(os.environ['rapp'])
from rayleigh_diagnostics import GridInfo, ReferenceState
from reference_tools import equation_coefficients
from get_eq import get_eq

fname = None
kappa_top = None

# Now change the defaults via command-line arguments
args = sys.argv[2:]
nargs = len(args)
for i in range(nargs):
    arg = args[i]
    if arg == '-fname':
        fname = args[i+1]
    elif arg == '-crb':
        fname = 'custom_reference_binary'

dirname = sys.argv[1]
eq = get_eq(dirname, fname=fname)
dsdr = eq.dsdr
print ("maxabs (dsdr) = %1.3e" %(np.max(np.abs(dsdr))))
