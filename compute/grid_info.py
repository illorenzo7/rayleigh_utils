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
sys.path.append(os.environ['raco'])
from compute_grid_info import compute_grid_info
dirname = sys.argv[1]

fname = 'grid_info'
args = sys.argv[2:]
nargs = len(args)
for i in range(nargs):
    arg = args[i]
    if arg == '-fname':
        fname = args[i+1]

# Get relevant info from main_input file
nt = get_parameter(dirname, 'n_theta')
use_extrema = get_parameter(dirname, 'use_extrema')

ncheby, domain_bounds = get_domain_bounds(dirname)
nr, nt, nphi, r, rw, tt, cost, sint, tw, phi, dphi =\
        compute_grid_info(domain_bounds, ncheby, nt,\
        use_extrema=use_extrema)

# Write the data
f = open(dirname + '/' + fname, 'wb')
sigpi = np.array(314, dtype=np.int32)
nr = np.array(nr, dtype=np.int32)
nt = np.array(nt, dtype=np.int32)
nphi = np.array(nphi, dtype=np.int32)
dphi = np.array(dphi, dtype=np.float64)

f.write(sigpi.tobytes())
f.write(nr.tobytes())
f.write(nt.tobytes())
f.write(nphi.tobytes())
f.write(r.tobytes())
f.write(rw.tobytes())
f.write(tt.tobytes())
f.write(cost.tobytes())
f.write(sint.tobytes())
f.write(tw.tobytes())
f.write(phi.tobytes())
f.write(dphi.tobytes())
f.close()
