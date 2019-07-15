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
#dirname = sys.argv[1]

def compute_grid_info(domain_bounds, ncheby, nt, use_extrema=False):
    if not isinstance(ncheby, tuple): # this is probably because
        # ncheby was specified as nr, not (nr,)
        ncheby = (ncheby,)
    ndomains = len(ncheby)
    nr = np.sum(ncheby)

    r = np.zeros(nr)
    rw = np.zeros(nr)
    tt = np.zeros(nt)
    tw = np.zeros(nt)

    # Compute the radial collocation points/weights
    ir_min, ir_max = 0, ncheby[0] - 1
    for idomain in range(ndomains):
        rmin, rmax = domain_bounds[idomain], domain_bounds[idomain+1]
        nr_loc = ncheby[idomain]
        x = np.zeros(nr_loc)
        r_loc = np.zeros(nr_loc)
        for ix in range(nr_loc):
            if use_extrema:
                x[ix] = np.cos(ix*np.pi/(nr_loc - 1))
            else:
                x[ix] = np.cos((ix + 0.5)*np.pi/(nr_loc))
        # Transform x --> r via an affine transformation
        xmin, xmax = np.min(x), np.max(x)
        r_loc = rmin + (x - xmin)*(rmax - rmin)/(xmax - xmin)
        int_scale = 3.*np.pi/(rmax**3. - rmin**3.)*nr_loc*\
                (rmax - rmin)/(xmax - xmin)
        rw_loc = int_scale * r_loc**2. * np.sqrt(1. - x**2.)
        rw_loc[0] *= 0.5 # These multiplications are only justified for
        rw_loc[-1] *= 0.5 # Guass-Lobatto (use_extrema = True)

        r[ir_min:ir_max + 1] = r_loc
        rw[ir_min:ir_max + 1] = rw_loc

        if idomain < ndomains - 1:
            ir_min += ncheby[idomain]
            ir_max += ncheby[idomain + 1]

    return (r, rw)
