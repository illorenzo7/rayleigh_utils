# This file has been re-purposed on 07/15/2019 to compute the grid
# information for a Rayleigh run using the domain bounds, 
# number of radial points in each domain, number of theta points,
# and whether use_extrema is True or False
# Computes the Chebyshev (radial) weights in accordance to the code
# version 0.9.1 as it is NOW, although probably these weights are
# incorrect at the boundary points (Nick has use_extrema = False)

import numpy as np

def compute_theta_grid(nt):
    tt = np.zeros(nt)
    tw = np.zeros(nt)
    # Now compute Legendre collocation points and weights (theta weights)
    coefs = np.zeros(nt + 1)
    coefs[-1] = 1
    cost = np.polynomial.legendre.legroots(coefs)

    coefs_np1 = np.zeros(nt + 2)
    coefs_np1[-1] = 1
    p_np1 = np.polynomial.legendre.legval(cost,coefs_np1)

    tw = (1. - cost**2.)/(nt + 1.)**2./p_np1**2.
    tt = np.arccos(cost)
    return tt, tw

def compute_r_grid(nr, r1, r2, rmin=None, rmax=None):
    # calculate colocation points and integration weights
    # for (sub)domain (r1, r2) in shell (rmin, rmax)

    # assume only one subdomain unless otherwise specified
    if rmin is None:
        rmin = r1
    if rmax is None:
        rmax = r2

    # initialize arrays
    rr = np.zeros(nr)
    rw = np.zeros(nr)

    # Compute the radial collocation points/weights
    x = np.zeros(nr)
    for ix in range(nr):
        x[ix] = np.cos((ix + 0.5)*np.pi/(nr))
        # x decreases from ~0 to ~ -1

    # Transform x --> r via an affine transformation
    x1, x2 = np.min(x), np.max(x)
    rr = r1 + (x - x1)*(r2 - r1)/(x2 - x1)
    int_scale = 3.*np.pi/((rmax**3. - rmin**3.)*nr)*(r2 - r1)/(x2 - x1)
    rw = int_scale*rr**2.*np.sqrt(1. - x**2.)
    rw[0] *= 0.5 # These multiplications are only justified for
    rw[-1] *= 0.5 # Guass-Lobatto (use_extrema = True)
    return rr, rw

def compute_grid_info(domain_bounds, ncheby, nt):
    ndomains = len(ncheby)
    nr = int(np.sum(ncheby))
    rmin, rmax = domain_bounds[0], domain_bounds[-1]

    rr = np.zeros(nr)
    rw = np.zeros(nr)

    # Compute the radial collocation points/weights
    ir2, ir1 = nr - ncheby[0], nr - 1
    for idomain in range(ndomains):
        r1, r2 = domain_bounds[idomain], domain_bounds[idomain+1]
        nr_loc = ncheby[idomain]

        rr_loc, rw_loc = compute_r_grid(nr_loc, r1, r2, rmin, rmax)
        rr[ir2:ir1 + 1] = rr_loc
        rw[ir2:ir1 + 1] = rw_loc
        
        if idomain < ndomains - 1:
            ir2 -= ncheby[idomain + 1]
            ir1 -= ncheby[idomain]

    # Now compute Legendre collocation points and weights (theta weights)
    tt, tw = compute_theta_grid(nt)
    return rr, rw, tt, tw
