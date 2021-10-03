# This file has been re-purposed on 07/15/2019 to compute the grid
# information for a Rayleigh run using the domain bounds, 
# number of radial points in each domain, number of theta points,
# and whether use_extrema is True or False
# Computes the Chebyshev (radial) weights in accordance to the code
# version 0.9.1 as it is NOW, although probably these weights are
# incorrect when use_extrema = False

import numpy as np
from math import factorial

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

def compute_r_grid(nr, rmin, rmax, use_extrema=False):
    r = np.zeros(nr)
    rw = np.zeros(nr)

    # Compute the radial collocation points/weights
    x = np.zeros(nr)
    for ix in range(nr):
        if use_extrema:
            x[ix] = np.cos(ix*np.pi/(nr - 1))
        else:
            x[ix] = np.cos((ix + 0.5)*np.pi/(nr))

    # Transform x --> r via an affine transformation
    xmin, xmax = np.min(x), np.max(x)
    r = rmin + (x - xmin)*(rmax - rmin)/(xmax - xmin)
    int_scale = 3.*np.pi/((rmax**3. - rmin**3.)*nr)*\
            (rmax - rmin)/(xmax - xmin)
    rw = int_scale * r**2. * np.sqrt(1. - x**2.)
    rw[0] *= 0.5 # These multiplications are only justified for
    rw[-1] *= 0.5 # Guass-Lobatto (use_extrema = True)
    return r, rw

def compute_grid_info(domain_bounds, ncheby, nt, use_extrema=False):
    ndomains = len(ncheby)
    nr = np.sum(ncheby)
    ri, ro = domain_bounds[0], domain_bounds[-1]

    r = np.zeros(nr)
    rw = np.zeros(nr)
    tt = np.zeros(nt)
    tw = np.zeros(nt)

    # Compute the radial collocation points/weights
    ir_min, ir_max = nr - ncheby[0], nr - 1
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
        int_scale = 3.*np.pi/((ro**3. - ri**3.)*nr_loc)*\
                (rmax - rmin)/(xmax - xmin)
        rw_loc = int_scale * r_loc**2. * np.sqrt(1. - x**2.)
        rw_loc[0] *= 0.5 # These multiplications are only justified for
        rw_loc[-1] *= 0.5 # Guass-Lobatto (use_extrema = True)

        r[ir_min:ir_max + 1] = r_loc
        rw[ir_min:ir_max + 1] = rw_loc
        
        if idomain < ndomains - 1:
            ir_min -= ncheby[idomain + 1]
            ir_max -= ncheby[idomain]

    # Now compute Legendre collocation points and weights (theta weights)
    coefs = np.zeros(nt + 1)
    coefs[-1] = 1
    cost = np.polynomial.legendre.legroots(coefs)

    coefs_np1 = np.zeros(nt + 2)
    coefs_np1[-1] = 1
    p_np1 = np.polynomial.legendre.legval(cost,coefs_np1)

    tw = (1. - cost**2.)/(nt + 1.)**2./p_np1**2.
    tt = np.arccos(cost)

    # Compute sint and stuff for phi-grid
    sint = np.sin(tt)
    nphi = 2*nt
    phi = 2*np.pi*np.arange(0, nphi)/nphi
    dphi = 2*np.pi/nphi + np.zeros(nphi)
    return nr, nt, nphi, r, rw, tt, cost, sint, tw, phi, dphi
