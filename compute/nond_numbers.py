import numpy as np
from scipy.integrate import cumtrapz
import sys, os
sys.path.append(os.environ['rapp'])
sys.path.append(os.environ['raco'])
from rayleigh_diagnostics import Shell_Avgs, GridInfo
from common import *

def get_numbers_input(dirname, r1='rmin', r2='rmax'):
    di = dotdict()
    rotation = get_parameter(dirname, 'rotation')
    magnetism = get_parameter(dirname, 'magnetism')

    # get non-rotating, non-magnetic numbers first:

    # aspect ratio
    r1, r2 = interpret_rvals(dirname, np.array([r1, r2]))
    di.r1, di.r2 = r1, r2
    di.aspect = r2/r1

    # density contrast
    eq = get_eq(dirname)
    gi = get_grid_info(dirname)
    rr = gi.rr
    ir1 = np.argmin(np.abs(rr - r1))
    ir2 = np.argmin(np.abs(rr - r2))
    di.dc = eq.rho[ir1]/eq.rho[ir2]

    # Prandtl number
    nu_volav = volav_in_radius(dirname, eq.nu, r1, r2)
    k_volav = volav_in_radius(dirname, eq.kappa, r1, r2)
    di.pr = nu_volav/k_volav

    # flux rayleigh number
    vol = get_vol(dirname, r1, r2)
    flux_rad = vol/(4*np.pi*eq.rr**2)*np.cumsum(eq.heat*gi.rw)
    flux_nonrad = eq.lum/(4*np.pi*eq.rr**2) - flux_rad
    flux_volav = volav_in_radius(dirname, flux_nonrad, r1, r2)
    rho_volav = volav_in_radius(dirname, eq.rho, r1, r2)
    tmp_volav = volav_in_radius(dirname, eq.tmp, r1, r2)
    kappa_volav = volav_in_radius(dirname, eq.kappa, r1, r2)
    grav_volav = volav_in_radius(dirname, np.abs(eq.grav), r1, r2)
    shell_depth = r2 - r1
    di.raf = grav_volav*flux_volav*shell_depth**4/(eq.c_p*rho_volav*tmp_volav*nu_volav*kappa_volav**2)

    return di
