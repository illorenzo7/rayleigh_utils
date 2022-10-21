import numpy as np
from scipy.integrate import cumtrapz
import sys, os
sys.path.append(os.environ['rapp'])
sys.path.append(os.environ['raco'])
from rayleigh_diagnostics import Shell_Avgs, GridInfo
from common import *

numbers_input_def = dotdict({
    "aspect": ("A", "r_1/r_2"),
    "nrho": ("N_rho", "ln(rho_1/rho_2)"),
    "dc": ("DC", "exp(N_rho)"),
    "pr": ("Pr", "nu/kappa"),
    "raf": ("Ra_F", "g*F*H^4/(c_p*rho*T*nu*kappa^2)"),
    "di": ("Di", "g*H/(c_p*T)"),
    "ek": ("Ek", "nu/(Om_0*H^2)"), 
    "ta": ("Ta", "1/Ek^2"),
    "buoy": ("B", "N^2/Om_0^2"),
    "sound": ("SC", "(c/H)^2/Om_0^2"),
    "prm": ("Pr_m", "nu/eta"),
    "ekm": ("Ek_m", "Ek/Pr_m")
    })


def get_numbers_input(dirname, r1='rmin', r2='rmax'):
    di = dotdict()
    rotation = get_parameter(dirname, 'rotation')
    magnetism = get_parameter(dirname, 'magnetism')

    # get non-rotating, non-magnetic numbers first:

    # aspect ratio
    r1, r2 = interpret_rvals(dirname, np.array([r1, r2]))
    di.aspect = r1/r2

    # density contrast
    eq = get_eq(dirname)
    gi = get_grid_info(dirname)
    rr = gi.rr
    ir1 = np.argmin(np.abs(rr - r1))
    ir2 = np.argmin(np.abs(rr - r2))
    di.dc = eq.rho[ir1]/eq.rho[ir2]
    di.nrho = np.log(di.dc)

    # Prandtl number
    nu_volav = volav_in_radius(dirname, eq.nu, r1, r2)
    kappa_volav = volav_in_radius(dirname, eq.kappa, r1, r2)
    di.pr = nu_volav/kappa_volav

    # flux rayleigh number
    vol = get_vol(dirname) # make sure to get the whole volume here...
    flux_rad = vol/(4*np.pi*eq.rr**2)*np.cumsum(eq.heat*gi.rw)
    flux_nonrad = eq.lum/(4*np.pi*eq.rr**2) - flux_rad
    flux_volav = volav_in_radius(dirname, flux_nonrad, r1, r2)
    rho_volav = volav_in_radius(dirname, eq.rho, r1, r2)
    tmp_volav = volav_in_radius(dirname, eq.tmp, r1, r2)
    grav_volav = volav_in_radius(dirname, np.abs(eq.grav), r1, r2)
    shell_depth = r2 - r1
    di.raf = grav_volav*flux_volav*shell_depth**4/(eq.c_p*rho_volav*tmp_volav*nu_volav*kappa_volav**2)

    # dissipation number
    di.di = grav_volav*shell_depth/(eq.c_p*tmp_volav)

    if rotation:
        # Ekman and Taylor
        di.ek = nu_volav/(eq.om0*shell_depth**2)
        di.ta = 1/di.ek**2

        # buoyancy
        nsq_volav = volav_in_radius(dirname, eq.nsq, r1, r2)
        di.buoy = nsq_volav/eq.om0**2

        # ratio of rotation period to sound crossing time (squared)
        dlnprs = eq.dlnrho + eq.dlntmp
        dprs = eq.prs*dlnprs
        drho = eq.rho*eq.dlnrho
        csq = dprs/drho
        csq_volav = volav_in_radius(dirname, csq, r1, r2)
        di.sound = (csq_volav/shell_depth*2)/eq.om0**2

    if magnetism:
        # magnetic Prandtl
        eta_volav = volav_in_radius(dirname, eq.eta, r1, r2)
        di.prm = nu_volav/eta_volav

        # "magnetic Ekman number"
        di.ekm = eta_volav/(eq.om0*shell_depth**2)

    return di
