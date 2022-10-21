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

def get_numbers_output(dirname, r1='rmin', r2='rmax'):
    di = dotdict()
    rotation = get_parameter(dirname, 'rotation')
    magnetism = get_parameter(dirname, 'magnetism')


    # get non-rotating, non-magnetic numbers first:

def nonD_numbers(dirname, rbcz=None):
    # all the nonD numbers (as functions of radius and in different zones)
    # we could ever want

    # Make empty dictionary for length_scale arrays
    di_out = dict([])

    # See if run is magnetic
    magnetism = get_parameter(dirname, 'magnetism')
    rotation = get_parameter(dirname, 'rotation')

    # get reference state
    eq = get_eq(dirname)
    rr = eq.radius
    nr = len(rr)
    #di_out['rr'] = rr
    #di_out['nr'] = nr

    di_amp = field_amp(dirname)
    di_len = length_scales(dirname)

    # get the reference state
    eq = get_eq(dirname)

    # get the Reynolds numbers
    shell_depth = di_len['shell_depth']
    hrho = di_len['L_rho']

    di_out['Re'] = di_amp['vamp']*shell_depth/eq.nu
    di_out['Re_fluc'] = di_amp['vfluc']*shell_depth/eq.nu
    di_out['Re_mean'] = di_amp['vmean']*shell_depth/eq.nu

    di_out['Rehrho'] = di_amp['vamp']*hrho/eq.nu
    di_out['Rehrho_fluc'] = di_amp['vfluc']*hrho/eq.nu
    di_out['Rehrho_mean'] = di_amp['vfluc']*hrho/eq.nu

    L_om = di_len['L_om']
    di_out['Revort'] = di_amp['vamp']*L_om/eq.nu
    di_out['Revort_fluc'] = di_amp['vfluc']*L_om/eq.nu
    di_out['Revort_mean'] = di_amp['vmean']*L_om/eq.nu

    # Read in the Shell_Spectra data
    datadir = dirname + '/data/'
    the_file = get_widest_range_file(datadir, 'Shell_Spectra')
    if the_file == '':
        have_spec = False
    else: 
        have_spec = True

    if have_spec:
        ir_spec = di_len['ir_spec']
        rr_spec = di_len['rr_spec']
        L_v = di_len['L_v']
        di_out['Respec'] = (di_amp['vamp']/eq.nu)[ir_spec]*L_v
        di_out['Respec_fluc'] = (di_amp['vfluc']/eq.nu)[ir_spec]*L_v
        di_out['Respec_mean'] = (di_amp['vmean']/eq.nu)[ir_spec]*L_v

    if magnetism: # magnetic Reynolds numbers Rm
        L_J = di_len['L_J']
        di_out['Rm'] = di_amp['vamp']*L_J/eq.eta
        di_out['Rm_fluc'] = di_amp['vfluc']*L_J/eq.eta
        di_out['Rm_mean'] = di_amp['vmean']*L_J/eq.eta

        if have_spec:
            L_B = di_len['L_B']
            di_out['Rmspec'] = (di_amp['vamp']/eq.eta)[ir_spec]*L_B
            di_out['Rmspec_fluc'] = (di_amp['vfluc']/eq.eta)[ir_spec]*L_B
            di_out['Rmspec_mean'] = (di_amp['vmean']/eq.eta)[ir_spec]*L_B

    if rotation: # Rossby numbers
        Om0 = 2*np.pi/compute_Prot(dirname)
        di_out['Ro'] = di_amp['vamp']/(2.0*Om0*shell_depth)
        di_out['Ro_fluc'] = di_amp['vfluc']/(2.0*Om0*shell_depth)
        di_out['Ro_mean'] = di_amp['vmean']/(2.0*Om0*shell_depth)

        di_out['Rohrho'] = di_amp['vamp']/(2.0*Om0*hrho)
        di_out['Rohrho_fluc'] = di_amp['vfluc']/(2.0*Om0*hrho)
        di_out['Rohrho_mean'] = di_amp['vmean']/(2.0*Om0*hrho)

        di_out['Rovort'] = di_amp['vamp']/(2.0*Om0*L_om)
        di_out['Rovort_fluc'] = di_amp['vfluc']/(2.0*Om0*L_om)
        di_out['Rovort_mean'] = di_amp['vmean']/(2.0*Om0*L_om)

        if have_spec:
            di_out['Rospec'] = (di_amp['vamp']/eq.eta)[ir_spec]/(2.0*Om0*L_v)
            di_out['Rospec_fluc'] = (di_amp['vfluc']/eq.eta)[ir_spec]/(2.0*Om0*L_v)
            di_out['Rospec_mean'] = (di_amp['vmean']/eq.eta)[ir_spec]/(2.0*Om0*L_v)

    # now compute the global average of all numbers
    gi = GridInfo(dirname + '/grid_info', '')
    rw = gi.rweights
    if not rbcz is None:
        irbcz = np.argmin(np.abs(rr/rsun - rbcz))
        if have_spec:
            irbcz_spec = np.argmin(np.abs(rr_spec/rsun - rbcz))
        if not (irbcz == 0 or irbcz == nr - 1):
            rwcz = rw[:irbcz+1]/np.sum(rw[:irbcz+1])
            rwrz = rw[irbcz+1:]/np.sum(rw[irbcz+1:])
        else:
            print ('nonD_numbers(): dude, you entered a stupid value for')
            print ('rbcz. you set rbcz = %1.3e' %rbcz)
            print ('it needs be in the range [%.3f, %.3f]' %(np.min(rr)/rsun, np.max(rr)/rsun))
            print ('resetting rbcz = None')
            rbcz = None

    all_keys = list(di_out.keys())
    for key in all_keys:
        if 'spec' in key:
            di_out[key + '_gav'] = np.mean(di_out[key])
        else:
            di_out[key + '_gav'] = np.sum(di_out[key]*rw)
        if not rbcz is None:
            if 'spec' in key:
                if not (irbcz_spec == 0 or irbcz_spec == len(rr_spec) - 1):
                    di_out[key + '_cz'] = np.mean(di_out[key][:irbcz_spec+1])
                    di_out[key + '_rz'] = np.mean(di_out[key][irbcz_spec+1:])
                else:
                    di_out[key + '_cz'] = di_out[key]
                    di_out[key + '_rz'] = di_out[key]
            else:
                di_out[key + '_cz'] = np.sum(di_out[key][:irbcz+1]*rwcz)
                di_out[key + '_rz'] = np.sum(di_out[key][irbcz+1:]*rwrz)
    # I think we got it all!
    return di_out

