import numpy as np
import sys, os
sys.path.append(os.environ['raco'])
from common import *

linebreaks_input = [3, 7, 11]
numbers_input_def = dotdict({
    "aspect": ("A", "r_1/r_2"),
    "nrho": ("N_rho", "ln(rho_1/rho_2)"),
    "dc": ("DC", "exp(N_rho)"),

    "pr": ("Pr", "nu/kappa"),
    "raf": ("Ra_F", "g*F*H^4/(c_p*rho*T*nu*kappa^2)"),
    "di": ("Di", "g*H/(c_p*T)"),
    "soundbuoy": ("CsN", "(c/H)^2/N^2)"),

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

    # sound crossing frequency to buoyancy

    # buoyancy
    nsq_volav = volav_in_radius(dirname, eq.nsq, r1, r2)

    # sound speed squared
    dlnprs = eq.dlnrho + eq.dlntmp
    dprs = eq.prs*dlnprs
    drho = eq.rho*eq.dlnrho
    csq = dprs/drho
    csq_volav = volav_in_radius(dirname, csq, r1, r2)
    di.soundbuoy = (csq_volav/shell_depth*2)/nsq_volav

    if rotation:
        # Ekman and Taylor
        di.ek = nu_volav/(eq.om0*shell_depth**2)
        di.ta = 1/di.ek**2

        # buoyancy
        di.buoy = nsq_volav/eq.om0**2

        # ratio of rotation period to sound crossing time (squared)
        di.sound = (csq_volav/shell_depth*2)/eq.om0**2

    if magnetism:
        # magnetic Prandtl
        eta_volav = volav_in_radius(dirname, eq.eta, r1, r2)
        di.prm = nu_volav/eta_volav

        # "magnetic Ekman number"
        di.ekm = eta_volav/(eq.om0*shell_depth**2)

    return di

# numbers groups
numbers_output_groups = ["Mach numbers", "Reynolds numbers", "vort. Reynolds num.", "KE fractions",\
        "Rossby numbers", "vort. Rossby num.", "DR fraction",\
        "mag. Reynolds num.", "mag. current Reyn.", "plasma beta", "ME fraction"]
numbers_output_ngroup = 4
numbers_output_ngroup_rot = 3
numbers_output_ngroup_mag = 4

linebreaks_output = [3, 6, 9, 11, 14, 17, 18, 21, 24, 25]
numbers_output_def = dotdict({
    "ma": ("Ma", "v/c"),
    "mamean": ("Ma_mean","<v>/c"),
    "mafluc": ("Ma_fluc","v'/c"),

    "re": ("Re", "v*H/nu"),
    "remean": ("Re_mean", "<v>*H/nu"),
    "refluc": ("Re_fluc", "v'*H/nu"),

    "revort": ("Re_vort", "v^2/(om*nu)"),
    "revortmean": ("Re_vort,mean", "<v>^2/(<om>*nu)"),
    "revortfluc": ("Re_vort,fluc", "v'^2/(om'*nu)"),

    "kemean": ("KE_mean", "<v>^2/v^2"),
    "kefluc": ("KE_fluc", "v'^2/v^2"),

    "ro": ("Ro", "v/(2*H*Om_0)"),
    "romean": ("Ro_mean", "<v>/(2*H*Om_0)"),
    "rofluc": ("Ro_fluc", "v'/(2*H*Om_0)"),

    "rovort": ("Ro_vort", "om/(2*Om_0)"),
    "rovortmean": ("Ro_vort,mean", "<om>/(2*Om_0)"),
    "rovortfluc": ("Ro_vort,fluc", "om'/(2*Om_0)"),

    "diffrot": ("DR", "(Om_eq - Om_60)/Om_0"),

    "rem": ("Re_m", "v*H/eta"),
    "remmean": ("Re_m,mean", "<v>*H/eta"),
    "remfluc": ("Re_m,fluc", "v'*H/eta"),

    "remcur": ("Re_m,cur", "v*(B/J)/eta"),
    "remcurmean": ("Re_m,cur,mean", "<v>*(<B>/<J>)/eta"),
    "remcurfluc": ("Re_m,cur,fluc", "v'*(B'/J')/eta"),

    "beta": ("plasma beta", "8*pi*P/B^2"),
    
    "me": ("ME", "(B^2/(8*pi)) / (rho*v^2/2)") })


def get_numbers_output(dirname, shell_depth=None, the_file=None, the_file_az=None):
    # get diagnostic numbers (e.g., Re and Ro) as functions of radius

    # dictionary for output
    di = dotdict()
    rotation = get_parameter(dirname, 'rotation')
    magnetism = get_parameter(dirname, 'magnetism')

    # get reference state
    eq = get_eq(dirname)
    rr = eq.rr
    if shell_depth is None:
        shell_depth = np.max(rr) - np.min(rr)

    # get field amplitudes and length scales
    di_amp = field_amp(dirname)

    # get non-rotating, non-magnetic numbers first:

    # get the Mach numbers
    dlnprs = eq.dlnrho + eq.dlntmp
    dprs = eq.prs*dlnprs
    drho = eq.rho*eq.dlnrho
    csq = dprs/drho

    di.ma = di_amp.v/np.sqrt(csq)
    di.mamean = di_amp.vmean/np.sqrt(csq)
    di.mafluc = di_amp.vfluc/np.sqrt(csq)

    # get the system Reynolds numbers
    di.re = di_amp['v']*shell_depth/eq.nu
    di.remean = di_amp['vmean']*shell_depth/eq.nu
    di.refluc = di_amp['vfluc']*shell_depth/eq.nu

    # get the vorticity ("real") Reynolds numbers
    di.revort = di_amp['v']**2/di_amp.om/eq.nu
    di.revortmean = di_amp['vmean']**2/di_amp.ommean/eq.nu
    di.revortfluc = di_amp['vfluc']**2/di_amp.omfluc/eq.nu

    # get ratios of KE in mean vs. fluc flows
    ke = eq.rho*di_amp.v**2/2
    kemean = eq.rho*di_amp.vmean**2/2
    kefluc = eq.rho*di_amp.vfluc**2/2

    di.kemean = kemean/ke
    di.kefluc = kefluc/ke

    # rotational numbers
    if rotation:
        om0 = eq.om0
        
        # get the system Rossby numbers
        di.ro = di_amp.v/(2.0*om0*shell_depth)
        di.romean = di_amp.vmean/(2.0*om0*shell_depth)
        di.rofluc = di_amp.vfluc/(2.0*om0*shell_depth)

        # get the vorticity ("real") Rossby numbers
        di.rovort = di_amp.om/(2.0*om0)
        di.rovortmean = di_amp.ommean/(2.0*om0)
        di.rovortfluc = di_amp.omfluc/(2.0*om0)

        # rotation contrast
        datadir = dirname + '/data/'
        if the_file_az is None:
            the_file_az = get_widest_range_file(datadir, 'AZ_Avgs')
        print ("get_numbers_output(): reading " + the_file_az)
        di_az = get_dict(the_file_az)
        vals_az = di_az['vals']
        lut_az = di_az['lut']
        vp_av = vals_az[:, :, lut_az[3]]

        # Get necessary grid info
        gi = get_grid_info(dirname)

        # Get differential rotation in the rotating frame. 
        rotrate = vp_av/gi.xx

        # rotation contrast between equator and 60 degrees
        latcut = 60
        iteq = np.argmin(np.abs(gi.tt_lat))
        itnorth = np.argmin(np.abs(gi.tt_lat - latcut))
        itsouth = np.argmin(np.abs(gi.tt_lat + latcut))
        roteq = rotrate[iteq, :]
        rotpol = 0.5*(rotrate[itnorth, :] + rotrate[itsouth, :])
        di.diffrot = (roteq - rotpol)/om0

    # magnetic numbers
    if magnetism:
        # system magnetic Reynolds numbers
        di.rem = di_amp['v']*shell_depth/eq.eta
        di.remmean = di_amp['vmean']*shell_depth/eq.eta
        di.remfluc = di_amp['vfluc']*shell_depth/eq.eta

        # current ("real") magnetic Reynolds numbers
        di.remcur = di_amp['v']*(di_amp.b/di_amp.j)/eq.eta
        di.remcurmean = di_amp['vmean']*(di_amp.bmean/di_amp.jmean)/eq.eta
        di.remcurfluc = di_amp['vfluc']*(di_amp.bfluc/di_amp.jfluc)/eq.eta

        # plasma beta
        pgas = eq.prs
        pmag = di_amp.b**2/(8*np.pi)
        di.beta = pgas/pmag

        # ratio of mag. energy to kin. energy
        di.me = pmag/ke

    return di