import numpy as np
import sys, os
sys.path.append(os.environ['raco'])
from common import *

# header info for input numbers
linebreaks_input = [4, 7, 13]
numbers_input_def = dotdict({
    "aspect": ("A", "r_1/r_2"),
    "nrho": ("N_rho", "ln(rho_1/rho_2)"),
    "dc": ("DC", "exp(N_rho)"),
    "di": ("Di", "g*H/(c_p*T)"),

    "pr": ("Pr", "nu/kappa"),
    "raf": ("Ra_F", "g*F*H^4/(c_p*rho*T*nu*kappa^2)"),
    "bvisc": ("B_visc", "N^2*H^4/nu^2"),

    "ek": ("Ek", "nu/(2*Om_0*H^2)"), 
    "ta": ("Ta", "1/Ek^2"),
    "rafmod": ("Ra_mod", "Ra_F*Ek^2/Pr"),
    "roc": ("Ro_c", "sqrt(Ra_mod)"),
    "brot": ("B_rot", "N^2/(2*Om_0)^2"),
    "sigma": ("sigma", "sqrt(Pr)*N/(2*Om_0)"),

    "prm": ("Pr_m", "nu/eta"),
    "ekm": ("Ek_m", "Ek/Pr_m")
    })


def get_numbers_input(dirname, r1='rmin', r2='rmax', verbose=False, diman=False, use2=True):
    di = dotdict()
    rotation = get_parameter(dirname, 'rotation')
    magnetism = get_parameter(dirname, 'magnetism')
    advect_reference_state = get_parameter(dirname, 'advect_reference_state')

    # get non-rotating, non-magnetic numbers first:

    # aspect ratio
    r1, r2 = interpret_rvals(dirname, np.array([r1, r2]))
    di.aspect = r1/r2

    # density contrast
    eq = get_eq(dirname, verbose=verbose)
    gi = get_grid_info(dirname, verbose=verbose)
    rr = gi.rr
    ir1, ir2 = inds_from_vals(rr, [r1, r2])
    di.dc = eq.rho[ir1]/eq.rho[ir2]
    di.nrho = np.log(di.dc)

    # how we get the rest depends on reference_type
    if eq.reference_type == 2 or diman:
        # get dimensional volume averages
        nu_volav = volav_in_radius(dirname, eq.nu, r1, r2)
        kappa_volav = volav_in_radius(dirname, eq.kappa, r1, r2)
        rho_volav = volav_in_radius(dirname, eq.rho, r1, r2)
        tmp_volav = volav_in_radius(dirname, eq.tmp, r1, r2)
        grav_volav = volav_in_radius(dirname, eq.grav, r1, r2)
        # this is g/c_p for dimensional anelastic
        vol = get_vol(dirname) # make sure to use the full volume
                    # to calculate the non-radiative heat flux vs radius
        lstar = vol*np.sum(eq.heat*gi.rw)
        flux_rad = vol/(4*np.pi*eq.rr**2)*np.cumsum(eq.heat*gi.rw)
        flux_nonrad = lstar/(4*np.pi*eq.rr**2) - flux_rad
        flux_volav = volav_in_radius(dirname, flux_nonrad, r1, r2)

        # use the local shell's volume for the luminositiy in the heating
        # integral
        vol_loc = (4.0*np.pi/3.0)*(r2**3 - r1**3)
        lum_volint = vol_loc*volav_in_radius(dirname, eq.heat, r1, r2)

        # buoyancy frequency
        nsq = eq.grav*eq.dsdr
        nsq_volav = volav_in_radius(dirname, nsq, r1, r2)

        # get numbers from ratios of volume integrals

        # Prandtl number
        di.pr = nu_volav/kappa_volav

        # flux rayleigh number
        shelldepth = r2 - r1
        di.raf = grav_volav*flux_volav*shelldepth**4/(rho_volav*tmp_volav*nu_volav*kappa_volav**2)

        # dissipation number
        di.di = grav_volav*shelldepth/(tmp_volav)

        # buoyancy number (viscous)
        if advect_reference_state:
            di.bvisc = nsq_volav*shelldepth**4/nu_volav**2
        else:
            di.bvisc = 0.0

        if rotation:
            # Ekman and Taylor
            if use2:
                di.ek = nu_volav/(2*eq.omega0*shelldepth**2)
            else:
                di.ek = nu_volav/(eq.omega0*shelldepth**2)
            di.ta = 1.0/di.ek**2

            # modified Rayleigh
            di.rafmod = di.raf*di.ek**2/di.pr
            di.roc = np.sqrt(di.rafmod)

            # buoyancy number (rotational)
            di.brot = di.bvisc*di.ek**2
            di.sigma = np.sqrt(di.brot*di.pr)

        if magnetism:
            # magnetic Prandtl
            eta_volav = volav_in_radius(dirname, eq.eta, r1, r2)
            di.prm = nu_volav/eta_volav

            if rotation:
                # "magnetic Ekman number"
                di.ekm = eta_volav/(eq.omega0*shelldepth**2)

    elif eq.reference_type == 1: # Boussinesq, get nonD from c's
        di.pr = 1./eq.constants[5]
        di.raf = eq.constants[1]*di.pr
        di.di = 0.
        di.bvisc = 0.0

        if rotation:
            di.ek = 1./eq.constants[2]
            di.ta = 1.0/di.ek**2
            di.rafmod = di.raf*di.ek**2/di.pr
            di.roc = np.sqrt(di.rafmod)
            di.brot = 0.
            di.sigma = 0.

        if magnetism:
            di.prm = 1.0/eq.constants[6]
            if rotation:
                di.ekm = di.ek/di.prm

    elif eq.reference_type == 5: # General anelastic, get nonD from c's
        di.pr = 1./eq.constants[5]
        di.raf = eq.constants[1]*di.pr
        di.di = eq.constants[7]*di.raf/di.pr
        di.bvisc = eq.constants[1]*\
                volav_in_radius(dirname, eq.grav*eq.dsdr, r1, r2)

        if rotation:
            di.ek = 1./eq.constants[0]
            di.ta = 1.0/di.ek**2
            di.rafmod = di.raf*di.ek**2/di.pr
            di.roc = np.sqrt(di.rafmod)
            di.brot = di.bvisc*di.ek**2
            di.sigma = np.sqrt(di.brot*di.pr)

        if magnetism:
            di.prm = 1.0/eq.constants[6]
            if rotation:
                di.ekm = di.ek/di.prm

    return di

# header info for output numbers
linebreaks_output = [3, 6, 9, 12, 15, 16, 19, 22, 23, 24]
numbers_output_def = dotdict({
    "re": ("Re", "v*H/nu"),
    "remean": ("Re_mean", "<v>*H/nu"),
    "refluc": ("Re_fluc", "v'*H/nu"),

    "revort": ("Re_vort", "v^2/(om*nu)"),
    "revortmean": ("Re_vort,mean", "<v>^2/(<om>*nu)"),
    "revortfluc": ("Re_vort,fluc", "v'^2/(om'*nu)"),

    "kedr": ("KE_DR", "rho<v_phi>^2"),
    "kemc": ("KE_MC", "rho<v_m>^2"),
    "kefluc": ("KE_fluc", "rho v'^2"),

    "ro": ("Ro", "v/(2*H*Om_0)"),
    "romean": ("Ro_mean", "<v>/(2*H*Om_0)"),
    "rofluc": ("Ro_fluc", "v'/(2*H*Om_0)"),

    "rovort": ("Ro_vort", "om/(2*Om_0)"),
    "rovortmean": ("Ro_vort,mean", "<om>/(2*Om_0)"),
    "rovortfluc": ("Ro_vort,fluc", "om'/(2*Om_0)"),

    #"diffrot": ("DR", "<|Om_eq - Om_60|^2>^(1/2)"),
    "diffrot": ("DR", "<Om_eq - Om_60>"),

    "rem": ("Re_m", "v*H/eta"),
    "remmean": ("Re_m,mean", "<v>*H/eta"),
    "remfluc": ("Re_m,fluc", "v'*H/eta"),

    "remcur": ("Re_m,cur", "v*(B/J)/eta"),
    "remcurmean": ("Re_m,cur,mean", "<v>*(<B>/<J>)/eta"),
    "remcurfluc": ("Re_m,cur,fluc", "v'*(B'/J')/eta"),

    "me": ("ME", "(B^2/(8*pi)) / (rho*v^2/2)"),

    "dtmp": ("Delta S", "S(r_1) - S(r_2)") })

def get_dr_contrast(dirname, r1='rmin', r2='rmax', lat1=0., lat2=60., the_file=None, verbose=False, alt=False, norms=False, ntheta=None):
    # rotation contrast
    if the_file is None:
        datadir = dirname + '/data/'
        the_file = get_widest_range_file(datadir, 'AZ_Avgs')
    if verbose:
        print ("get_dr_contrast(): reading " + the_file)
    di = get_dict(the_file)
    vals = di['vals']
    lut = di['lut']
    vp_av = vals[:, :, lut[3]]

    # get grid info
    gi = get_grid_info(dirname, ntheta=ntheta)

    # get background stuff
    eq = get_eq(dirname)

    # Get differential rotation in the rotating frame. 
    rotrate = vp_av/gi.xx

    if alt:
        # spherical rms. differential rotation
        dr_contrast = np.sqrt(np.sum(rotrate**2*gi.tw_2d, axis=0))
    else:
        nthalf = gi.nt//2
        tt_lat_sym = gi.tt_lat[nthalf:]
        rotrate_sym = (rotrate[nthalf:, :] + (rotrate[nthalf-1::-1, :]))/2
        ilat1, ilat2 = inds_from_vals(tt_lat_sym, [lat1, lat2])
        dr_contrast = rotrate_sym[ilat1, :] - rotrate_sym[ilat2, :]

    if norms: # don't do RMS of contrast, just average
        out = volav_in_radius(dirname, dr_contrast, r1, r2)/eq.omega0
    else:
        out = volav_in_radius(dirname, dr_contrast**2, r1, r2)**0.5/eq.omega0
    return out

def get_numbers_output(dirname, r1='rmin', r2='rmax', the_file=None, the_file_az=None, verbose=False, shelldepth=None):
    # get diagnostic numbers (e.g., Re and Ro), quantities vol. avg.'d 
    # between r1 and r2
   
    # desired shell
    r1, r2 = interpret_rvals(dirname, np.array([r1, r2]))
    if shelldepth is None:
        shelldepth = r2 - r1


    # dictionary for output, rotation, magnetism
    di = dotdict()
    rotation = get_parameter(dirname, 'rotation')
    magnetism = get_parameter(dirname, 'magnetism')

    # need input numbers for some things
    di_input = get_numbers_input(dirname, r1, r2)

    # get reference state
    eq = get_eq(dirname, verbose=verbose)
    rr = eq.rr

    ir1, ir2 = inds_from_vals(rr, [r1, r2])

    # get shell averaged data for some things
    datadir = dirname + '/data/'
    if the_file is None:
        the_file = get_widest_range_file(datadir, 'Shell_Avgs')
    if verbose:
        print ("get_numbers_output(): reading " + the_file)
    di_shav = get_dict(the_file)
    vals = di_shav['vals']
    lut = di_shav['lut']

    # get field amplitudes
    di_amp_vsr = field_amp(dirname, the_file=the_file, verbose=verbose) # this one contains full radial profiles
    di_amp = dotdict(di_amp_vsr) # profiles averaged between (r1, r2)
    # average each radial profile between (r1, r2):
    for key, profile in di_amp.items():
        if not key in ['rr', 'iter1', 'iter2']:
            # rr, iter1, iter2 shouldn't be averaged
            
            # see if it's a thermal variable or field variable
            # (the field variables should be averaged in quadrature)
            fieldvar = False
            if key[0] in ['v', 'b', 'j']:
                fieldvar = True
            if key[:2] == 'om':
                fieldvar = True

            if fieldvar:
                volav_sq = volav_in_radius(dirname, profile**2, r1, r2)
                di_amp[key] = np.sqrt(volav_sq)
            else:
                di_amp[key] = volav_in_radius(dirname, profile, r1, r2)

    # get non-rotating, non-magnetic numbers first:

    # get the system Reynolds numbers
    nu_volav = volav_in_radius(dirname, eq.constants[4]*eq.functions[2,:], r1, r2)
    di.re = di_amp.v*shelldepth/nu_volav
    di.remean = di_amp.vmean*shelldepth/nu_volav
    di.refluc = di_amp.vfluc*shelldepth/nu_volav

    # get the vorticity ("real") Reynolds numbers
    di.revort = di_amp.v**2/di_amp.om/nu_volav
    di.revortmean = di_amp.vmean**2/di_amp.ommean/nu_volav
    di.revortfluc = di_amp.vfluc**2/di_amp.omfluc/nu_volav

    # get estimated (and real) potential energy
    tmp_vsr = vals[:, 0, lut[501]]
    dtmp = tmp_vsr[ir1] - tmp_vsr[ir2] # achieved temperature difference across shell

    # achieved potential energy across shell
    grav_volav = volav_in_radius(dirname, eq.grav, r1, r2)
    # this is g/c_p for dimensional anelastic
    rho_volav = volav_in_radius(dirname, eq.rho, r1, r2)
    kappa_volav = volav_in_radius(dirname, eq.kappa, r1, r2)

    # get ratios of KE in mean vs. fluc flows
    ke = eq.rho*di_amp_vsr.v**2/2
    kedr = eq.rho*di_amp_vsr.vpmean**2/2
    kemc = eq.rho*di_amp_vsr.vpolmean**2/2
    kefluc = eq.rho*di_amp_vsr.vfluc**2/2

    ke_volav = volav_in_radius(dirname, ke, r1, r2)
    kedr_volav = volav_in_radius(dirname, kedr, r1, r2)
    kemc_volav = volav_in_radius(dirname, kemc, r1, r2)
    kefluc_volav = volav_in_radius(dirname, kefluc, r1, r2)

    di.kedr = kedr_volav#%/ke_volav
    di.kemc = kemc_volav#%/ke_volav
    di.kefluc = kefluc_volav#/ke_volav


    # rotational numbers
    if rotation:
        omega0 = eq.omega0
        
        # get the system Rossby numbers
        di.ro = di_amp.v/(2.0*omega0*shelldepth)
        di.romean = di_amp.vmean/(2.0*omega0*shelldepth)
        di.rofluc = di_amp.vfluc/(2.0*omega0*shelldepth)

        # get the vorticity ("real") Rossby numbers
        di.rovort = di_amp.om/(2.0*omega0)
        di.rovortmean = di_amp.ommean/(2.0*omega0)
        di.rovortfluc = di_amp.omfluc/(2.0*omega0)

        di.diffrot = get_dr_contrast(dirname, r1=r1, r2=r2, the_file=the_file_az, verbose=verbose, norms=True)

    # magnetic numbers
    if magnetism:
        # system magnetic Reynolds numbers
        eta_volav = volav_in_radius(dirname, eq.eta, r1, r2)
        di.rem = di_amp.v*shelldepth/eta_volav
        di.remmean = di_amp.vmean*shelldepth/eta_volav
        di.remfluc = di_amp.vfluc*shelldepth/eta_volav

        # current ("real") magnetic Reynolds numbers
        di.remcur = di_amp.v*(di_amp.b/di_amp.j)/eta_volav
        di.remcurmean = di_amp.vmean*(di_amp.bmean/di_amp.jmean)/eta_volav
        di.remcurfluc = di_amp.vfluc*(di_amp.bfluc/di_amp.jfluc)/eta_volav

        # plasma beta
        pgas_volav = volav_in_radius(dirname, eq.prs, r1, r2)
        pmag = di_amp_vsr.b**2/(8*np.pi)
        pmag_volav = volav_in_radius(dirname, pmag, r1, r2)
        #di.beta = pgas_volav/pmag_volav

        # ratio of mag. energy to kin. energy
        di.me = pmag_volav#/ke_volav

    di.dtmp = dtmp

    return di
