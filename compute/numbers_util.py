import numpy as np
import sys, os
sys.path.append(os.environ['raco'])
from common import *

linebreaks_input = [3, 8, 12]
numbers_input_def = dotdict({
    "aspect": ("A", "r_1/r_2"),
    "nrho": ("N_rho", "ln(rho_1/rho_2)"),
    "dc": ("DC", "exp(N_rho)"),

    "pr": ("Pr", "nu/kappa"),
    "raf": ("Ra_F", "g*F*H^4/(c_p*rho*T*nu*kappa^2)"),
    "di": ("Di", "g*H/(c_p*T)"),
    "bvisc": ("B_visc", "N^2*H^4/nu^2"),
    "he": ("He", "L_star/(F*H^2)"),

    "ek": ("Ek", "nu/(Om_0*H^2)"), 
    #"ta": ("Ta", "1/Ek^2"),
    "rafmod": ("Ra_mod", "Ra_F*Ek^2/Pr"),
    "brot": ("B_rot", "N^2/(Om_0)^2"),

    "prm": ("Pr_m", "nu/eta"),
    "ekm": ("Ek_m", "Ek/Pr_m")
    })


def get_numbers_input(dirname, r1='rmin', r2='rmax', verbose=False):
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
    ir1 = np.argmin(np.abs(rr - r1))
    ir2 = np.argmin(np.abs(rr - r2))
    di.dc = eq.rho[ir1]/eq.rho[ir2]
    di.nrho = np.log(di.dc)

    # get the following volume averages, to be used for all reference_types
    nu_volav = volav_in_radius(dirname, eq.nu, r1, r2)
    kappa_volav = volav_in_radius(dirname, eq.kappa, r1, r2)
    rho_volav = volav_in_radius(dirname, eq.rho, r1, r2)
    tmp_volav = volav_in_radius(dirname, eq.tmp, r1, r2)
    grav_volav = volav_in_radius(dirname, eq.grav, r1, r2)

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

    # B_visc requires some special attention
    if eq.reference_type in [2, 4]:
        nsq_volav = volav_in_radius(dirname, eq.nsq, r1, r2)
    else:
        if advect_reference_state:
            bvisc = get_parameter(dirname, "Buoyancy_Number_Visc")
            pr = get_parameter(dirname, "Prandtl_Number")
            if bvisc is None:
                brot = get_parameter(dirname, "Buoyancy_Number_Rot")
                ek = get_parameter(dirname, "Ekman_Number")
                bvisc = brot/ek**2
            ra = get_parameter(dirname, "Rayleigh_Number")
            if ra is None:
                ramod = get_parameter(dirname, "Modified_Rayleigh_Number")
                ek = get_parameter(dirname, "Ekman_Number")
                ra = ramod*pr/ek**2

            nsq_volav = ra/(pr*bvisc)*volav_in_radius(eq.grav*eq.dsdr, eq.nsq, r1, r2)

    if eq.reference_type in [2, 4]:
        # Prandtl number
        di.pr = nu_volav/kappa_volav

        # flux rayleigh number
        shell_depth = r2 - r1
        di.raf = grav_volav*flux_volav*shell_depth**4/(eq.c_p*rho_volav*tmp_volav*nu_volav*kappa_volav**2)

        # dissipation number
        di.di = grav_volav*shell_depth/(eq.c_p*tmp_volav)

        # non-D heating integral
        if eq.heating_type > 0:
            di.he = lum_volint/(flux_volav*shell_depth**2)
        else:
            di.he = 0.0

        # buoyancy number (viscous)
        if advect_reference_state:
            di.bvisc = nsq_volav*shell_depth**4/nu_volav**2
        else:
            di.bvisc = 0.0

        if rotation:
            # Ekman and Taylor
            di.ek = nu_volav/(eq.om0*shell_depth**2)
            #di.ta = 1.0/di.ek**2

            # modified Rayleigh
            di.rafmod = di.raf*di.ek**2/di.pr

            # buoyancy number (rotational)
            di.brot = di.bvisc*di.ek**2

        if magnetism:
            # magnetic Prandtl
            eta_volav = volav_in_radius(dirname, eq.eta, r1, r2)
            di.prm = nu_volav/eta_volav

            # "magnetic Ekman number"
            di.ekm = eta_volav/(eq.om0*shell_depth**2)
    else: # reference_type = 1,3,5
        # get the original ("global") numbers from main_input,
        # then adjust them in proportion to the variuos volume averages
        # (this assumes we have non-dimensionalized using volume averages)

        # also need non-D of time (won't work for Boussinesq)
        nd_time_visc = get_parameter(dirname, "ND_Time_Visc")
        nd_time_rot = get_parameter(dirname, "ND_Time_Rot")

        # Prandtl number
        pr = get_parameter(dirname, "Prandtl_Number")
        di.pr = pr * (nu_volav/kappa_volav)

        # flux Rayleigh number
        ra = get_parameter(dirname, "Rayleigh_Number")
        if ra is None:
            ramod = get_parameter(dirname, "Modified_Rayleigh_Number")
            ek = get_parameter(dirname, "Ekman_Number")
            ra = ramod*pr/ek**2

        # flux rayleigh number
        di.raf = ra *  ( grav_volav*flux_volav/(rho_volav*tmp_volav*nu_volav*kappa_volav**2) )
        # also remember flux_volav included factor of 1/Pr or Ek/Pr
        if nd_time_visc:
            di.raf *= di.pr
        elif nd_time_rot:
            ek = get_parameter(dirname, "Ekman_Number")
            di.raf *= (di.pr/ek)

        # dissipation number (only for reference type 5)
        # not Boussinesq and I don't use ref type 3
        if eq.reference_type == 5:
            if nd_time_visc:
                di_loc = eq.constants[7] * (ra/pr)
            elif nd_time_rot:
                di_loc = eq.constants[7] * ramod
            di.di = di_loc * (grav_volav/tmp_volav)

        # non-D heating integral
        if eq.heating_type > 0:
            di.he = lum_volint/(flux_volav)
        else:
            di.he = 0.0

        if advect_reference_state:
            # buoyancy number (viscous)
            bvisc = get_parameter(dirname, "Buoyancy_Number_Visc")
            if bvisc is None:
                brot = get_parameter(dirname, "Buoyancy_Number_Rot")
                ek = get_parameter(dirname, "Ekman_Number")
                bvisc = brot/ek**2
            di.bvisc = bvisc * (nsq_volav*nu_volav**2)
        else:
            di.bvisc = 0.0

        if rotation:
            # Ekman and Taylor
            ek = get_parameter(dirname, "Ekman_Number")
            di.ek = ek * nu_volav
            #di.ta = 1.0/di.ek**2

            # modified Rayleigh
            di.rafmod = di.raf*di.ek**2/di.pr

            # buoyancy number (rotational)
            di.brot = di.bvisc*di.ek**2

        if magnetism:
            # magnetic Prandtl
            eta_volav = volav_in_radius(dirname, eq.eta, r1, r2)
            di.prm = nu_volav/eta_volav

            # "magnetic Ekman number"
            ekm = di.ek / di.prm
            di.ekm = ekm * eta_volav

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


def get_numbers_output(dirname, r1='rmin', r2='rmax', the_file=None, the_file_az=None, verbose=False):
    # get diagnostic numbers (e.g., Re and Ro), quantities vol. avg.'d 
    # between r1 and r2
   
    # desired shell
    r1, r2 = interpret_rvals(dirname, np.array([r1, r2]))
    shell_depth = r2 - r1

    # dictionary for output, rotation, magnetism
    di = dotdict()
    rotation = get_parameter(dirname, 'rotation')
    magnetism = get_parameter(dirname, 'magnetism')

    # get reference state
    eq = get_eq(dirname, verbose=verbose)
    rr = eq.rr

    # get field amplitudes
    di_amp_vsr = field_amp(dirname, verbose=verbose) # this one contains full radial profiles
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

    # get the Mach numbers
    dlnprs = eq.dlnrho + eq.dlntmp
    dprs = eq.prs*dlnprs
    drho = eq.rho*eq.dlnrho
    csq = dprs/drho
    csq_volav = volav_in_radius(dirname, csq, r1, r2)

    di.ma = di_amp.v/np.sqrt(csq_volav)
    di.mamean = di_amp.vmean/np.sqrt(csq_volav)
    di.mafluc = di_amp.vfluc/np.sqrt(csq_volav)

    # get the system Reynolds numbers
    nu_volav = volav_in_radius(dirname, eq.nu, r1, r2)
    di.re = di_amp.v*shell_depth/nu_volav
    di.remean = di_amp.vmean*shell_depth/nu_volav
    di.refluc = di_amp.vfluc*shell_depth/nu_volav

    # get the vorticity ("real") Reynolds numbers
    di.revort = di_amp.v**2/di_amp.om/nu_volav
    di.revortmean = di_amp.vmean**2/di_amp.ommean/nu_volav
    di.revortfluc = di_amp.vfluc**2/di_amp.omfluc/nu_volav

    # get ratios of KE in mean vs. fluc flows
    ke = eq.rho*di_amp_vsr.v**2/2
    kemean = eq.rho*di_amp_vsr.vmean**2/2
    kefluc = eq.rho*di_amp_vsr.vfluc**2/2
    
    ke_volav = volav_in_radius(dirname, ke, r1, r2)
    kemean_volav = volav_in_radius(dirname, kemean, r1, r2)
    kefluc_volav = volav_in_radius(dirname, kefluc, r1, r2)

    di.kemean = kemean_volav/ke_volav
    di.kefluc = kefluc_volav/ke_volav

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
        if verbose:
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
        diffrot = (roteq - rotpol)/om0
        di.diffrot = volav_in_radius(dirname, diffrot, r1, r2)

    # magnetic numbers
    if magnetism:
        # system magnetic Reynolds numbers
        eta_volav = volav_in_radius(dirname, eq.eta, r1, r2)
        di.rem = di_amp.v*shell_depth/eta_volav
        di.remmean = di_amp.vmean*shell_depth/eta_volav
        di.remfluc = di_amp.vfluc*shell_depth/eta_volav

        # current ("real") magnetic Reynolds numbers
        di.remcur = di_amp.v*(di_amp.b/di_amp.j)/eta_volav
        di.remcurmean = di_amp.vmean*(di_amp.bmean/di_amp.jmean)/eta_volav
        di.remcurfluc = di_amp.vfluc*(di_amp.bfluc/di_amp.jfluc)/eta_volav

        # plasma beta
        pgas_volav = volav_in_radius(dirname, eq.prs, r1, r2)
        pmag = di_amp_vsr.b**2/(8*np.pi)
        pmag_volav = volav_in_radius(dirname, pmag, r1, r2)
        di.beta = pgas_volav/pmag_volav

        # ratio of mag. energy to kin. energy
        di.me = pmag_volav/ke_volav

    return di
