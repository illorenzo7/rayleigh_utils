# Author: Loren Matilsky
# Created: 12/19/2022
#
# Description: Script to plot 
# each term in the total energy equation in the meridional plane,
# each term's latitudinally (spherically) integrated profile vs. radius,
# and the total integrated rate of change of the energy.

# Plots the total energy equation, KE eq., ME eq., and heat eq. separately
# to choose only one plot, specify
# -tot, -ke, -inte, -me

import numpy as np
import pickle
import matplotlib as mpl
mpl.use('TkAgg')
import matplotlib.pyplot as plt
plt.rcParams['mathtext.fontset'] = 'dejavuserif'
csfont = {'fontname':'DejaVu Serif'}
import sys, os
sys.path.append(os.environ['rapp'])
sys.path.append(os.environ['raco'])
sys.path.append(os.environ['rapl'])
from azav_util import plot_azav
from common import *
from rayleigh_diagnostics import GridInfo
from read_inner_vp import read_inner_vp
from read_eq_vp import read_eq_vp

# Get directory name and stripped_dirname for plotting purposes
dirname = sys.argv[1]
dirname_stripped = strip_dirname(dirname)

# domain bounds
ncheby, domain_bounds = get_domain_bounds(dirname)
ri = np.min(domain_bounds)
ro = np.max(domain_bounds)
d = ro - ri

# Directory with data and plots, make the plotting directory if it doesn't
# already exist    
datadir = dirname + '/data/'

# Read command-line arguments (CLAs)
magnetism = None
saveplot = True
plotcontours = True
plotlatlines = True
minmax = None
linthresh = None
linscale = None
minmaxrz = None
linthreshrz = None
linscalerz = None
AZ_Avgs_file = get_widest_range_file(datadir, 'AZ_Avgs')
forced = False
rvals = []
rbcz = None
symlog = False
crudeint = False # by default use exact latitudinal weights and instead 
# if crudeint = True, use a "crude weight"
tag = ''
plot_az = True
plot_av = True
av_ext = '.png'

# determine magnetism from main_input
magnetism = get_parameter(dirname, 'magnetism')
keys = ['toten_tote', 'toten_ke', 'toten_inte']
if magnetism:
    keys.append('toten_me')

# need stellar luminosity
lstar = get_lum(dirname)

plotdir = None

args = sys.argv[2:]
nargs = len(args)
for i in range(nargs):
    arg = args[i]
    if arg == '-plotdir':
        plotdir = args[i+1]
    if arg == '-minmax':
        minmax = float(args[i+1]), float(args[i+2])
    elif arg == '-minmaxrz':
        minmaxrz = float(args[i+1]), float(args[i+2])
    elif arg == '-rbcz':
        rbcz = float(args[i+1])
        rvals.append(rbcz)
    elif arg == '-nocontour':
        plotcontours = False
    elif arg == '-usefile':
        AZ_Avgs_file = args[i+1]
        AZ_Avgs_file = AZ_Avgs_file.split('/')[-1]
    elif arg == '-forced':
        forced = True
        # check if force_econs is True
        # If False, this work will need to get included in total energy eq
        # If True, energy is conserved, so term in heat equation
        # balances out the work term
        try:
            force_econs = get_parameter(dirname, 'force_econs')
        except:
            force_econs = False
    elif arg == '-rvals':
        rvals_str = args[i+1].split()
        for rval_str in rvals_str:
            rvals.append(float(rval_str))
    elif arg == '-symlog':
        symlog = True
    elif arg == '-linthresh':
        linthresh = float(args[i+1])
    elif arg == '-linscale':
        linscale = float(args[i+1])
    elif arg == '-linthreshrz':
        linthreshrz = float(args[i+1])
    elif arg == '-linscalerz':
        linscalerz = float(args[i+1])
    elif arg == '-nolats':
        plotlatlines = False
    elif arg == '-crudeint':
        crudeint = True
    elif arg == '-tag':
        tag = '_' + args[i+1]
    elif arg == '-ke': # just KE equation
        keys = ['toten_ke']
    elif arg == '-inte': # just heat equation
        keys = ['toten_inte']
    elif arg == '-me': # just ME equation
        keys = ['toten_me']
    elif arg == '-tote':
        keys = ['toten_tote']
    elif arg == '-av':
        plot_az = False
    elif arg == '-az':
        plot_av = False
    elif arg == '-pdf':
        av_ext = '.pdf'

# Get the data:
print ('Getting energy production terms from ' +\
        datadir + AZ_Avgs_file)
di = get_dict(datadir + AZ_Avgs_file)

# Will need integration weights for a few things
gi = GridInfo(dirname + '/grid_info')
rr = gi.radius
rw = gi.rweights
nr = gi.nr
nt = gi.ntheta
if crudeint:
    tt = gi.theta
    tw = np.zeros_like(tt)
    tw[1:-1] = 0.5*(tt[:-2] - tt[2:])
    tw[0] = tw[1]
    tw[-1] = tw[-2]
    tw /= np.sum(tw)
    print ("crudeint = True, using sloppily defined lat. weights")
else:
    tw = gi.tweights
    print ("crudeint = False, using exact lat. weights")
    print ("To use crude weights, specify -crudeint")

# will need shell volume for integrated terms
ri, ro = np.min(rr), np.max(rr)
volume = 4./3.*np.pi*(ro**3. - ri**3.)

# get integration weights
rw_2d = rw.reshape((1, nr))
tw_2d = tw.reshape((nt, 1))
if not rbcz is None:
    irbcz = np.argmin(np.abs(rr - rbcz))
    rbcz_real = rr[irbcz]
    volume_cz = 4./3.*np.pi*(ro**3. - rbcz_real**3.)
    volume_rz = 4./3.*np.pi*(rbcz_real**3. - ri**3.)
    rw_cz = np.copy(rw[:irbcz])
    rw_rz = np.copy(rw[irbcz:])
    rw_cz /= np.sum(rw_cz)
    rw_rz /= np.sum(rw_rz)
    rw_cz_2d = rw_cz.reshape((1, len(rw_cz)))
    rw_rz_2d = rw_rz.reshape((1, len(rw_rz)))

iter1, iter2 = di['iter1'], di['iter2']
vals = di['vals']
lut = di['lut']

# Get the time range in sec
t1 = translate_times(iter1, dirname, translate_from='iter')['val_sec']
t2 = translate_times(iter2, dirname, translate_from='iter')['val_sec']

# Get the baseline time unit
rotation = get_parameter(dirname, 'rotation')
if rotation:
    time_unit = compute_Prot(dirname)
    time_label = r'$\rm{P_{rot}}$'
else:
    time_unit = compute_tdt(dirname)
    time_label = r'$\rm{TDT}$'

if plotdir is None:
    plotdir = dirname + '/plots/'
    if not os.path.isdir(plotdir):
        os.makedirs(plotdir)

# Get necessary grid info
#rr = di['rr']
cost = di['cost']
sint = di['sint']
tt_lat = di['tt_lat']
xx = di['xx']

# Get reference state stuff 
eq = get_eq(dirname)
rho = eq.density
T = eq.temperature
dlnT = eq.dlnT
dsdr = eq.dsdr
rhoT = rho*T
grav = eq.gravity

rho_2d = rho.reshape((1, nr))
T_2d = T.reshape((1, nr))
dsdr_2d = dsdr.reshape((1, nr))
rhoT_2d = rhoT.reshape((1, nr))
dlnT_2d = dlnT.reshape((1, nr))
grav_2d = grav.reshape((1, nr))

#========================
# KINETIC ENERGY EQUATION
#========================

# Advection work on KE
work_ke_advec = vals[:, :, lut[1910]]

# Pressure work on KE
work_pressure = vals[:, :, lut[1901]]

# Buoyancy work on KE 
# need this first
S_00 = np.sum(vals[:, :, lut[501]]*tw_2d, axis=0)
vr = vals[:, :, lut[1]]
vrS_00 = vr*S_00.reshape((1, nr))
try:
    work_buoy = vals[:, :, lut[1904]]
    print ("buoy_work = 1904 output in AZ_Avgs")
    # Will need this for "negligible work" and possibly buoyancy work
    vrS_lnot0 = c_P*work_buoy/(rho_2d*grav_2d)
    rhoTvrS = rhoT_2d*(vrS_lnot0 + vrS_00)
except: # hopefully they output 1440
    rhoTvrS = vals[:, :, lut[1440]]
    vrS = rhoTvrS/rhoT
    work_buoy = (rho*grav/c_P).reshape((1, nr))*(vrS - vrS_00)
    print ("buoy_work = 1904 not output in AZ_Avgs")
    print ("getting buoyancy work from 1440 = rhoTvrS")
    print ("subtracting l = 0 part using sum(tw*(501=S))")

# viscous work on kinetic energy
work_visc_on_ke = vals[:, :, lut[1907]]

# total ke work
work_ke = work_ke_advec + work_pressure + work_buoy + work_visc_on_ke

# J X B work on KE, if magnetic
if magnetism:
    work_mag = vals[:, :, lut[1916]]
    work_ke += work_mag

# Possibly work from forcing function
if forced:
    mean_vp = vals[:, :, lut[3]]
    vp2 = vals[:, :, lut[416]]
    fluc_vp2 = vals[:, :, lut[416]]
    mean_vp2 = vp2 - fluc_vp2
    tacho_r = get_parameter(dirname, 'tacho_r')
    print ("read tacho_r = %1.2e" %tacho_r)
    tacho_dr = get_parameter(dirname, 'tacho_dr')
    tacho_tau = get_parameter(dirname, 'tacho_tau')
    work_forcing = np.zeros((nt, nr))
    eq = get_eq(dirname)
    rho = eq.rho

    if os.path.exists(dirname + '/eq_vp'): 
        print ("eq_vp file exists, so I assume you have a forcing function which\n quartically matches on to a CZ differential rotation\n with viscous-torque-free buffer zone")
        tacho_r2 = get_parameter(dirname, 'tacho_r2')
        i_tacho_r = np.argmin(np.abs(rr - tacho_r))
        print ("read tacho_r2 = %1.2e" %tacho_r2)
        eq_vp = read_eq_vp(dirname + '/eq_vp', nt, nr)
        for it in range(nt):
            for ir in range(nr):
                if rr[ir] <= tacho_r2:
                    if rr[ir] > tacho_r:
                        # Here is where the DR is forced differentially
                        # (a "buffer zone" to reduce shear)
                        desired_vp = eq_vp[it, ir]
                    elif rr[ir] > tacho_r - tacho_dr*rr[0]:
                        # Here is where the DR is forced to match 
                        # quartically from differential to solid-body
                        desired_vp = eq_vp[it, i_tacho_r]*(1.0 - ( (rr[ir] - tacho_r)/(tacho_dr*rr[0]) )**2)**2
                    else:
                        desired_vp = 0.0
                    work_forcing[it, ir] = -rho[ir]*(mean_vp2[it, ir] -\
                            desired_vp*mean_vp[it, ir])/tacho_tau
                else:
                    work_forcing[it, ir] = 0.
    else:
        forcing_coeff = -rho/tacho_tau*0.5*(1.0 - np.tanh((rr - tacho_r)/(tacho_dr*rr[0])))
        work_forcing = forcing_coeff.reshape((1, nr))*mean_vp2
    work_ke += work_forcing

#========================
# HEAT (ENTROPY) EQUATION
#========================

# Advection of entropy
work_thermal_advec = -vals[:, :, lut[1401]] # 1401 = + rho*T*v dot grad S

# Stable gradient (advection of reference entropy)
vr = vals[:, :, lut[1]]
work_thermal_advec_ref = -rhoT_2d*dsdr_2d*vr

# heat by conduction
work_cond = vals[:, :, lut[1421]]

# heat from heating function (representing radiation)
work_rad = vals[:, :, lut[1434]] # Q(r)

# (irreversible) viscous heating
work_visc_on_inte = rhoT_2d*vals[:, :, lut[1435]] 
# 1435 is off by rho*T compared to what it's supposed to be
# (i.e., it's the viscous heating as in the entropy equation, not energy eq)

# total work on internal energy
work_inte = work_thermal_advec + work_thermal_advec_ref + work_cond +\
        work_rad + work_visc_on_inte

# If magnetism get Joule heating
if magnetism:
    work_joule = rhoT_2d*vals[:, :, lut[1436]]
    # joule heating also off by rho*T
    work_inte += work_joule

#=========================
# MAGNETIC ENERGY EQUATION
#=========================
if magnetism:
    work_induct = 1.0/(4.0*np.pi)*(vals[:, :, lut[2019]]) 
    work_idiff = 1.0/(4.0*np.pi)*vals[:, :, lut[2043]]
    # total magnetic energy work
    work_me = work_induct + work_idiff
    try: # decompose induction work further if possible if possible
        work_ishear = 1.0/(4.0*np.pi)*vals[:, :, lut[2025]]
        work_ishear_mmm = 1.0/(4.0*np.pi)*vals[:, :, lut[2034]]
        work_ishear_ppp = 1.0/(4.0*np.pi)*vals[:, :, lut[2040]]
        work_iadvec = 1.0/(4.0*np.pi)*vals[:, :, lut[2026]]
        work_icomp = 1.0/(4.0*np.pi)*vals[:, :, lut[2027]]
        have_all_induct_terms = True
        print ("I have all the magnetic induction work terms:")
        print ("shear (tot, mmm, ppp), advec., and comp.")
    except:
        have_all_induct_terms = False
        print ("Missing one or more induction work terms; include")
        print ("2025,2034,2040,2026,2027:")
        print ("[shear (tot, mmm, ppp), advec., and comp.] next time")

#======================================
# SPHERICALLY/VOLUME AVERAGED EQUATIONS
#======================================

# first collect all terms by category, making Latex titles 
# and plain-text labels along the way

# kinetic energy
ke_terms = [work_ke_advec, work_pressure, work_buoy, work_visc_on_ke,\
        work_ke]
ke_titles = [r'$\left\langle-\overline{\rho}\mathbf{u}\cdot\nabla\frac{u^2}{2}\right\rangle$', r'$-\nabla\cdot\langle P\mathbf{u}\rangle$', r'$\overline{\rho}g\left\langle u_r\frac{S}{c_p}\right\rangle$', r'$\langle\mathbf{u}\cdot(\nabla\cdot\mathbf{D})\rangle$', r'$\frac{\partial}{\partial t}\left\langle\frac{1}{2}\overline{\rho}u^2\right\rangle$']
ke_labels = ['ke adv', 'press', 'buoy', 'visc work', 'd(ke)/dt']

# internal energy
inte_terms = [work_thermal_advec, work_thermal_advec_ref, work_cond,\
        work_rad, work_visc_on_inte, work_inte]
inte_titles = [r'$\left\langle-\overline{\rho}\overline{T}\mathbf{u}\cdot\nabla S\right\rangle$', r'$-\overline{\rho}\overline{T}\frac{dS}{dr}\langle u_r\rangle$', r'$\nabla\cdot[\kappa\overline{\rho}\overline{T}\nabla\langle\ S\rangle]$', r'$Q(r)$', r'$\langle\mathbf{D}:\nabla\mathbf{u}\rangle$', r'$\overline{\rho}\overline{T}\frac{\partial\langle S\rangle}{\partial t}$']
inte_labels = ['S adv', 'ref adv', 'cond', 'rad', 'visc heat', 'd(inte)/dt']

# add more terms in case of forcing/magnetism
if forced:
    # kinetic energy
    ke_terms.insert(4, work_forcing)
    ke_titles.insert(4, r'$\langle\mathbf{v}\cdot\mathbf{f}\rangle$')
    ke_labels.insert(4, 'forcing')

    # internal energy
    if force_econs:
        inte_terms(5, -work_forcing)
        ke_titles.insert(5, r'$-\langle\mathbf{v}\cdot\mathbf{f}\rangle$')
        ke_labels.insert(5, 'forcing')

if magnetism:
    # kinetic energy
    ke_terms.insert(4, work_mag)
    ke_titles.insert(4, r'$\frac{1}{4\pi}\langle\mathbf{u}\cdot[(\nabla\times\mathbf{B})\times\mathbf{B}]\rangle$')
    ke_labels.insert(4, 'mag work')

    # internal energy
    inte_terms.insert(5, work_joule)
    inte_titles.insert(5, r'$\frac{\eta}{4\pi}\langle(\nabla\times\mathbf{B})^2\rangle$')
    inte_labels.insert(5, 'mag heat')

    # magnetic energy
    me_terms = [work_induct, work_idiff, work_me]
    me_titles = [r'$\frac{1}{4\pi}\langle\mathbf{B}\cdot\nabla\times(\mathbf{v}\times\mathbf{B})\rangle$', r'$\frac{1}{4\pi}\langle\mathbf{B}\cdot\nabla\times(\eta(r)\nabla\times\mathbf{B})\rangle$', r'$\frac{\partial}{\partial t}\left\langle\frac{B^2}{8\pi}\right\rangle$']
    me_labels = ['induct', 'idiff', 'd(me)/dt']
    if have_all_induct_terms:
        me_terms.insert(1, work_ishear)
        me_terms.insert(2, work_ishear_mmm)
        me_terms.insert(3, work_ishear_ppp)
        me_terms.insert(4, work_iadvec)
        me_terms.insert(5, work_icomp)
        me_titles.insert(1, r'$\frac{1}{4\pi}\langle\mathbf{B}\cdot[\mathbf{B}\cdot\nabla\mathbf{v}]\rangle$')
        me_titles.insert(2, r'$\frac{1}{4\pi}\langle\mathbf{B}\rangle\cdot[\langle mathbf{B}\rangle\cdot\nabla\langle\mathbf{v}\rangle]$')
        me_titles.insert(3, r'$\frac{1}{4\pi}\langle\mathbf{B}^\prime\cdot[\mathbf{B}^\prime\cdot\nabla\mathbf{v}^\prime]\rangle$')
        me_titles.insert(4, r'$-\frac{1}{4\pi}\langle\mathbf{B}\cdot[\mathbf{v}\cdot\nabla\mathbf{B}]\rangle$')
        me_titles.insert(5, r'$-\frac{1}{4\pi}\langle(\nabla\cdot\mathbf{v})B^2\rangle$')
        me_labels.insert(1, 'ishear')
        me_labels.insert(2, 'ishear (mmm)')
        me_labels.insert(3, 'ishear (ppp)')
        me_labels.insert(4, 'iadvec')
        me_labels.insert(5, 'icomp')

#====================================================
# TOTAL ENERGY EQUATION (negative divergences mostly)
#====================================================
work_visc = work_visc_on_inte + work_visc_on_ke
work_enth = work_pressure + work_thermal_advec
# Negligible work (probably second-order but still worried about it)
work_negligible = rhoTvrS*dsdr_2d/c_P
# "extra" term for T_ref to go inside convective derivative
work_extra_advec = rhoTvrS*dlnT_2d 
work_inte_advec = work_thermal_advec - work_extra_advec
work_tot = work_ke + work_inte

tote_terms = [work_ke_advec, work_enth, work_cond, work_rad, work_visc,\
        work_thermal_advec_ref, work_buoy, work_extra_advec,\
        work_negligible, work_tot]
tote_titles =\
[r'$-\nabla\cdot\left\langle\frac{1}{2}\overline{\rho}u^2\right\rangle$',\
r'$-\nabla\cdot\langle (\overline{\rho}\overline{T}S + P)\mathbf{u}\rangle$',\
r'$\nabla\cdot[\kappa\overline{\rho}\overline{T}\nabla\langle\ S \rangle]$',\
r'$Q(r)$',\
r'$\nabla\cdot\langle\mathbf{D}\cdot\mathbf{u}\rangle$',\
r'$-\overline{\rho}\overline{T}\frac{d\overline{S}}{dr}\langle u_r\rangle$',\
r'$\overline{\rho}g\left\langle u_r\frac{S}{c_p}\right\rangle$',\
r'$\overline{\rho}\frac{d\overline{T}}{dr}\langle u_r S \rangle$',\
r'$\overline{\rho}\overline{T}\frac{d\overline{S}}{dr}\left\langle u_r\frac{S}{c_p}\right\rangle$',\
r'$\frac{\partial}{\partial t}\left\langle TOTE\right\rangle$']
tote_labels = ['ke adv', 'enth', 'cond', 'rad', 'visc', 'ref adv', 'buoy',\
        'extra', 'small', 'tot']

if magnetism:
    work_Poyn = work_mag + work_induct
    tote_terms[-1] += work_me
    tote_terms.insert(9, work_Poyn)
    tote_titles.insert(9, r'$-\frac{1}{4\pi}\nabla\cdot\left\langle[\eta\nabla\times\mathbf{B}-\mathbf{u}\times\mathbf{B}]\times\mathbf{B}\right\rangle$')
    tote_labels.insert(9, 'Poyn')

# now collect all terms for plotting
all_terms = {'toten_tote': tote_terms, 'toten_ke': ke_terms, 'toten_inte': inte_terms}
all_titles = {'toten_tote': tote_titles, 'toten_ke': ke_titles, 'toten_inte': inte_titles}
all_labels = {'toten_tote': tote_labels, 'toten_ke': ke_labels, 'toten_inte': inte_labels}
plot_labels = {'toten_tote': 'total energy equation', 'toten_ke': 'kinetic energy equation', 'toten_inte': 'heat equation'}

if magnetism:
    all_terms['toten_me'] = me_terms
    all_titles['toten_me'] = me_titles
    all_labels['toten_me'] = me_labels
    plot_labels['toten_me'] = 'magnetic energy equation'

# general figure parameters for AZ plot
fs = 12.
small_fs = 9.
lw = 1.
units = r'$\rm{erg}\ \rm{cm}^{-3}\ \rm{s}^{-1}$'

# loop through keys and make plots
for key in keys:
    # terms to plot
    print ("plotting ", key, " terms")
    terms = all_terms[key]
    nplots = len(terms)
    titles = all_titles[key]
    labels = all_labels[key]

    # set up AZ_Avg figure axes
    ncol = 4
    nrow = np.int(np.ceil(nplots/ncol))
    fig_width_inches = 8.5 # 8.5 x 11 paper
    # default margin
    margin_inches = 1./8.
     # wider top margin to accommodate metadata   
    margin_top_inches = 1.
    # margin to accommodate just subplot titles   
    margin_subplot_top_inches = 1./4.
    # larger subplot bottom margin to make room for colorbar(s)
    margin_subplot_bottom_inches = 0.75*(2 - (rbcz is None)) 

    subplot_width_inches = (fig_width_inches - (ncol+1)*margin_inches)/ncol
    # aspect ratio for azavg plot = 2 x 1
    subplot_height_inches = 2*subplot_width_inches
    fig_height_inches = margin_top_inches + nrow*(subplot_height_inches +\
            margin_subplot_top_inches + margin_subplot_bottom_inches)

    # unitless dimensions
    margin_x = margin_inches/fig_width_inches
    margin_y = margin_inches/fig_height_inches
    margin_top = margin_top_inches/fig_height_inches
    margin_subplot_top = margin_subplot_top_inches/fig_height_inches
    margin_subplot_bottom = margin_subplot_bottom_inches/fig_height_inches
    subplot_width = subplot_width_inches/fig_width_inches
    subplot_height = subplot_height_inches/fig_height_inches

    # Generate figure for AZ Avgs plot
    fig = plt.figure(figsize=(fig_width_inches, fig_height_inches))

    # 1D averages figure
    av_nrow = 2 + 2*(not rbcz is None)
    av_fig_width_inches = 8.5 
    av_width_inches = 5.
    av_height_inches = 3.5
    av_margin_inches = 3./4.
    av_margin_top_inches = 1.
    #av_margin_subplot_top_inches = 3./4.
    av_fig_height_inches = av_margin_top_inches +\
            av_nrow*(av_height_inches + av_margin_inches)
    av_width = av_width_inches/av_fig_width_inches
    av_height = av_height_inches/av_fig_height_inches
    av_margin_x = av_margin_inches/av_fig_width_inches
    av_margin_y = av_margin_inches/av_fig_height_inches
    av_margin_top = av_margin_top_inches/av_fig_height_inches
    fig_av = plt.figure(figsize=(av_fig_width_inches,\
            av_fig_height_inches))
    ax_shav = fig_av.add_axes((av_margin_x, 1. - av_margin_top -\
            av_height, av_width, av_height))
    ax_rav = fig_av.add_axes((av_margin_x, 1. - av_margin_top -\
            av_height - 1.*(av_height + av_margin_y), av_width, av_height))
    # collect all av axes
    axes_shav = [ax_shav]
    axes_rav = [ax_rav]
    if not rbcz is None:
        ax_shav_rz = ax_shav.twinx()
        axes_shav.append(ax_shav_rz)
        ax_rav_cz = fig_av.add_axes((av_margin_x, 1. - av_margin_top -\
                av_height - 2.*(av_height + av_margin_y), av_width,\
                av_height))
        axes_rav.append(ax_rav_cz)
        ax_rav_rz = fig_av.add_axes((av_margin_x, 1. - av_margin_top -\
                av_height - 3.*(av_height + av_margin_y), av_width,\
                av_height))
        axes_rav.append(ax_rav_rz)

    # loop over terms and add associated plot to axes
    for iplot in range(nplots):
        # compute relevant terms
        term = terms[iplot]
        shav_term = np.sum(term*tw_2d, axis=0)
        integrated_term = volume*np.sum(shav_term*rw)
        rav_term = np.sum(term*rw_2d, axis=1)
        if not rbcz is None:
            rav_term_cz = np.sum(term[:, :irbcz]*rw_cz_2d, axis=1)
            rav_term_rz = np.sum(term[:, irbcz:]*rw_rz_2d, axis=1)
            integrated_term_cz = volume_cz*np.sum(shav_term[:irbcz]*rw_cz)
            integrated_term_rz = volume_rz*np.sum(shav_term[irbcz:]*rw_rz)

        # plot azav
        ax_left = margin_x + (iplot%ncol)*(subplot_width + margin_x)
        ax_bottom = 1 - margin_top - margin_subplot_top - subplot_height -\
                (iplot//ncol)*(margin_subplot_top + subplot_height +\
                margin_subplot_bottom)
        ax = fig.add_axes((ax_left, ax_bottom, subplot_width,\
                subplot_height))
        plot_azav (term, rr, cost, fig=fig, ax=ax, units=units,\
               minmax=minmax, plotcontours=plotcontours, rvals=rvals,\
               minmaxrz=minmaxrz, rbcz=rbcz, symlog=symlog,\
        linthresh=linthresh, linscale=linscale, linthreshrz=linthreshrz,\
        linscalerz=linscalerz, plotlatlines=plotlatlines)
        ax.set_title(titles[iplot], verticalalignment='bottom', **csfont)

        # plot shav
        shav_label = labels[iplot] + ': ' +\
                sci_format(integrated_term/lstar, 3)
        if rbcz is None:
            ax_shav.plot(rr/rsun, shav_term, label=shav_label, linewidth=lw)
            ax_rav.plot(tt_lat, rav_term, linewidth=lw)
        else:
            ax_shav.plot(rr[:irbcz]/rsun, shav_term[:irbcz],\
                    label=shav_label, linewidth=lw)
            ax_shav_rz.plot(rr[irbcz:]/rsun, shav_term[irbcz:],\
                    label=shav_label, linewidth=lw)
            ax_rav.plot(tt_lat, rav_term, linewidth=lw)
            rav_cz_label = labels[iplot] + ': ' +\
                    sci_format(integrated_term_cz/lstar, 3)
            rav_rz_label = labels[iplot] + ': ' +\
                    sci_format(integrated_term_rz/lstar, 3)
            ax_rav_cz.plot(tt_lat, rav_term_cz, label=rav_cz_label,\
                    linewidth=lw)
            ax_rav_rz.plot(tt_lat, rav_term_rz, label=rav_rz_label,\
                    linewidth=lw)

    # label av plot x axes
    ax_shav.set_xlabel(r'$r/R_\odot$', fontsize=fs, **csfont)
    for ax in axes_rav:
        ax.set_xlabel('latitude (deg)', fontsize=fs, **csfont)
    # make av plot legends
    ax_shav.legend(bbox_to_anchor=(1.05, 1), loc=2, fontsize=small_fs,\
            labelspacing=0.5, title='(volume integral)/L')
    if not rbcz is None:
        ax_rav_cz.legend(bbox_to_anchor=(1.05, 1), loc=2, fontsize=small_fs,            labelspacing=0.5, title='(vol. integral CZ)/L')
        ax_rav_rz.legend(bbox_to_anchor=(1.05, 1), loc=2, fontsize=small_fs,            labelspacing=0.5, title='(vol. integral RZ)/L')
    # make av plot titles
    ax_shav.set_title("spherical average", fontsize=fs, **csfont )
    ax_rav.set_title("radial average", fontsize=fs, **csfont )
    if not rbcz is None:
        ax_rav_cz.set_title("radial avg. CZ", fontsize=fs, **csfont)
        ax_rav_rz.set_title("radial avg. RZ", fontsize=fs, **csfont)
    # if rbcz is desired centralize the zero point of the y axis
        for ax in [ax_shav, ax_shav_rz]:
            ymin, ymax = ax.get_ylim()
            ymaxabs = max(np.abs(ymin), np.abs(ymax))
            ax.set_ylim(-ymaxabs, ymaxabs)
    # mark desired radii on shav plot
    for rval in rvals:
        if rbcz is None:
            ax = ax_shav
        else:
            if rval <= rbcz:
                ax = ax_shav_rz
            else:
                ax = ax_shav
        ymin, ymax = ax.get_ylim()
        yvals = np.linspace(ymin, ymax, 100)
        ax.plot(np.zeros(100) + rval/rsun, yvals, 'k--', linewidth=0.5*lw)
    # mark zero points on all plots
    for ax in axes_shav + axes_rav:
        xmin, xmax = ax.get_xlim()
        xvals = np.linspace(xmin, xmax, 100)
        ax.plot(xvals, np.zeros(100), 'k--', linewidth=0.5*lw)

    # shav ticks (mostly everywhere, deal with split axes)
    if rbcz is None:
        plt.sca(ax_shav)
        plt.minorticks_on()
        plt.tick_params(top=True, right=True, direction='in',\
                which='both')
    else:
        plt.sca(ax_shav)
        plt.minorticks_on()
        plt.tick_params(top=True, left=False, right=True, direction='in',\
                which='both')
        ax_shav.yaxis.tick_right()
        plt.sca(ax_shav_rz)
        plt.minorticks_on()
        plt.tick_params(top=True, left=True, right=False, direction='in',\
                which='both')
        ax_shav_rz.yaxis.tick_left()
    # rav ticks (everywhere)
    for ax in axes_rav:
        plt.sca(ax)
        ax.set_xlabel('latitude (deg)', fontsize=fs, **csfont)
        plt.minorticks_on()
        plt.tick_params(top=True, right=True, direction='in', which='both')

    # put averaging interval in metadata
    if rotation:
        time_string = ('t = %.1f to %.1f ' %(t1/time_unit, t2/time_unit))\
                + time_label + (r'$\ (\Delta t = %.1f\ $'\
                %((t2 - t1)/time_unit)) + time_label + ')'
    else:
        time_string = ('t = %.3f to %.3f ' %(t1/time_unit, t2/time_unit))\
                + time_label + (r'$\ (\Delta t = %.3f\ $'\
                %((t2 - t1)/time_unit)) + time_label + ')'

    # metadata for azav
    fig.text(margin_x, 1 - 0.1*margin_top, dirname_stripped,\
             ha='left', va='top', fontsize=fs, **csfont)
    fig.text(margin_x, 1 - 0.3*margin_top, plot_labels[key] +\
            ' (zonal average)', ha='left', va='top',\
            fontsize=fs, **csfont)
    fig.text(margin_x, 1 - 0.5*margin_top, time_string,\
             ha='left', va='top', fontsize=fs, **csfont)

    # same metadata for av plot
    fig_av.text(av_margin_x, 1 - 0.1*av_margin_top, dirname_stripped,\
             ha='left', va='top', fontsize=fs, **csfont)
    fig_av.text(av_margin_x, 1 - 0.3*av_margin_top,\
            plot_labels[key] + ' (1D average)', ha='left', va='top',\
            fontsize=fs, **csfont)
    fig_av.text(av_margin_x, 1 - 0.5*av_margin_top, time_string,\
             ha='left', va='top', fontsize=fs, **csfont)

    # save both plots
    # azav
    savefile = plotdir + dirname_stripped + '_' + key + '_' +\
            str(iter1).zfill(8) + '_' + str(iter2).zfill(8) + tag + '.png'
    print ('Saving azav plot at\n' + savefile)
    fig.savefig(savefile, dpi=300)
    plt.close(fig)
    # av
    savefile_av = plotdir + dirname_stripped + '_' + key + '_av_' +\
            str(iter1).zfill(8) + '_' + str(iter2).zfill(8) + tag + av_ext
    print ('Saving 1D av line plot at\n' + savefile_av)
    fig_av.savefig(savefile_av)
    plt.close(fig_av)
