# Author: Loren Matilsky
# Created: 12/30/2019
# This script plots each term in the total energy equation in the 
# meridional plane,
# each term's latitudinally (spherically) integrated profile vs. radius,
# and the total integrated rate of change of the energy.
# Plots the total energy equation, KE eq., ME eq., and heat eq. separately
# to choose only one plot, specify
# -tot, -ke, -inte, -me
# Analyzes Rayleigh run directory indicated by [dirname]. To use an
# AZ_Avgs file different than the one associated with the longest 
# averaging range, use
# -usefile [complete name of desired AZ_Avgs file]
# Saves plot in
# [dirname]_toten_[first iter]_[last iter].png

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
from common import get_widest_range_file, strip_dirname, get_dict, rsun,\
        c_P, get_lum, sci_format
from get_parameter import get_parameter
from rayleigh_diagnostics import GridInfo
from get_eq import get_eq
from time_scales import compute_Prot, compute_tdt
from translate_times import translate_times
from read_inner_vp import read_inner_vp
from read_eq_vp import read_eq_vp

# Get directory name and stripped_dirname for plotting purposes
dirname = sys.argv[1]
dirname_stripped = strip_dirname(dirname)

# Directory with data and plots, make the plotting directory if it doesn't
# already exist    
datadir = dirname + '/data/'
plotdir = dirname + '/plots/'
if not os.path.isdir(plotdir):
    os.makedirs(plotdir)

# Read command-line arguments (CLAs)
magnetism = None
showplot = False
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
rvals = None
rbcz = None
symlog = False
crudeint = False # by default use exact latitudinal weights and instead 
# if crudeint = True, use a "crude weight"
tag = ''

# determine magnetism from main_input
magnetism = get_parameter(dirname, 'magnetism')
keys = ['toten_tote', 'toten_ke', 'toten_inte']
if magnetism:
    keys.append('toten_me')

# need stellar luminosity
lstar = get_lum(dirname)

args = sys.argv[2:]
nargs = len(args)
for i in range(nargs):
    arg = args[i]
    if arg == '-minmax':
        minmax = float(args[i+1]), float(args[i+2])
    elif arg == '-minmaxrz':
        minmaxrz = float(args[i+1]), float(args[i+2])
    elif arg == '-rbcz':
        rbcz = float(args[i+1])
    elif arg == '-show':
        showplot = True
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
        rvals = []
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
        showplot = True
    elif arg == '-inte': # just heat equation
        keys = ['toten_inte']
        showplot = True
    elif arg == '-me': # just ME equation
        keys = ['toten_me']
        showplot = True
    elif arg == '-tote':
        keys = ['toten_tote']
        showplot = True

# Get the data:
print ('Getting energy production terms from ' +\
        datadir + AZ_Avgs_file)
di = get_dict(datadir + AZ_Avgs_file)

# Will need integration weights for a few things
gi = GridInfo(dirname + '/grid_info')
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

rw_2d = rw.reshape((1, nr))
tw_2d = tw.reshape((nt, 1))

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

# Get necessary grid info
rr = di['rr']
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
    work_buoy = (rho*grav/c_P).rshape((1, nr))*(vrS - vrS_00)
    print ("buoy_work = 1904 not output in AZ_Avgs")
    print ("getting buoyancy work from 1440 = rhoTvrS")
    print ("subtracting l = 0 part using sum (tw * (501 = S))")

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

# compute the spherical averages and integrated terms

# shell volume
ri, ro = np.min(rr), np.max(rr)
volume = 4./3.*np.pi*(ro**3. - ri**3.)

# kinetic energy
ke_shav_terms = []
ke_integrated_terms = []
for term in ke_terms:
    shav_term = np.sum(term*tw_2d, axis=0)
    integrated_term = volume*np.sum(shav_term*rw)
    ke_shav_terms.append(shav_term)
    ke_integrated_terms.append(integrated_term)

# internal energy
inte_shav_terms = []
inte_integrated_terms = []
for term in inte_terms:
    shav_term = np.sum(term*tw_2d, axis=0)
    integrated_term = volume*np.sum(shav_term*rw)
    inte_shav_terms.append(shav_term)
    inte_integrated_terms.append(integrated_term)

# magnetic energy
if magnetism:
    me_shav_terms = []
    me_integrated_terms = []
    for term in inte_terms:
        shav_term = np.sum(term*tw_2d, axis=0)
        integrated_term = volume*np.sum(shav_term*rw)
        me_shav_terms.append(shav_term)
        me_integrated_terms.append(integrated_term)

#====================================================
# TOTAL ENERGY EQUATION (negative divergences mostly)
#====================================================
work_visc = work_visc_on_inte + work_visc_on_ke
work_enth = work_pressure + work_thermal_advec
# Negligible work (probably second-order but still worried about it)
work_negligible = rhoTvrS*dsdr_2d/c_P
# "extra" term for T_ref togoes inside convective derivative
work_extra_advec = -rhoTvrS*dlnT_2d 
work_inte_advec = work_thermal_advec + work_extra_advec
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
r'$-\overline{\rho}\frac{d\overline{T}}{dr}\langle u_r S \rangle$',\
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

# spherically and volume averaged terms
tote_shav_terms = []
tote_integrated_terms = []
for term in tote_terms:
    shav_term = np.sum(term*tw_2d, axis=0)
    integrated_term = volume*np.sum(shav_term*rw)
    tote_shav_terms.append(shav_term)
    tote_integrated_terms.append(integrated_term)

# now collect all terms for plotting
all_terms = {'toten_tote': tote_terms, 'toten_ke': ke_terms, 'toten_inte': inte_terms}
all_shav_terms = {'toten_tote': tote_shav_terms, 'toten_ke': ke_shav_terms, 'toten_inte': inte_shav_terms}
all_integrated_terms = {'toten_tote': tote_integrated_terms, 'toten_ke': ke_integrated_terms, 'toten_inte': inte_integrated_terms}
all_titles = {'toten_tote': tote_titles, 'toten_ke': ke_titles, 'toten_inte': inte_titles}
all_labels = {'toten_tote': tote_labels, 'toten_ke': ke_labels, 'toten_inte': inte_labels}
plot_labels = {'toten_tote': 'total energy equation', 'toten_ke': 'kinetic energy equation', 'toten_inte': 'heat equation'}

if magnetism:
    all_terms['toten_me'] = me_terms
    all_shav_terms['toten_me'] = me_shav_terms
    all_integrated_terms['toten_me'] = me_integrated_terms
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
    print ("plotting ", key)
    terms = all_terms[key]
    nplots = len(terms)
    shav_terms = all_shav_terms[key]
    integrated_terms = all_integrated_terms[key]
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

    # spherical average figure
    shav_nrow = 1
    shav_fig_width_inches = 8.5 
    shav_width_inches = 5.
    shav_height_inches = 3.5
    shav_margin_inches = 0.6
    shav_margin_top_inches = 1.
    shav_fig_height_inches = shav_margin_top_inches +\
            shav_nrow*(shav_height_inches + shav_margin_inches)
    shav_width = shav_width_inches/shav_fig_width_inches
    shav_height = shav_height_inches/shav_fig_height_inches
    shav_margin_x = shav_margin_inches/shav_fig_width_inches
    shav_margin_y = shav_margin_inches/shav_fig_height_inches
    shav_margin_top = shav_margin_top_inches/shav_fig_height_inches
    fig_shav = plt.figure(figsize=(shav_fig_width_inches,\
            shav_fig_height_inches))
    ax_shav = fig_shav.add_axes((shav_margin_x, shav_margin_y, shav_width,\
            shav_height))
    if not rbcz is None:
        ax_shav_rz = ax_shav.twinx()
        irbcz = np.argmin(np.abs(rr - rbcz))

    for iplot in range(nplots):
        ax_left = margin_x + (iplot%ncol)*(subplot_width + margin_x)
        ax_bottom = 1 - margin_top - margin_subplot_top - subplot_height -\
                (iplot//ncol)*(margin_subplot_top + subplot_height +\
                margin_subplot_bottom)
        ax = fig.add_axes((ax_left, ax_bottom, subplot_width,\
                subplot_height))
        plot_azav (terms[iplot], rr, cost, fig=fig, ax=ax, units=units,\
               minmax=minmax, plotcontours=plotcontours, rvals=rvals,\
               minmaxrz=minmaxrz, rbcz=rbcz, symlog=symlog,\
        linthresh=linthresh, linscale=linscale, linthreshrz=linthreshrz,\
        linscalerz=linscalerz, plotlatlines=plotlatlines)
        ax.set_title(titles[iplot], verticalalignment='bottom', **csfont)

        # plot spherically averaged work
        shav_label = labels[iplot] + ': ' +\
                sci_format(integrated_terms[iplot]/lstar, 3) + r'$L_*$'
        if rbcz is None:
            ax_shav.plot(rr/rsun, shav_terms[iplot], label=shav_label,\
                linewidth=lw)
        else:
            ax_shav.plot(rr[:irbcz]/rsun, shav_terms[iplot][:irbcz],\
                    label=shav_label, linewidth=lw)
            ax_shav_rz.plot(rr[irbcz:]/rsun, shav_terms[iplot][irbcz:],\
                    label=shav_label, linewidth=lw)

    # fix up some stuff for the shav work line plot
    ax_shav.set_xlabel(r'$r/R_\odot$', fontsize=fs, **csfont)
    #ax_shav.set_xlim((ri/rsun, ro/rsun))
    ax_shav.legend(bbox_to_anchor=(1.05, 1), loc=2, fontsize=small_fs,\
            labelspacing=0.5, title='volume integral')
    ax_shav.set_title("cgs units: erg/cm^3", fontsize=fs,\
            **csfont )
    # mark zero line
    ax_shav.plot(rr/rsun, np.zeros_like(rr), 'k--', linewidth=0.5*lw)
    # ticks (mostly) everywhere
    if rbcz is None:
        plt.sca(ax_shav)
        plt.minorticks_on()
        plt.tick_params(top=True, right=True, direction='in', which='both')
    else:
        plt.sca(ax_shav)
        plt.minorticks_on()
        plt.tick_params(top=True, left=False, right=True, direction='in',\
                which='both')
        #ax_shav.yaxis.set_label_position('right')
        ax_shav.yaxis.tick_right()
        plt.sca(ax_shav_rz)
        plt.minorticks_on()
        plt.tick_params(top=True, left=True, right=False, direction='in',\
                which='both')
        ax_shav_rz.yaxis.tick_left()
        #ax_shav_rz.set_xlim((ri/rsun, ro/rsun))
        # centralize the zero point
        for ax in [ax_shav, ax_shav_rz]:
            ymin, ymax = ax.get_ylim()
            ymaxabs = max(np.abs(ymin), np.abs(ymax))
            ax.set_ylim(-ymaxabs, ymaxabs)
        # mark the chosen rbcz
        ax_shav_rz.plot(np.zeros(100) + rr[irbcz]/rsun,\
                np.linspace(-ymaxabs, ymaxabs, 100), 'k--',\
                linewidth=0.5*lw)
        ax_shav_rz.plot(rr/rsun, np.zeros_like(rr), 'k--', linewidth=0.5*lw)

    # Label averaging interval
    if rotation:
        time_string = ('t = %.1f to %.1f ' %(t1/time_unit, t2/time_unit))\
                + time_label + (r'$\ (\Delta t = %.1f\ $'\
                %((t2 - t1)/time_unit)) + time_label + ')'
    else:
        time_string = ('t = %.3f to %.3f ' %(t1/time_unit, t2/time_unit))\
                + time_label + (r'$\ (\Delta t = %.3f\ $'\
                %((t2 - t1)/time_unit)) + time_label + ')'

    # Put some metadata in upper left of AZ_Avgs plot
    fig.text(margin_x, 1 - 0.1*margin_top, dirname_stripped,\
             ha='left', va='top', fontsize=fs, **csfont)
    fig.text(margin_x, 1 - 0.3*margin_top, plot_labels[key] +\
            ' (zonal average)', ha='left', va='top',\
            fontsize=fs, **csfont)
    fig.text(margin_x, 1 - 0.5*margin_top, time_string,\
             ha='left', va='top', fontsize=fs, **csfont)

    # ...and Sph Avgs plot
    fig_shav.text(shav_margin_x, 1 - 0.1*shav_margin_top, dirname_stripped,\
             ha='left', va='top', fontsize=fs, **csfont)
    fig_shav.text(shav_margin_x, 1 - 0.3*shav_margin_top,\
            plot_labels[key] + ' (spherical average)', ha='left', va='top',\
            fontsize=fs, **csfont)
    fig_shav.text(shav_margin_x, 1 - 0.5*shav_margin_top, time_string,\
             ha='left', va='top', fontsize=fs, **csfont)

    savefile = plotdir + dirname_stripped + '_' + key + '_' +\
            str(iter1).zfill(8) + '_' + str(iter2).zfill(8) + tag + '.png'

    print ('Saving AZ_Avgs plot at ' + savefile)
    fig.savefig(savefile, dpi=300)
    plt.close(fig)

    savefile_shav = plotdir + dirname_stripped + '_' + key + '_shav_' +\
            str(iter1).zfill(8) + '_' + str(iter2).zfill(8) + tag + '.pdf'

    print ('Saving spherically averaged plot at ' + savefile_shav)
    fig_shav.savefig(savefile_shav)
    plt.close(fig_shav)
