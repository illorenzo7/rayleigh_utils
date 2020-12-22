# Author: Loren Matilsky
# Created: 12/30/2019
# This script plots each term in the total energy equation in the 
# meridional plane,
# each term's latitudinally (spherically) integrated profile vs. radius,
# and the total integrated rate of change of the energy
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
showplot = True
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
    elif arg == '-noshow':
        showplot = False
    elif arg == '-nosave':
        saveplot = False
    elif arg == '-nocontour':
        plotcontours = False
    elif arg == '-usefile':
        AZ_Avgs_file = args[i+1]
        AZ_Avgs_file = AZ_Avgs_file.split('/')[-1]
    elif arg == '-forced':
        forced = True
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
    elif arg == '-mag':
        magnetism = bool(args[i+1])
    elif arg == '-crudeint':
        crudeint = True
    elif arg == '-tag':
        tag = '_' + args[i+1]

# by default determine magnetism from main_input
if magnetism is None:
    magnetism = get_parameter(dirname, 'magnetism')

# Get the data:
print ('Getting total energy production terms from ' +\
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
dsdr = eq.dsdr
rhoT = rho*T
grav = eq.gravity

rho_2d = rho.reshape((1, nr))
T_2d = T.reshape((1, nr))
dsdr_2d = dsdr.reshape((1, nr))
rhoT_2d = rhoT.reshape((1, nr))

# Calculate the negative divergences

#========================
# KINETIC ENERGY EQUATION
#========================

# Advection work on KE
work_ke_advec = vals[:, :, lut[1910]] #

# Pressure work on KE
work_pressure = vals[:, :, lut[1901]]

# Will need this for "negligible work" and possibly buoyancy work
rhoTvrS = vals[:, :, lut[1440]]

# Buoyancy work on KE (in future calculate 1904)
try:
    work_buoy = vals[:, :, lut[1904]]
    print ("buoy_work = 1904 output in AZ_Avgs")
except:
    vrS = rhoTvrS/rhoT
    S_00 = np.sum(vals[:, :, lut[501]]*tw_2d, axis=0)
    vr = vals[:, :, lut[1]]
    vrS_00 = vr*S_00.reshape((1, nr))
    work_buoy = (rho*grav/c_P).rshape((1, nr))*(vrS - vrS_00)
    print ("buoy_work = 1904 not output in AZ_Avgs")
    print ("getting buoyancy work from 1440 = rhoTvrS")
    print ("subtracting l = 0 part using sum (tw * (501 = S))")

# Viscous work on kinetic energy
work_visc_on_ke = vals[:, :, lut[1907]]

# total ke work
work_ke = work_ke_advec + work_pressure + work_buoy + work_visc_on_ke

# J X B work on KE, if magnetic
if magnetism:
    work_jcrossb = vals[:, :, lut[1916]]
    work_ke += work_jcrossb

# Possibly work from forcing function
if forced:
    force_econs = False # If False, this work will need to get 
                # included in total
                # If True, energy is conserved, so term in heat equation
                # balances out the work term
    try:
        force_econs = get_parameter(dirname, 'force_econs')
    except:
        pass
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
work_thermal_advec_ref = -rhoT_2d*vr*dsdr_2d

# heat by conduction
work_cond = vals[:, :, lut[1421]]

# heat from heating function (representing radiation)
work_rad = vals[:, :, lut[1434]] # Q(r)

# Irreversible (viscous) heating
work_visc_on_inte = rhoT_2d*vals[:, :, lut[1435]] 
# (irreversible) viscous heating
# 1435 THIS IS OFF BY rho*T compared to what it's supposed to be
# (i.e., it's the viscous heating as in the entropy equation, not energy eq)

# total work on internal energy
work_inte = work_thermal_advec + work_thermal_advec_ref + work_cond +\
        work_rad + work_visc_on_inte

# If magnetism get Joule heating
if magnetism:
    work_joule = rhoT_2d*vals[:, :, lut[1436]]
    # joule heating also off by rho*T
    work_inte += work_joule

# total work
work_tot = work_ke + work_inte

# Not sure actually if I should plot this now...
work_dsdr_negligible = rhoTvrS*dsdr_2d/c_P
#print ("maxabs (small dsdr work) = ", np.max(np.abs(work_dsdr_negligible)))

# magnetic energy equation
if magnetism:
    work_induct = 1.0/(4.0*np.pi)*(vals[:, :, lut[2019]]) 
    work_idiff = 1.0/(4.0*np.pi)*vals[:, :, lut[2043]]
    # total magnetic energy work
    work_me = work_induct + work_idiff
    # total work
    work_tot += work_me

# get spherically averaged work

# kinetic energy equation
work_ke_advec_r = np.sum(work_ke_advec*tw_2d, axis=0)
work_pressure_r = np.sum(work_pressure*tw_2d, axis=0)
work_buoy_r = np.sum(work_buoy*tw_2d, axis=0)
work_visc_on_ke_r = np.sum(work_visc_on_ke*tw_2d, axis=0)
if magnetism:
    work_jcrossb_r = np.sum(work_jcrossb*tw_2d, axis=0)
work_ke_r = np.sum(work_ke*tw_2d, axis=0)

# heat equation
work_thermal_advec_r = np.sum(work_thermal_advec*tw_2d, axis=0)
work_thermal_advec_ref_r = np.sum(work_thermal_advec_ref*tw_2d, axis=0)
work_cond_r = np.sum(work_cond*tw_2d, axis=0)
work_rad_r = np.sum(work_rad*tw_2d, axis=0)
work_visc_on_inte_r = np.sum(work_visc_on_inte*tw_2d, axis=0)
if magnetism:
    work_joule_r = np.sum(work_joule*tw_2d, axis=0)
if forced:
    work_forcing_r = np.sum(work_forcing*tw_2d, axis=0)
# not sure where this should go...
work_dsdr_negligible_r = np.sum(work_dsdr_negligible*tw_2d, axis=0)
work_inte_r = np.sum(work_inte*tw_2d, axis=0)

# total spherically averaged work
work_tot_r = work_ke_r + work_inte_r

# magnetic energy equation
if magnetism:
    work_induct_r = np.sum(work_induct*tw_2d, axis=0)
    work_idiff_r = np.sum(work_idiff*tw_2d, axis=0)
    work_me_r = np.sum(work_me*tw_2d, axis=0)
    work_tot_r += work_me_r

# get volume-integrated work (erg/s)
ri, ro = np.min(rr), np.max(rr)
volume = 4./3.*np.pi*(ro**3. - ri**3.)

# kinetic energy 
integrated_ke_advec = volume*np.sum(work_ke_advec_r*rw)
integrated_pressure = volume*np.sum(work_pressure_r*rw)
integrated_buoy = volume*np.sum(work_buoy_r*rw)
integrated_visc_on_ke = volume*np.sum(work_visc_on_ke_r*rw)
if magnetism:
    integrated_jcrossb = volume*np.sum(work_jcrossb_r*rw)
integrated_ke = volume*np.sum(work_ke_r*rw)

# internal energy
integrated_thermal_advec = volume*np.sum(work_thermal_advec_r*rw)
integrated_thermal_advec_ref = volume*np.sum(work_thermal_advec_ref_r*rw)
integrated_cond = volume*np.sum(work_cond_r*rw)
integrated_rad = volume*np.sum(work_rad_r*rw)
integrated_visc_on_inte = volume*np.sum(work_visc_on_inte_r*rw)
if magnetism:
    integrated_joule = volume*np.sum(work_joule_r*rw)
integrated_inte = volume*np.sum(work_inte_r*rw)

# total volume-integrated work
integrated_tot = integrated_ke + integrated_inte

# magnetic energy
if magnetism:
    integrated_induct = volume*np.sum(work_induct_r*rw)
    integrated_idiff = volume*np.sum(work_idiff_r*rw)
    integrated_me = volume*np.sum(work_me_r*rw)
    integrated_tot += integrated_me

# collect the terms to plot
work_terms = [work_ke_advec, work_pressure, work_buoy, work_visc_on_ke,\
        work_ke]
shav_work_terms = [work_ke_advec_r, work_pressure_r, work_buoy_r,\
        work_visc_on_ke_r, work_ke_r]
integrated_terms = [integrated_ke_advec, integrated_pressure,\
        integrated_buoy, integrated_visc_on_ke, integrated_ke]

titles = [r'$\left\langle-\overline{\rho}\mathbf{u}\cdot\nabla\frac{u^2}{2}\right\rangle$', r'$-\nabla\cdot\langle P\mathbf{u}\rangle$', r'$\overline{\rho}g\left\langle u_r\frac{S}{c_p}\right\rangle$', r'$\langle\mathbf{u}\cdot(\nabla\cdot\mathbf{D})\rangle$', r'$\frac{\partial}{\partial t}\left\langle\frac{1}{2}\overline{\rho}u^2\right\rangle$']
simple_labels = ['ke_advec', 'pressure', 'buoy', 'visc_on_ke', 'tot ke']

units = r'$\rm{erg}\ \rm{cm}^{-3}\ \rm{s}^{-1}$'

if forced:
    work_terms.insert(4, work_forcing)
    shav_work_terms.insert(4, work_forcing_r)
    integrated_terms.insert(4, integrated_forcing_r)
    titles.insert(4, r'$\mathbf{v}\cdot\mathbf{f}$')
    simple_labels.insert(4, 'forcing')

if magnetism:
    work_terms.insert(4, work_jcrossb)
    shav_work_terms.insert(4, work_jcrossb_r)
    integrated_terms.insert(4, integrated_jcrossb)
    titles.insert(4, r'$\frac{1}{4\pi}\langle\mathbf{u}\cdot[(\nabla\times\mathbf{B})\times\mathbf{B}]\rangle$')
    simple_labels.insert(4, 'u dot jcrossb')

# Set up figure from scratch
fig_width_inches = 11.5 # sideways paper
margin_inches = 1./8. # margin width in inches (for both x and y) and 
    # horizontally in between figures
margin_top_inches = 1 # wider top margin to accommodate subplot titles AND metadata
margin_bottom_inches = 0.75*(2 - (rbcz is None)) 
    # larger bottom margin to make room for colorbar(s)
margin_subplot_top_inches = 1/4 # margin to accommodate just subplot titles
nplots = 5 + magnetism + forced
nrow = 2
#ncol = np.int(np.ceil((nplots + 3)/nrow)) # +3 for room for line plot
#ncol = 4 # put three plots per row
#nrow = np.int(np.ceil(nplots/3)) + 1 # one more row for spherical avg plot

#subplot_width_inches = (fig_width_inches - (ncol + 1)*margin_inches)/ncol
subplot_width_inches = (fig_width_inches - (nplots + 1)*margin_inches)/nplots
    # Make the subplot width so that ncol subplots fit together side-by-side
    # with margins in between them and at the left and right.
subplot_height_inches = 2*subplot_width_inches # Each subplot should have an
    # aspect ratio of y/x = 2/1 to accommodate meridional planes. 
fig_height_inches = margin_top_inches + nrow*(subplot_height_inches +\
        margin_subplot_top_inches + margin_bottom_inches)
fig_aspect = fig_height_inches/fig_width_inches

# "Margin" in "figure units"; figure units extend from 0 to 1 in BOTH 
# directions, so unitless dimensions of margin will be different in x and y
# to force an equal physical margin
margin_x = margin_inches/fig_width_inches
margin_y = margin_inches/fig_height_inches
margin_top = margin_top_inches/fig_height_inches
margin_bottom = margin_bottom_inches/fig_height_inches
margin_subplot_top = margin_subplot_top_inches/fig_height_inches

# Subplot dimensions in figure units
subplot_width = subplot_width_inches/fig_width_inches
subplot_height = subplot_height_inches/fig_height_inches

# Generate figure of the correct dimensions
fig = plt.figure(figsize=(fig_width_inches, fig_height_inches))

# axes to hold spherically averaged line plot
#n_on_second_row = nplots - ncol
#ax_r_width = 2.*subplot_width # space for left label + make pretty
ax_r_width= 3*subplot_width
ax_r_left = 0.15
#ax_r_center = 0.5*(1. + margin_x +\
#        n_on_second_row*(subplot_width + margin_x))
#ax_r_left = ax_r_center - 0.5*ax_r_width
ax_r = fig.add_axes((ax_r_left, margin_bottom,\
        ax_r_width, subplot_height))
# linewidth for line plot
lw = 1.
# stuff for labels showing total rates of change
#fs = 6
#offset = 1.0/4.0/fig_height_inches

for iplot in range(nplots):
    #ax_left = margin_x + (iplot%ncol)*(subplot_width + margin_x)
    ax_left = margin_x + iplot*(subplot_width + margin_x)
    #ax_bottom = 1 - margin_top - subplot_height - margin_subplot_top -\
            #margin_bottom)
    ax_bottom = 1. - margin_top - subplot_height - margin_subplot_top 
    # plot work in meridional plane
    ax = fig.add_axes((ax_left, ax_bottom, subplot_width, subplot_height))
    plot_azav (work_terms[iplot], rr, cost, fig=fig, ax=ax, units=units,\
           minmax=minmax, plotcontours=plotcontours, rvals=rvals,\
           minmaxrz=minmaxrz, rbcz=rbcz, symlog=symlog,\
    linthresh=linthresh, linscale=linscale, linthreshrz=linthreshrz,\
    linscalerz=linscalerz, plotlatlines=plotlatlines)
    ax.set_title(titles[iplot], verticalalignment='bottom', **csfont)

    # plot spherically averaged work
    ax_r.plot(rr/rsun, shav_work_terms[iplot],\
            label=simple_labels[iplot] + ': ' + sci_format(integrated_terms[iplot]/lstar, 3) + r'$L_*$', linewidth=lw)

    # print total rates of change for the works
    #fig.text(0.75, margin_bottom + iplot*offset, simple_labels[iplot] +\
    #        ' total dE/dt = %1.3e erg/s' %integrated_terms[iplot],\
    #        fontsize=fs)

# fix up some stuff for the shav work terms
ax_r.set_xlabel(r'$r/R_\odot$')
ax_r.set_ylabel("spherically avg'd terms")
ax_r.set_xlim((ri/rsun, ro/rsun))
# make room for legend
ymin, ymax = ax_r.get_ylim()
buff_frac = 0.3
buff = buff_frac*(ymax - ymin)
ymin -= buff
ax_r.set_ylim((ymin, ymax))
leg = ax_r.legend(loc = 'lower left', fontsize=7, ncol=1)
# ticks everywhere
plt.sca(ax_r)
plt.minorticks_on()
plt.tick_params(top=True, right=True, direction='in', which='both')

# Label averaging interval
if rotation:
    time_string = ('t = %.1f to %.1f ' %(t1/time_unit, t2/time_unit))\
            + time_label + (r'$\ (\Delta t = %.1f\ $'\
            %((t2 - t1)/time_unit)) + time_label + ')'
else:
    time_string = ('t = %.3f to %.3f ' %(t1/time_unit, t2/time_unit))\
            + time_label + (r'$\ (\Delta t = %.3f\ $'\
            %((t2 - t1)/time_unit)) + time_label + ')'

# Put some metadata in upper left
fsize = 12
fig.text(margin_x, 1 - 0.1*margin_top, dirname_stripped,\
         ha='left', va='top', fontsize=fsize, **csfont)
fig.text(margin_x, 1 - 0.3*margin_top, 'Kinetic energy production terms (zonally averaged)',\
         ha='left', va='top', fontsize=fsize, **csfont)
fig.text(margin_x, 1 - 0.5*margin_top, time_string,\
         ha='left', va='top', fontsize=fsize, **csfont)

savefile = plotdir + dirname_stripped + '_keq_' + str(iter1).zfill(8) +\
    '_' + str(iter2).zfill(8) + tag + '.png'

if saveplot:
    print ('Saving plot at ' + savefile)
    plt.savefig(savefile, dpi=300)
if showplot:
    plt.show()
plt.close()
