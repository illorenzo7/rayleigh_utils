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
from common import get_widest_range_file, strip_dirname, get_dict, rsun, c_P
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
minmax2 = None
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
novisc = True # don't include viscosity in total by default 
crudeint = False # by default use exact latitudinal weights and instead 
# if crudeint = True, use a "crude weight"
tag = ''

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
    elif arg == '-visc':
        novisc = False
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
work_KE_advec = vals[:, :, lut[1910]] #

# Pressure work on KE
work_pressure = vals[:, :, lut[1901]]

# Buoynacy work on KE (in future calculate 1904)
rhoTvrS = vals[:, :, lut[1440]]
vrS = rhoTvrS/rhoT
S_00 = np.sum(vals[:, :, lut[501]]*tw_2d, axis=0)
vr = vals[:, :, lut[1]]
vrS_00 = vr*S_00.reshape((1, nr))
work_buoy = (rho*grav/c_P).rshape((1, nr))*(vrS - vrS_00)

# Viscous work on kinetic energy
work_visc_on_KE = vals[:, :, lut[1907]]

# J X B work on KE, if available
if magnetism:
    work_jcrossb = vals[:, :, lut[1916]]

#========================
# HEAT (ENTROPY) EQUATION
#========================

# Advection of entropy
work_thermal_advec = -vals[:, :, lut[1401]] # 1401 = + rho*T*v dot grad S

# Stable gradient (advection of reference entropy)
work_thermal_advec_ref = -rhoT_2d*vr*dsdr_2d

# heat by conduction
work_cond = vals[:, :, lut[1421]]

# heat from heating function (representing radiation)
work_rad = vals[:, :, lut[1434]] # Q(r)

# Irreversible (viscous) heating
work_visc_on_intE = rhoT_2d*vals[:, :, lut[1435]] 
# (irreversible) viscous heating
# THIS IS OFF BY rho*T compared to what it's supposed to be (i.e., it's 
# the viscous heating as in the energy equation

# If magnetism get Joule heating

work_enth = work_thermal_advec + work_pressure

print ("std (visc work on KE): ", np.std(work_visc_on_KE))
print ("std (visc_heating: visc work on intE): ", np.std(work_visc_on_intE))
work_visc = work_visc_on_KE + work_visc_on_intE


work_dsdr_negligible = rhoTvrS*dsdr_2d/c_P
print ("maxabs (small dsdr work) = ", np.max(np.abs(work_dsdr_negligible)))
if magnetism:
    work_mag = vals[:, :, lut[1436]]  +\
            1.0/(4.0*np.pi)*(vals[:, :, lut[2019]] + vals[:, :, lut[2043]])
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

    work_forcing_shav = np.sum(work_forcing*tw.reshape((nt, 1)), axis=0)

# Calculate the integrated work
fourpi = 4.0*np.pi

work_KE_r = np.sum(work_KE*tw_2d, axis=0)
work_enth_r = np.sum(work_enth*tw_2d, axis=0)
work_cond_r = np.sum(work_cond*tw_2d, axis=0)
work_rad_r = np.sum(work_rad*tw_2d, axis=0)
work_visc_r = np.sum(work_visc*tw_2d, axis=0)
work_dsdr_r = np.sum(work_dsdr*tw_2d, axis=0)
work_dsdr_negligible_r = np.sum(work_dsdr_negligible*tw_2d, axis=0)
if magnetism:
    work_mag_r = np.sum(work_mag*tw_2d, axis=0)
if forced:
    work_forcing_r = np.sum(work_forcing*tw_2d, axis=0)

work_tot = work_KE + work_enth + work_cond + work_rad + work_dsdr +\
        work_dsdr_negligible
work_tot_r = work_KE_r + work_enth_r + work_cond_r + work_rad_r +\
        work_dsdr_r + work_dsdr_negligible_r
if novisc:
    print ("not including visc. work in total")
    print ("to include it, specify -visc")
else:
    work_tot += work_visc
    work_tot_r += work_visc_r

if magnetism:
    work_tot += work_mag
    work_tot_r += work_mag_r
if forced:
    if not force_econs:
        print("force_econs = False so including forcing work in total")
        work_tot += work_forcing
        work_tot_r += work_forcing_r
    else:
        print("force_econs = True so not including forcing work in total")

ri, ro = np.min(rr), np.max(rr)
integrated_KE = fourpi/3*(ro**3 - ri**3)*np.sum(work_KE_r*rw)
integrated_enth = fourpi/3*(ro**3 - ri**3)*np.sum(work_enth_r*rw)
integrated_cond = fourpi/3*(ro**3 - ri**3)*np.sum(work_cond_r*rw)
integrated_rad = fourpi/3*(ro**3 - ri**3)*np.sum(work_rad_r*rw)
integrated_visc = fourpi/3*(ro**3 - ri**3)*np.sum(work_visc_r*rw)
integrated_dsdr = fourpi/3*(ro**3 - ri**3)*np.sum(work_dsdr_r*rw)
integrated_dsdr_negligible = fourpi/3*(ro**3 - ri**3)*np.sum(work_dsdr_negligible_r*rw)
integrated_total = fourpi/3*(ro**3 - ri**3)*np.sum(work_tot_r*rw)
if magnetism:
    integrated_mag = fourpi/3*(ro**3 - ri**3)*np.sum(work_mag_r*rw)
if forced:
    integrated_forcing = fourpi/3*(ro**3 - ri**3)*np.sum(work_forcing_r*rw)

#max_sig = max(np.std(torque_rs), np.std(torque_mc), np.std(torque_visc))

# Set up the actual figure from scratch
fig_width_inches = 7. # TOTAL figure width, in inches
    # (i.e., 8x11.5 paper with 1/2-inch margins)
margin_inches = 1./8. # margin width in inches (for both x and y) and 
    # horizontally in between figures
margin_top_inches = 1 # wider top margin to accommodate subplot titles AND metadata
margin_bottom_inches = 0.75*(2 - (rbcz is None)) 
    # larger bottom margin to make room for colorbar(s)
margin_subplot_top_inches = 1/4 # margin to accommodate just subplot titles
nplots = 8 + magnetism + forced
ncol = 3 # put three plots per row
nrow = np.int(np.ceil(nplots/3)) + 1 # one more row for spherical avg plot

subplot_width_inches = (fig_width_inches - (ncol + 1)*margin_inches)/ncol
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

work_terms = [work_KE, work_enth, work_cond, work_rad, work_visc,\
        work_dsdr, work_dsdr_negligible, work_tot]
titles = [r'$-\nabla\cdot\mathbf{F}_{\rm{KE}}$',\
r'$-\nabla\cdot\mathbf{F}_{\rm{enth}}$',\
r'$-\nabla\cdot\mathbf{F}_{\rm{cond}}$',\
r'$-\nabla\cdot\mathbf{F}_{\rm{rad}}$',\
r'$-\nabla\cdot\mathbf{F}_{\rm{visc}}$',\
r'$-\overline{\rho}\overline{T}\mathbf{v}\cdot\nabla\overline{S}$',\
r'$\overline{\rho}\overline{T}(S/c_{\rm{p}})\mathbf{v}\cdot\nabla\overline{S}$',\
r'$\partial(\overline{\rho}w)/\partial t$']

units = r'$\rm{erg}\ \rm{cm}^{-3}\ \rm{s}^{-1}$'

if forced:
    work_terms.insert(7, work_forcing)
    titles.insert(7, r'$\mathbf{v}\cdot\mathbf{f}$')

if magnetism:
    work_terms.insert(7, work_mag)
    titles.insert(7, r'$-\nabla\cdot\mathbf{F}_{\rm{Poynting}}$')

# Generate the actual figure of the correct dimensions
fig = plt.figure(figsize=(fig_width_inches, fig_height_inches))

for iplot in range(nplots):
    ax_left = margin_x + (iplot%ncol)*(subplot_width + margin_x)
    ax_bottom = 1 - margin_top - subplot_height - margin_subplot_top -\
            (iplot//ncol)*(subplot_height + margin_subplot_top +\
            margin_bottom)
    ax = fig.add_axes((ax_left, ax_bottom, subplot_width, subplot_height))
    plot_azav (work_terms[iplot], rr, cost, fig=fig, ax=ax, units=units,\
           minmax=minmax, plotcontours=plotcontours, rvals=rvals,\
           minmaxrz=minmaxrz, rbcz=rbcz, symlog=symlog,\
    linthresh=linthresh, linscale=linscale, linthreshrz=linthreshrz,\
    linscalerz=linscalerz, plotlatlines=plotlatlines)
    ax.set_title(titles[iplot], verticalalignment='bottom', **csfont)

# Plot the latitudinally averaged divergences
ax_r = fig.add_axes((margin_x + 0.75*subplot_width, margin_bottom,\
        1.5*subplot_width, subplot_height))
lw = 0.5
ax_r.plot(rr/rsun, work_KE_r, label='KE', linewidth=lw)
ax_r.plot(rr/rsun, work_enth_r, label='enth', linewidth=lw)
ax_r.plot(rr/rsun, work_cond_r, label='cond', linewidth=lw)
ax_r.plot(rr/rsun, work_rad_r, label='rad', linewidth=lw)
ax_r.plot(rr/rsun, work_visc_r, label='visc', linewidth=lw)
ax_r.plot(rr/rsun, work_dsdr_r, label='dsdr', linewidth=lw)
ax_r.plot(rr/rsun, work_dsdr_negligible_r, label='dsdr\n(small)', linewidth=lw)
if magnetism:
    ax_r.plot(rr/rsun, work_mag_r, label='mag', linewidth=lw)
if forced:
    ax_r.plot(rr/rsun, work_forcing_r, label='forcing', linewidth=lw)

ax_r.plot(rr/rsun, work_tot_r, label='tot', linewidth=lw)
ax_r.set_xlabel(r'$r/R_\odot$')
ax_r.set_ylabel("spherically avg'd terms")
ax_r.set_xlim((ri/rsun, ro/rsun))
if not minmax2 is None:
    ax_r.set_ylim((minmax2[0], minmax2[1]))

leg = ax_r.legend(loc = (-0.4, 0), fontsize=7)

# Output total rates of change of energy due to various terms
fs = 6
offset = 1.0/4.0/fig_height_inches
fig.text(0.75, margin_bottom + 0*offset, 'total dE/dt = %1.3e erg/s' %integrated_total, fontsize=fs)
fig.text(0.75, margin_bottom + 1*offset, 'dsdr (small) dE/dt = %1.3e erg/s' %integrated_dsdr_negligible, fontsize=fs)
fig.text(0.75, margin_bottom + 2*offset, 'dsdr dE/dt = %1.3e erg/s' %integrated_dsdr, fontsize=fs)
fig.text(0.75, margin_bottom + 3*offset, 'visc dE/dt = %1.3e erg/s' %integrated_visc, fontsize=fs)
fig.text(0.75, margin_bottom + 4*offset, 'rad dE/dt = %1.3e erg/s' %integrated_rad, fontsize=fs)
fig.text(0.75, margin_bottom + 5*offset, 'cond dE/dt = %1.3e erg/s' %integrated_cond, fontsize=fs)
fig.text(0.75, margin_bottom + 6*offset, 'enth dE/dt = %1.3e erg/s' %integrated_enth, fontsize=fs)
fig.text(0.75, margin_bottom + 7*offset, 'KE dE/dt = %1.3e erg/s' %integrated_KE, fontsize=fs)
current_offset = 7*offset
if magnetism:
    current_offset += offset
    fig.text(0.75, margin_bottom + current_offset, 'mag dE/dt = %1.3e erg/s' %integrated_mag, fontsize=fs)
if forced:
    current_offset += offset
    fig.text(0.75, margin_bottom + current_offset, 'forcing dE/dt = %1.3e erg/s' %integrated_forcing, fontsize=fs)

# Get ticks everywhere
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
fig.text(margin_x, 1 - 0.3*margin_top, 'Total E production terms (zonally averaged)',\
         ha='left', va='top', fontsize=fsize, **csfont)
fig.text(margin_x, 1 - 0.5*margin_top, time_string,\
         ha='left', va='top', fontsize=fsize, **csfont)

savefile = plotdir + dirname_stripped + '_toten_' + str(iter1).zfill(8) +\
    '_' + str(iter2).zfill(8) + tag + '.png'

if saveplot:
    print ('Saving plot at ' + savefile)
    plt.savefig(savefile, dpi=300)
if showplot:
    plt.show()
plt.close()
