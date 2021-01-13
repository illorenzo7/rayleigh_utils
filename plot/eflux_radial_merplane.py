# Author: Loren Matilsky
# Created: 01/29/2019
# This script plots the radial energy fluxes in the meridional plane 
# (convective (enthalpy) -- BOTH mean and fluc, conductive,
# radiative heating, kinetic energy,
# viscous, and Poynting flux (if present) 
# Computes the mean enthalpy flux "manually" from <v_r> and <T>
# In the future, just add quantity codes 1458,1459,1460!
# ...for the Rayleigh run directory indicated by [dirname]. 
# To use an AZ_Avgs file
# different than the one assosciated with the longest averaging range, use
# -usefile [complete name of desired AZ_Avgs file]
# Saves plot in
# [dirname]_eflux_radial_merplane_[first iter]_[last iter].png

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
from compute_grid_info import compute_theta_grid

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
plotdir = dirname + '/plots/'
if not os.path.isdir(plotdir):
    os.makedirs(plotdir)

# Read command-line arguments (CLAs)
showplot = True
saveplot = True
plotcontours = True
plotlatlines = True
plotboundary = True
minmax = None
linthresh = None
linscale = None
minmaxrz = None
linthreshrz = None
linscalerz = None
the_file = get_widest_range_file(datadir, 'AZ_Avgs')
forced = False
rvals = []
rbcz = None
symlog = False

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
    elif arg == '-usefile':
        the_file = args[i+1]
        the_file = the_file.split('/')[-1]
    elif arg == '-forced':
        forced = True
    elif arg == '-depths':
        strings = args[i+1].split()
        for st in strings:
            rval = ro - float(st)*d
            rvals.append(rval)
    elif arg == '-depthscz':
        rm = domain_bounds[1]
        dcz = ro - rm
        strings = args[i+1].split()
        for st in strings:
            rval = ro - float(st)*dcz
            rvals.append(rval)
    elif arg == '-depthsrz':
        rm = domain_bounds[1]
        drz = rm - ri
        strings = args[i+1].split()
        for st in strings:
            rval = rm - float(st)*drz
            rvals.append(rval)
    elif arg == '-rvals':
        rvals = []
        strings = args[i+1].split()
        for st in strings:
            rval = float(st)*rsun
            rvals.append(rval)
    elif arg == '-rvalscm':
        rvals = []
        strings = args[i+1].split()
        for st in strings:
            rval = float(st)
            rvals.append(rval)
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
    elif arg == '-nocontour':
        plotcontours = False
    elif arg == '-nobound':
        plotboundary = False
    elif arg == '-nolat':
        plotlatlines = False

# See if magnetism is "on"
magnetism = get_parameter(dirname, 'magnetism')

# Get AZ_Avgs file
print ('Getting radial energy fluxes (zonally averaged) from ' +\
        datadir + the_file)
di = get_dict(datadir + the_file)

# Get Shell_Avgs file
Shell_Avgs_file = get_widest_range_file(datadir, 'Shell_Avgs')
print ('Getting spherical averages from ' +\
        datadir + Shell_Avgs_file)
di_sph = get_dict(datadir + Shell_Avgs_file)

# Zonal average data
vals = di['vals']
lut = di['lut']
nq = di['nq']
iter1, iter2 = di['iter1'], di['iter2']

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

# Spherical average data
vals_sph = di_sph['vals']
lut_sph = di_sph['lut']
nq_sph = di_sph['nq']

# Get necessary grid info
rr = di['rr']
cost = di['cost']
sint = di['sint']
tt = di['tt']
tt_lat = di['tt_lat']
nr, nt = len(rr), len(tt)

# Get zonally averaged fluxes
efr_enth = vals[:, :, lut[1455]]
efr_cond = vals[:, :, lut[1470]]
efr_visc = -vals[:, :, lut[1935]]
efr_ke = vals[:, :, lut[1923]]

# get spherically averaged fluxes
efr_enth_sph = vals_sph[:, lut[1455]]
efr_cond_sph = vals_sph[:, lut[1470]]
efr_visc_sph = -vals_sph[:, lut[1935]]
efr_ke_sph = vals_sph[:, lut[1923]]


# Check to see if enthalpy flux from the fluctuating flows 
# was already output
if lut[1458] < nq and lut_sph[1458] < nq_sph:
    efr_enth_fluc = vals[:, :, lut[1458]]
    efr_enth_mean = efr_enth - efr_enth_fluc
    efr_enth_fluc_sph = vals_sph[:, lut[1458]]
    efr_enth_mean_sph = efr_enth_sph - efr_enth_fluc_sph
else: # do the Reynolds decomposition "by hand"
    # Compute the enthalpy flux from mean flows (MER. CIRC.)
    eq = get_eq(dirname)
    rho = (eq.density).reshape((1, nr))
    ref_prs = (eq.pressure).reshape((1, nr))
    ref_temp = (eq.temperature).reshape((1, nr))
    prs_spec_heat = 3.5e8

    gamma = 5./3.
    vr_av = vals[:, :, lut[1]]
    entropy_av = vals[:, :, lut[501]]
    prs_av = vals[:, :, lut[502]]

    # Calculate mean temp. from EOS
    temp_av = ref_temp*((1.-1./gamma)*(prs_av/ref_prs) +\
            entropy_av/prs_spec_heat)

    # And, finally, the enthalpy flux from mean/fluc flows
    efr_enth_mean = rho*prs_spec_heat*vr_av*temp_av
    efr_enth_fluc = efr_enth - efr_enth_mean

    # Finally finally, the spherically averaged parts
    tt, tw = compute_theta_grid(nt)
    tw_2d = tw.reshape((nt, 1))
    efr_enth_mean_sph = np.sum(efr_enth_mean*tw_2d, axis=0)
    efr_enth_fluc_sph = np.sum(efr_enth_fluc*tw_2d, axis=0)

# Subtract off the spherically averaged parts
efr_enth = efr_enth - efr_enth_sph
efr_enth_fluc = efr_enth_fluc - efr_enth_fluc_sph
efr_enth_mean = efr_enth_mean - efr_enth_mean_sph
efr_cond = efr_cond - efr_cond_sph
efr_ke = efr_ke - efr_ke_sph
efr_visc = efr_visc - efr_visc_sph

efr_tot = efr_enth + efr_cond + efr_visc + efr_ke

max_sig = max(np.std(efr_enth), np.std(efr_cond), np.std(efr_visc),\
        np.std(efr_ke))

if magnetism:
    efr_Poynt = vals[:, :, lut[2001]]       
    efr_Poynt_sph = vals_sph[:, lut[2001]]
    efr_Poynt = efr_Poynt - efr_Poynt_sph
    max_sig = max(max_sig, np.std(efr_Poynt))
    efr_tot += efr_Poynt
max_sig = max(max_sig, np.std(efr_tot)) 

if minmax is 'same': 
    minmax = -3.*max_sig, 3.*max_sig

# Set up the actual figure from scratch
fig_width_inches = 7. # TOTAL figure width, in inches
    # (i.e., 8x11.5 paper with 1/2-inch margins)
margin_inches = 1./8. # margin width in inches (for both x and y) and 
    # horizontally in between figures
margin_top_inches = 1 # wider top margin to accommodate subplot titles AND metadata
margin_bottom_inches = 0.75*(2 - (rbcz is None)) 
    # larger bottom margin to make room for colorbar(s)
margin_subplot_top_inches = 1/4 # margin to accommodate just subplot titles
nplots = 7 + magnetism
ncol = 3 # put three plots per row
nrow = np.int(np.ceil(nplots/3))

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

efr_terms = [efr_enth_mean, efr_enth_fluc, efr_enth, efr_cond, efr_visc,\
        efr_ke, efr_tot]

titles = [r'$(\mathbf{\mathcal{F}}_{\rm{enth,mm}})_r$',\
          r'$(\mathbf{\mathcal{F}}_{\rm{enth,pp}})_r$',\
          r'$(\mathbf{\mathcal{F}}_{\rm{enth}})_r$',\
          r'$(\mathbf{\mathcal{F}}_{\rm{cond}})_r$',\
          r'$(\mathbf{\mathcal{F}}_{\rm{visc}})_r$',\
          r'$(\mathbf{\mathcal{F}}_{\rm{KE}})_r$',
          r'$(\mathbf{\mathcal{F}}_{\rm{tot}})_r$']
units = r'$\rm{g}\ \rm{s}^{-3}$'

if magnetism:
    # Insert the magnetism plot to just before the last
    # plot (total flux)
    efr_terms.insert(7, efr_Poynt)
    titles.insert(7, r'$(\mathbf{\mathcal{F}}_{\rm{Poynt}})_r$')

# Generate the actual figure of the correct dimensions
fig = plt.figure(figsize=(fig_width_inches, fig_height_inches))

for iplot in range(nplots):
    ax_left = margin_x + (iplot%ncol)*(subplot_width + margin_x)
    ax_bottom = 1 - margin_top - subplot_height - margin_subplot_top -\
            (iplot//ncol)*(subplot_height + margin_subplot_top +\
            margin_bottom)
    ax = fig.add_axes((ax_left, ax_bottom, subplot_width, subplot_height))
    plot_azav (efr_terms[iplot] - np.mean(efr_terms[iplot], axis=0),\
            rr, cost, fig=fig, ax=ax, units=units,\
           minmax=minmax, plotcontours=plotcontours, rvals=rvals,\
           minmaxrz=minmaxrz, rbcz=rbcz, symlog=symlog,\
    linthresh=linthresh, linscale=linscale, linthreshrz=linthreshrz,\
    linscalerz=linscalerz, plotlatlines=plotlatlines, plotboundary=plotboundary)
    ax.set_title(titles[iplot], verticalalignment='bottom', **csfont)

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
fig.text(margin_x, 1 - 0.3*margin_top, 'Radial energy flux, ' +\
        r'$l\neq0$', ha='left', va='top', fontsize=fsize, **csfont)
fig.text(margin_x, 1 - 0.5*margin_top, time_string,\
         ha='left', va='top', fontsize=fsize, **csfont)

savefile = plotdir + dirname_stripped + '_eflux_radial_merplane_' +\
    str(iter1).zfill(8) + '_' + str(iter2).zfill(8) + '.png'

if saveplot:
    print ('Saving radial energy fluxes (in the meridional plane) at ' +\
       savefile)
    plt.savefig(savefile, dpi=300)
if showplot:
    plt.show()
plt.close()
