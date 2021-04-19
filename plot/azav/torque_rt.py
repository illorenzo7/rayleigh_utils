# Author: Loren Matilsky
# Created: 03/07/2021
# This script plots the axial torques in the meridional plane (viscous, 
# Meridional Circ., Reynolds stress, and Maxwell torques (mean and turbulent) 
# Broken into contributions from radial/latitudinal amom flux separately
# if applicablein the meridional plane 
# ...for the Rayleigh run directory indicated by [dirname]. To use an AZ_Avgs file
# different than the one associated with the longest averaging range, use
# -usefile [complete name of desired AZ_Avgs file]
# Saves plot in
# [dirname]_torque_[first iter]_[last iter].png

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
nadd = 0
nsubset = []
torques_to_add = []
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
AZ_Avgs_file = get_widest_range_file(datadir, 'AZ_Avgs')
forced = False
rvals = []
rbcz = None
symlog = False
tag = ''

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
    elif arg == '-noshow':
        showplot = False
    elif arg == '-nosave':
        saveplot = False
    elif arg == '-usefile':
        AZ_Avgs_file = args[i+1]
        AZ_Avgs_file = AZ_Avgs_file.split('/')[-1]
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
    elif arg == '-add':
        loc_list = args[i+1].split()
        nadd += 1
        for j in range(len(loc_list)):
            loc_list[j] = int(loc_list[j])
        torques_to_add += loc_list
        nsubset.append(len(loc_list))
    elif arg == '-tag':
        tag = '_' + args[i+1]

# See if magnetism is "on"
try:
    magnetism = get_parameter(dirname, 'magnetism')
except:
    magnetism = False # if magnetism wasn't specified, it must be "off"

# Get the torques:
print ('Getting amom fluxes from ' + datadir + AZ_Avgs_file)
di = get_dict(datadir + AZ_Avgs_file)

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
rr = di['rr']
rr_2d = di['rr_2d']
cost = di['cost']
cost_2d = di['cost_2d']
tt = np.arccos(cost)
sint = di['sint']
sint_2d = di['sint_2d']
cott = cost_2d/sint_2d
tt_lat = di['tt_lat']
xx = di['xx']
nr, nt = di['nr'], di['nt']

torque_rs = -vals[:, :, lut[1801]]
torque_mc = -vals[:, :, lut[1802]] + vals[:, :, lut[1803]]
torque_visc = vals[:, :, lut[1804]]
if magnetism:
    torque_mm = vals[:, :, lut[1805]]
    torque_ms = vals[:, :, lut[1806]]
 
f_rs_r = vals[:, :, lut[1807]]
f_rs_t = vals[:, :, lut[1808]] 

f_mc_r = vals[:, :, lut[1809]] + vals[:, :, lut[1811]]
f_mc_t = vals[:, :, lut[1810]] + vals[:, :, lut[1812]]

f_v_r = vals[:, :, lut[1813]]
f_v_t = vals[:, :, lut[1814]]

torque_rs_r = -(drad(f_rs_r, rr) + 2.*f_rs_r/rr)
torque_rs_t = -( (1./rr)*(dth(f_rs_t, tt) + cott*f_rs_t) )

torque_mc_r = -(drad(f_mc_r, rr) + 2.*f_mc_r/rr)
torque_mc_t = -( (1./rr)*(dth(f_mc_t, tt) + cott*f_mc_t) )

torque_visc_r = -(drad(f_v_r, rr) + 2.*f_v_r/rr)
torque_visc_t = -( (1./rr)*(dth(f_v_t, tt) + cott*f_v_t) )

if magnetism:
    f_ms_r = vals[:, :, lut[1815]]
    f_ms_t = vals[:, :, lut[1816]]
    f_mm_r = vals[:, :, lut[1817]]
    f_mm_t = vals[:, :, lut[1818]]
    torque_ms_r = -(drad(f_ms_r, rr) + 2.*f_ms_r/rr)
    torque_ms_t = -( (1./rr)*(dth(f_ms_t, tt) + cott*f_ms_t) )
    torque_mm_r = -(drad(f_mm_r, rr) + 2.*f_mm_r/rr)
    torque_mm_t = -( (1./rr)*(dth(f_mm_t, tt) + cott*f_mm_t) )

# Set up the actual figure from scratch
fig_width_inches = 7. # TOTAL figure width, in inches
    # (i.e., 8x11.5 paper with 1/2-inch margins)
margin_inches = 1./8. # margin width in inches (for both x and y) and 
    # horizontally in between figures
margin_bottom_inches = 0.75*(2 - (rbcz is None)) 
    # larger bottom margin to make room for colorbar(s)
margin_top_inches = 1 # wider top margin to accommodate subplot titles AND metadata
margin_subplot_top_inches = 1/4 # margin to accommodate just subplot titles
nplots = 5
ncol = 3 # put three plots per row
nrow = np.int(np.ceil(nplots/ncol))

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

terms = {'rs': [torque_rs_r, torque_rs_t, torque_rs],\
        'mc': [torque_mc_r, torque_mc_t, torque_mc],\
        'visc': [torque_visc_r, torque_visc_t, torque_visc]}
base_titles = {'rs': r'$\tau_{\rm{rs}}$', 'mc':  r'$\tau_{\rm{mc}}$', 'visc': r'$\tau_{\rm{v}}$'}

if magnetism:
    terms['ms'] = [torque_ms_r, torque_ms_t, torque_ms]
    terms['mm'] = [torque_mm_r, torque_mm_t, torque_mm]
    base_titles['ms'] = r'$\tau_{\rm{ms}}$'
    base_titles['mm'] = r'$\tau_{\rm{mm}}$'

units = r'$\rm{g}\ \rm{cm}^{-1}\ \rm{s}^{-2}$'

for key in terms.keys():
    torque_r, torque_t, torque_exact = terms[key]
    torque_sum = torque_r + torque_t
    diff = torque_sum - torque_exact
    torques = [torque_r, torque_t, torque_sum, torque_exact, diff]

    base_title = base_titles[key]
    titles = [base_title + r'$_r$', base_title + r'$_\theta$', 'sum',\
            base_title, 'sum - ' + base_title]

    # Generate the actual figure of the correct dimensions
    fig = plt.figure(figsize=(fig_width_inches, fig_height_inches))

    for iplot in range(nplots):
        ax_left = margin_x + (iplot%ncol)*(subplot_width + margin_x)
        ax_bottom = 1 - margin_top - subplot_height - margin_subplot_top -\
                (iplot//ncol)*(subplot_height + margin_subplot_top +\
                margin_bottom)
        ax = fig.add_axes((ax_left, ax_bottom, subplot_width, subplot_height))
        plot_azav (torques[iplot], rr, cost, fig=fig, ax=ax, units=units,\
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
    fig.text(margin_x, 1 - 0.3*margin_top, 'Torque balance (radial and latitudinal separate)',\
             ha='left', va='top', fontsize=fsize, **csfont)
    fig.text(margin_x, 1 - 0.5*margin_top, time_string,\
             ha='left', va='top', fontsize=fsize, **csfont)

    savefile = plotdir + dirname_stripped + '_torque_rt_' + str(iter1).zfill(8) +\
        '_' + str(iter2).zfill(8) + tag + '_' + key + '.png'

    if saveplot:
        print ('Saving torques at ' + savefile)
        plt.savefig(savefile, dpi=300)
    plt.close()
