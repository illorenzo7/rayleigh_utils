# Author: Loren Matilsky
# Created: 12/19/2022
#
# Description: Script to plot the mean shear (dOm/dr and 1/r DOm/dth)
# in meridional plane

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
from cla_util import *

# Get directory name and stripped_dirname for plotting purposes
args = sys.argv
clas = read_clas(args)
dirname = clas['dirname']
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
the_file = get_widest_range_file(datadir, 'AZ_Avgs')
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
        the_file = args[i+1]
        the_file = the_file.split('/')[-1]
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

# Get the terms:
print ('Getting terms from ' + the_file)
di = get_dict(the_file)

iter1, iter2 = get_iters_from_file(the_file)
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
di_grid = get_grid_info(dirname)
rr = di_grid['rr']
rr_2d = di_grid['rr_2d']
sint_2d = di_grid['sint_2d']
cost_2d = di_grid['cost_2d']
nr = di_grid['nr']
cost = di_grid['cost']
tt_lat = di_grid['tt_lat']
tt = di_grid['tt']
xx = di_grid['xx']

# Get <v_phi>
vp = vals[:, :, lut[3]]

# get angular momentum viscous fluxes, and thus derivatives of <v_phi>
amom_visc_r = vals[:, :, lut[1813]]
amom_visc_t = vals[:, :, lut[1814]]

eq = get_eq(dirname)
nu = eq.nu.reshape((1, nr))
rho = eq.rho.reshape((1, nr))
mu = rho*nu
prefactor = -1./(rho*nu*rr_2d**2.*sint_2d**2.)

# get diffrot and its derivs
Om = vp/(rr_2d*sint_2d)
dOmdr = prefactor*amom_visc_r
dOmdt = prefactor*amom_visc_t

# get <v_phi> derivs
dvpdr = 1./rr_2d*(vp - amom_visc_r/mu/sint_2d)
dvpdt = 1./rr_2d/sint_2d*(cost_2d*vp - amom_visc_t/mu)

terms = [Om, dOmdr, dOmdt, np.sqrt(dOmdr**2. + dOmdt**2.),\
        vp, dvpdr, dvpdt, np.sqrt(dvpdr**2. + dvpdt**2.)]

# field units and labels
units = [r'$\rm{s^{-1}}$', r'$\rm{s^{-1}\ cm^{-1}}$', r'$\rm{s^{-1}\ cm^{-1}}$', r'$\rm{s^{-1}\ cm^{-1}}$', r'$\rm{cm\ s^{-1}}$', r'$\rm{s^{-1}}$', r'$\rm{s^{-1}}$', r'$\rm{s^{-1}}$']
titles = ['<Om>', 'd<Om>/dr', 'd<Om>/dT', 'magnitude', '<v_phi>', 'd<v_phi>/dr', 'd<v_phi>/dT', 'magnitude']

# Set up the actual figure from scratch
fig_width_inches = 7. # TOTAL figure width, in inches
    # (i.e., 8x11.5 paper with 1/2-inch margins)
margin_inches = 1./8. # margin width in inches (for both x and y) and 
    # horizontally in between figures
margin_bottom_inches = 0.75*(2 - (rbcz is None)) 
    # larger bottom margin to make room for colorbar(s)
margin_top_inches = 1 # wider top margin to accommodate subplot titles AND metadata
margin_subplot_top_inches = 1/4 # margin to accommodate just subplot titles
nplots = len(terms)
ncol = 4 # put three plots per row
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

# Generate the actual figure of the correct dimensions
fig = plt.figure(figsize=(fig_width_inches, fig_height_inches))

for iplot in range(nplots):
    ax_left = margin_x + (iplot%ncol)*(subplot_width + margin_x)
    ax_bottom = 1 - margin_top - subplot_height - margin_subplot_top -\
            (iplot//ncol)*(subplot_height + margin_subplot_top +\
            margin_bottom)
    ax = fig.add_axes((ax_left, ax_bottom, subplot_width, subplot_height))
    plot_azav (terms[iplot], rr, cost, fig=fig, ax=ax, units=units[iplot],\
           minmax=minmax, plotcontours=plotcontours, rvals=rvals,\
           minmaxrz=minmaxrz, rbcz=rbcz, symlog=symlog,\
    linthresh=linthresh, linscale=linscale, linthreshrz=linthreshrz,\
    linscalerz=linscalerz, plotlatlines=plotlatlines,\
    plotboundary=plotboundary)

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
fig.text(margin_x, 1 - 0.3*margin_top, 'Diff. Rot. (derivs)',\
         ha='left', va='top', fontsize=fsize, **csfont)
fig.text(margin_x, 1 - 0.5*margin_top, time_string,\
         ha='left', va='top', fontsize=fsize, **csfont)

# save the figure
plotdir = my_mkdir(clas['plotdir'])
savefile = plotdir + 'diffrot_derivs' + clas['tag'] + '-' + str(iter1).zfill(8) + '_' + str(iter2).zfill(8) + '.png'

if clas['saveplot']:
    print ('saving figure at ' + savefile)
    plt.savefig(savefile, dpi=300)
if clas['showplot']:
    plt.show()
plt.close()
