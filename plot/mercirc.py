# Author: Loren Matilsky
# Created: 05/14/2018
# This script plots the meridional circulation cells the meridional plane 
# for the Rayleigh run directory indicated by [dirname], using the AZ_Avgs
# data. To use an AZ_Avgs file different than the one associated with the 
# longest averaging range, run with option
# -usefile [complete name of desired vavg file]
# Saves plot in
# [dirname]_diffrot_[first iter]_[last iter].npy

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
from azav_util import plot_azav, streamfunction
from common import *

# Get directory name and stripped_dirname for plotting purposes
dirname = sys.argv[1]
dirname_stripped = strip_dirname(dirname)

# domain bounds
ncheby, domain_bounds = get_domain_bounds(dirname)
ri = np.min(domain_bounds)
ro = np.max(domain_bounds)
d = ro - ri

# Split dirname_stripped into two lines if it is very long
if len(dirname_stripped) > 25:
    dirname_stripped_title = dirname_stripped[:25] + '\n' +\
            dirname_stripped[25:]
else:
    dirname_stripped_title = dirname_stripped

# Get density
eq = get_eq(dirname)
rho = eq.density

# Directory with data and plots, make the plotting directory if it doesn't
# already exist    
datadir = dirname + '/data/'
plotdir = dirname + '/plots/'
if not os.path.isdir(plotdir):
    os.makedirs(plotdir)

# Set defaults
minmax = None
linthresh = None
linscale = None
minmaxrz = None
linthreshrz = None
linscalerz = None
AZ_Avgs_file = get_widest_range_file(datadir, 'AZ_Avgs')
rbcz = None
symlog = False
plotcontours = True
plotlatlines = True
plotboundary = True
rvals = []

# Read in CLAs (if any) to change default variable ranges and other options
args = sys.argv[2:]
nargs = len(args)
for i in range(nargs):
    arg = args[i]
    if arg == '-minmax':
        minmax = float(args[i+1]), float(args[i+2])
    elif arg == '-minmaxrz':
        minmaxrz = float(args[i+1]), float(args[i+2])
    elif arg == '-nlevs':
        my_nlevs = int(args[i+1])
    elif arg == '-usefile':
        AZ_Avgs_file = args[i+1]
        AZ_Avgs_file = AZ_Avgs_file.split('/')[-1]
    elif arg == '-nocontour':
        plotcontours = False
    elif arg == '-nobound':
        plotboundary = False
    elif arg == '-nolat':
        plotlatlines = False
    elif arg == '-rbcz':
        rbcz = float(args[i+1])
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
    elif arg == '-depths':
        rvals = []
        strings = args[i+1].split()
        for st in strings:
            rval = ro - float(st)*d
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

# Read in AZ_Avgs data
print ('Getting data from ' + datadir + AZ_Avgs_file + ' ...')
di = get_dict(datadir + AZ_Avgs_file)

vals = di['vals']
lut = di['lut']
iter1, iter2 = di['iter1'], di['iter2']
rr = di['rr']
tt = di['tt']
cost = di['cost']
sint = di['sint']

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

vr_av, vt_av, vp_av = vals[:, :, lut[1]], vals[:, :, lut[2]],\
        vals[:, :, lut[3]]

# Compute the mass flux
rhovm = rho*np.sqrt(vr_av**2 + vt_av**2)

# Compute the streamfunction
psi = streamfunction(rho*vr_av, rho*vt_av, rr, cost)

# Make CCW negative and CW positive
rhovm *= np.sign(psi)

# Create plot
subplot_width_inches = 2.5
subplot_height_inches = 5.
margin_inches = 1./8.
margin_top_inches = 1.5 # larger top margin to make room for titles
margin_bottom_inches = 0.75*(2 - (rbcz is None)) 
    # larger bottom margin to make room for colorbar(s)

fig_width_inches = subplot_width_inches + 2*margin_inches
fig_height_inches = subplot_height_inches + margin_top_inches +\
        margin_bottom_inches

fig_aspect = fig_height_inches/fig_width_inches
margin_x = margin_inches/fig_width_inches
margin_y = margin_inches/fig_height_inches
margin_top = margin_top_inches/fig_height_inches
margin_bottom = margin_bottom_inches/fig_height_inches
subplot_width = subplot_width_inches/fig_width_inches
subplot_height = subplot_height_inches/fig_height_inches

fig = plt.figure(figsize=(fig_width_inches, fig_height_inches))
ax = fig.add_axes((margin_x, margin_bottom, subplot_width, subplot_height))

# Plot mass flux
plot_azav (rhovm, rr, cost, fig=fig, ax=ax,\
    units = r'$\rm{g}\ \rm{cm}^{-2}\ \rm{s}^{-1}$', plotcontours=False,\
    minmax=minmax, minmaxrz=minmaxrz, rbcz=rbcz, symlog=symlog,\
    linthresh=linthresh, linscale=linscale, linthreshrz=linthreshrz,\
    linscalerz=linscalerz, plotlatlines=plotlatlines, rvals=rvals,\
    plotboundary=plotboundary)

# Plot streamfunction contours, if desired
if plotcontours:
    lilbit = 0.01
    maxabs = np.max(np.abs(psi))
    levels = (-maxabs/2., -maxabs/4., -lilbit*maxabs, 0., lilbit*maxabs,\
            maxabs/4., maxabs/2.)
    plot_azav (psi, rr, cost, fig=fig, ax=ax, plotfield=False,\
        levels=levels, symlog=symlog, plotlatlines=plotlatlines,\
        plotboundary=plotboundary)

# Label averaging interval
if rotation:
    time_string = ('t = %.1f to %.1f ' %(t1/time_unit, t2/time_unit))\
            + time_label + '\n' + (r'$\ (\Delta t = %.1f\ $'\
            %((t2 - t1)/time_unit)) + time_label + ')'
else:
    time_string = ('t = %.3f to %.3f ' %(t1/time_unit, t2/time_unit))\
            + time_label + (r'$\ (\Delta t = %.3f\ $'\
            %((t2 - t1)/time_unit)) + time_label + ')'

# Make title
fsize = 12
fig.text(margin_x, 1 - 1/8*margin_top, dirname_stripped_title,\
         ha='left', va='top', fontsize=fsize, **csfont)
fig.text(margin_x, 1 - 3/8*margin_top,\
         r'$|\langle\overline{\rho}\mathbf{v}_m\rangle|$',\
         ha='left', va='top', fontsize=fsize, **csfont)
fig.text(margin_x, 1 - 5/8*margin_top,\
        time_string, ha='left', va='top', fontsize=fsize, **csfont)

savefile = plotdir + dirname_stripped + '_mercirc_' + str(iter1).zfill(8) +\
    '_' + str(iter2).zfill(8) + '.png'
print ('Saving plot at %s ...' %savefile)
plt.savefig(savefile, dpi=300)
plt.show()
