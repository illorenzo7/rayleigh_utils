# Author: Loren Matilsky
# Created: 04/06/2019
# This script plots mean magnetic field components in the meridional plane 
# for the Rayleigh run directory indicated by [dirname]. 
# Uses an instantaneous AZ_Avgs file unless told otherwise
# Saves plot in
# [dirname]_B_azav_[iter].npy
# or [dirname]_B_azav_[first iter]_[last iter].npy if a time average was 
# specified

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
from rayleigh_diagnostics import AZ_Avgs

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

the_file = get_widest_range_file(datadir, 'AZ_Avgs')
rvals = []
rbcz = None

args = sys.argv[2:]
nargs = len(args)
for i in range(nargs):
    arg = args[i]
    if arg == '-minmax':
        minmax = float(args[i+1]), float(args[i+2])
    elif arg == '-rbcz':
        rbcz = float(args[i+1])
    elif arg == '-noshow':
        showplot = False
    elif arg == '-nosave':
        saveplot = False
    elif arg == '-nocontour':
        plotcontours = False
    elif arg == '-nobound':
        plotboundary = False
    elif arg == '-nolat':
        plotlatlines = False
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
       
# Read in AZ_Avgs data
print ('Getting data from ' + datadir + the_file)
di = get_dict(datadir + the_file)

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

# Grid info
rr = di['rr']
cost = di['cost']
sint = di['sint']
ro = di['ro']

br, bt = vals[:, :, lut[801]], vals[:, :, lut[802]]

# Compute the magnitude of B_m
bm = np.sqrt(br**2 + bt**2)

# Compute the streamfunction
psi = streamfunction(br, bt, rr, cost)

# Make CCW negative and CW positive
bm *= np.sign(psi)

if minmax is None:
    std = np.std(bm)
    nstd = 3.
    minmax = -nstd*std, nstd*std

# Create plot
subplot_width_inches = 2.5
subplot_height_inches = 5.
margin_inches = 1/8
margin_top_inches = 1.25 # larger top margin to make room for titles
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

# Plot B_m magnitude
plot_azav (bm, rr, cost, fig=fig, ax=ax,\
    units = r'$\rm{G}$', plotcontours=False,\
    minmax=minmax, plotlatlines=plotlatlines)

# Plot streamfunction contours
#maxabs = np.max(np.abs(psi))
maxabs = 3*np.std(psi)
levels = np.linspace(-maxabs, maxabs, 20)
plot_azav (psi, rr, cost, fig=fig, ax=ax, plotfield=False, levels=levels)

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
line_height = 1./4./fig_height_inches
fig.text(margin_x, 1 - margin_y, dirname_stripped,\
         ha='left', va='top', fontsize=fsize, **csfont)
fig.text(margin_x, 1 - margin_y - line_height,\
         r'$|\langle\mathbf{B}_m\rangle|$',\
         ha='left', va='top', fontsize=fsize, **csfont)
fig.text(margin_x, 1 - margin_y - 2*line_height,\
        time_string, ha='left', va='top', fontsize=fsize, **csfont)

savefile = plotdir + dirname_stripped + '_Bm_' + str(iter1).zfill(8) +\
    '_' + str(iter2).zfill(8) + '.png'
print ('Saving plot at %s' %savefile)
if saveplot:
    plt.savefig(savefile, dpi=300)
if showplot:
    plt.show()
