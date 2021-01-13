# Author: Loren Matilsky
# Created: 04/06/2019
# This script plots mean velocity components in the meridional plane 
# for the Rayleigh run directory indicated by [dirname]. 
# Uses an instantaneous AZ_Avgs file unless told otherwise
# Saves plot in
# [dirname]_v_azav_[iter].npy
# or [dirname]_v_azav_[first iter]_[last iter].npy if a time average was 
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
from azav_util import plot_azav
from common import *
from rayleigh_diagnostics import AZ_Avgs

# Get directory name and stripped_dirname for plotting purposes
dirname = sys.argv[1]
dirname_stripped = strip_dirname(dirname)

# Directory with data and plots, make the plotting directory if it doesn't
# already exist    
datadir = dirname + '/data/'
plotdir = dirname + '/plots/'
if (not os.path.isdir(plotdir)):
    os.makedirs(plotdir)

radatadir = dirname + '/AZ_Avgs/'

# Read command-line arguments (CLAs)
showplot = True
saveplot = True
plotcontours = True
plotlatlines = True
plotboundary = True

minmaxvr = None
minmaxvt = None
minmaxvp = None

minmaxvrrz = None
minmaxvtrz = None
minmaxvprz = None

the_file = get_widest_range_file(datadir, 'AZ_Avgs')
rvals = []
rbcz = None

args = sys.argv[2:]
nargs = len(args)
for i in range(nargs):
    arg = args[i]
    if arg == '-minmaxvr':
        minmaxvr = float(args[i+1]), float(args[i+2])
    elif arg == '-minmaxvt':
        minmaxvt = float(args[i+1]), float(args[i+2])
    elif arg == '-minmaxvp':
        minmaxvp = float(args[i+1]), float(args[i+2])
    elif arg == '-minmaxvrrz':
        minmaxvrrz = float(args[i+1]), float(args[i+2])
    elif arg == '-minmaxvtrz':
        minmaxvtrz = float(args[i+1]), float(args[i+2])
    elif arg == '-minmaxvprz':
        minmaxvprz = float(args[i+1]), float(args[i+2])
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
    elif arg == '-rvals':
        rvals_str = args[i+1].split()
        rvals = []
        for rval_str in rvals_str:
            rvals.append(float(rval_str))
       
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

vr, vt, vp = vals[:, :, lut[1]]/100., vals[:, :, lut[2]]/100.,\
        vals[:, :, lut[3]]/100.

# Set up the actual figure from scratch
fig_width_inches = 7. # TOTAL figure width, in inches
    # (i.e., 8x11.5 paper with 1/2-inch margins)
margin_inches = 1./8. # margin width in inches (for both x and y) and 
    # horizontally in between figures
margin_bottom_inches = 0.75*(2 - (rbcz is None)) 
    # larger bottom margin to make room for colorbar(s)
margin_top_inches = 1 # wider top margin to accommodate subplot titles AND metadata
margin_subplot_top_inches = 1/4 # margin to accommodate just subplot titles
ncol = 3 # put three plots per row
nrow = 1

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

field_components = [vr, vt, vp]
titles = [r'$v_r$', r'$v_\theta$', r'$v_\phi$']
units = r'$\rm{m}\ \rm{s}^{-1}$'
minmax = [minmaxvr, minmaxvt, minmaxvp]
minmaxrz = [minmaxvrrz, minmaxvtrz, minmaxvprz]

# Generate the actual figure of the correct dimensions
fig = plt.figure(figsize=(fig_width_inches, fig_height_inches))

for iplot in range(3):
    ax_left = margin_x + (iplot%ncol)*(subplot_width + margin_x)
    ax_bottom = 1 - margin_top - subplot_height - margin_subplot_top -\
            (iplot//ncol)*(subplot_height + margin_subplot_top +\
            margin_bottom)
    ax = fig.add_axes((ax_left, ax_bottom, subplot_width, subplot_height))
    plot_azav (field_components[iplot], rr, cost, fig=fig, ax=ax,\
        units=units, minmax=minmax[iplot], rbcz=rbcz,\
        minmaxrz=minmaxrz[iplot], plotlatlines=plotlatlines,\
           plotcontours=plotcontours, plotboundary=plotboundary)
    ax.set_title(titles[iplot], va='bottom', **csfont)

    # Mark radii if desired
    if not rvals is None:
        for rval in rvals:
            rval_n = rval/ro
            plt.plot(rval_n*sint, rval_n*cost, 'k--', linewidth=0.5)

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
fig.text(margin_x, 1 - 0.3*margin_top, 'Velocity field (zonally averaged)',\
         ha='left', va='top', fontsize=fsize, **csfont)

fig.text(margin_x, 1 - 0.5*margin_top, time_string,\
         ha='left', va='top', fontsize=fsize, **csfont)

savename = dirname_stripped + '_v_azav_' + str(iter1).zfill(8) + '_' +\
        str(iter2).zfill(8) + '.png'
if saveplot:
    print ("Saving azimuthal average of velocity field in " + plotdir +\
            savename + ' ...')
    plt.savefig(plotdir + savename, dpi=300)
if showplot:
    plt.show()
plt.close()
