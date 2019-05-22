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
from binormalized_cbar import MidpointNormalize
import sys, os
sys.path.append(os.environ['rapp'])
sys.path.append(os.environ['co'])
from azav_util import plot_azav
from common import get_widest_range_file, strip_dirname, get_file_lists,\
        get_desired_range, get_dict
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
minmax = None
AZ_Avgs_file = get_widest_range_file(datadir, 'AZ_Avgs')

args = sys.argv[2:]
nargs = len(args)
for i in range(nargs):
    arg = args[i]
    if arg == '-minmax':
        min_vr, max_vr  = float(args[i+1]), float(args[i+2])
        min_vt, max_vt  = float(args[i+3]), float(args[i+4])
        min_vp, max_vp  = float(args[i+5]), float(args[i+6])
    elif arg == '-noshow':
        showplot = False
    elif arg == '-nosave':
        saveplot = False
    elif arg == '-nocontour':
        plotcontours = False
    elif (arg == '-usefile'):
        AZ_Avgs_file = args[i+1]
        AZ_Avgs_file = AZ_Avgs_file.split('/')[-1]
       
# Read in AZ_Avgs data
print ('Getting data from ' + datadir + AZ_Avgs_file + ' ...')
di = get_dict(datadir + AZ_Avgs_file)

iter1, iter2 = di['iter1'], di['iter2']
vals = di['vals']
lut = di['lut']

# Grid info
rr = di['rr']
cost = di['cost']
sint = di['sint']

vr, vt, vp = vals[:, :, lut[1]]/100, vals[:, :, lut[2]]/100,\
        vals[:, :, lut[3]]

if minmax is None:
    nstd = 5.
    min_vr, max_vr = -nstd*np.std(vr), nstd*np.std(vr)
    min_vt, max_vt = -nstd*np.std(vt), nstd*np.std(vt)
    min_vp, max_vp = -nstd*np.std(vp), nstd*np.std(vp)

# Set up the actual figure from scratch
fig_width_inches = 7. # TOTAL figure width, in inches
    # (i.e., 8x11.5 paper with 1/2-inch margins)
margin_inches = 1./8. # margin width in inches (for both x and y) and 
    # horizontally in between figures
margin_top_inches = 1. # wider top margin to accommodate subplot titles AND metadata
margin_subplot_top_inches = 1. # margin to accommodate just subplot titles
ncol = 3 # put three plots per row
nrow = 1

subplot_width_inches = (fig_width_inches - (ncol + 1)*margin_inches)/ncol
    # Make the subplot width so that ncol subplots fit together side-by-side
    # with margins in between them and at the left and right.
subplot_height_inches = 2*subplot_width_inches # Each subplot should have an
    # aspect ratio of y/x = 2/1 to accommodate meridional planes. 
fig_height_inches = nrow*subplot_height_inches + margin_top_inches +\
    (nrow - 1)*margin_subplot_top_inches + margin_inches 
    # Room for titles on each row and a regular margin on the bottom
fig_aspect = fig_height_inches/fig_width_inches

# "Margin" in "figure units"; figure units extend from 0 to 1 in BOTH 
# directions, so unitless dimensions of margin will be different in x and y
# to force an equal physical margin
margin_x = margin_inches/fig_width_inches
margin_y = margin_inches/fig_height_inches
margin_top = margin_top_inches/fig_height_inches
margin_subplot_top = margin_subplot_top_inches/fig_height_inches

# Subplot dimensions in figure units
subplot_width = subplot_width_inches/fig_width_inches
subplot_height = subplot_height_inches/fig_height_inches

field_components = [vr, vt, vp]
titles = [r'$v_r$', r'$v_\theta$', r'$v_\phi$']
units = r'$\rm{m}\ \rm{s}^{-1}$'
my_mins = [min_vr, min_vt, min_vp]
my_maxes = [max_vr, max_vt, max_vp]

# Generate the actual figure of the correct dimensions
fig = plt.figure(figsize=(fig_width_inches, fig_height_inches))

for iplot in range(3):
    ax_left = margin_x + (iplot%ncol)*(subplot_width + margin_x)
    ax_bottom = 1 - margin_top - subplot_height - \
            (iplot//ncol)*(subplot_height + margin_subplot_top)
    ax = fig.add_axes((ax_left, ax_bottom, subplot_width, subplot_height))
    plot_azav (field_components[iplot], rr, cost, sint, fig=fig, ax=ax,\
        units=units, minmax=(my_mins[iplot], my_maxes[iplot]),\
           norm=MidpointNormalize(0.), plotcontours=plotcontours)
    ax.set_title(titles[iplot], verticalalignment='bottom', **csfont)

# Put some metadata in upper left
fsize = 12
fig.text(margin_x, 1 - 0.1*margin_top, dirname_stripped,\
         ha='left', va='top', fontsize=fsize, **csfont)
fig.text(margin_x, 1 - 0.3*margin_top, 'Velocity field (zonally averaged)',\
         ha='left', va='top', fontsize=fsize, **csfont)
iter1, iter2 = di['iter1'], di['iter2']
iter_string = str(iter1).zfill(8) + ' to ' + str(iter2).zfill(8) 

fig.text(margin_x, 1 - 0.5*margin_top, iter_string,\
         ha='left', va='top', fontsize=fsize, **csfont)

savename = dirname_stripped + '_v_azav_' + str(iter1).zfill(8) + '_' +\
        str(iter2).zfill(8) + '.png'
if saveplot:
    print ("Saving azimuthal average of B field in " + plotdir +\
            savename + ' ...')
    plt.savefig(plotdir + savename, dpi=300)
if showplot:
    plt.show()
plt.close()
