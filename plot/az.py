# Author: Loren Matilsky
# Created: 06/22/2019
# This script plots any user-specified quantitities (or by default, the
# three velocity components) from the AZ_Avgs, in the meridional plane 
# ...for the Rayleigh run directory indicated by [dirname]. To use an 
# AZ_Avgs file  different than the one associated with the longest 
# averaging time interval, use
# -usefile [complete name of desired AZ_Avgs file]
# Does not save the plot by default

import numpy as np
import pickle
import matplotlib as mpl
mpl.use('TkAgg')
import matplotlib.pyplot as plt
plt.rcParams['mathtext.fontset'] = 'dejavuserif'
csfont = {'fontname':'DejaVu Serif'}
import sys, os
sys.path.append(os.environ['rapp'])
sys.path.append(os.environ['co'])
sys.path.append(os.environ['pl'])
from azav_util import plot_azav
from common import get_widest_range_file, get_dict, strip_dirname

# Get directory name and stripped_dirname for plotting purposes
dirname = sys.argv[1]
dirname_stripped = strip_dirname(dirname)

# Directory with data and plots, make the plotting directory if it doesn't
# already exist    
datadir = dirname + '/data/'
plotdir = dirname + '/plots/'
if (not os.path.isdir(plotdir)):
    os.makedirs(plotdir)

# Set defaults and read command-line arguments (CLAs)
showplot = True
saveplot = False
plotcontours = True
plotlatlines = False
rvals = None # radii to mark on the meridional plane
minmax = None # if specified, must give minmax pair for each quantity 
posdef = None # 1 for each quantity
logscale = None # 1 for each quantity
qv = [1, 2, 3] # by default plot the velocity components
ncol = 3 # in the figure, put three plots per row
AZ_Avgs_file = get_widest_range_file(datadir, 'AZ_Avgs')

args = sys.argv[2:]
nargs = len(args)
for i in range(nargs):
    arg = args[i]
    if arg == '-qv':
        qv_str = args[i+1].split()
        print(qv_str)
        qv = []
        for qv_str_elem in qv_str:
            qv.append(int(qv_str_elem))
    elif arg == '-minmax':
        minmax_str = args[i+1].split()
        minmax = []
        for j in range(len(minmax_str)):
            minmax[j] = float(minmax_str[j])
    elif arg == '-logscale':
        logscale_str = args[i+1].split()
        logscale = []
        for j in range(len(logscale_str)):
            logscale[j] = bool(logscale_str[j])
    elif arg == '-posdef':
        posdef_str = args[i+1].split()
        posdef = []
        for j in range(len(posdef_str)):
            posdef[j] = bool(posdef_str[j])
    elif arg == '-noshow':
        showplot = False
    elif arg == '-save':
        saveplot = True
        savename = args[i+1]
    elif arg == '-nocontour':
        plotcontours = False
    elif arg == '-plotlat':
        plotlatlines = True
    elif arg == '-rvals':
        rvals_str = args[i+1].split()
        rvals = []
        for j in range(len(rvals_str)):
            rvals[j] = float(rvals_str[j])
    elif arg == '-ncol':
        ncol = int(args[i+1])
    elif arg == '-usefile':
        AZ_Avgs_file = args[i+1]
        AZ_Avgs_file = AZ_Avgs_file.split('/')[-1]

# Get the AZ_Avg data
print ('AZ_Avgs data from ' + datadir + AZ_Avgs_file + ' ...')
di = get_dict(datadir + AZ_Avgs_file)

iter1, iter2 = di['iter1'], di['iter2']
vals = di['vals']
lut = di['lut']

# Get necessary grid info
rr = di['rr']
cost = di['cost']
sint = di['sint']
tt_lat = di['tt_lat']
xx = di['xx']

ind_pp = lut[1801]
ind_mm = lut[1802]
ind_cor = lut[1803]
ind_visc = lut[1804]

torque_rs, torque_mc, torque_visc = -vals[:, :, ind_pp],\
        -vals[:, :, ind_mm] + vals[:, :, ind_cor],\
        vals[:, :, ind_visc]
torque_tot = torque_rs + torque_mc + torque_visc

max_sig = max(np.std(torque_rs), np.std(torque_mc), np.std(torque_visc))

# Set up the actual figure from scratch
fig_width_inches = 7. # TOTAL figure width, in inches
    # (i.e., 8x11.5 paper with 1/2-inch margins)
margin_inches = 1./8. # margin width in inches (for both x and y) and 
    # horizontally in between figures
margin_top_inches = 1. # wider top margin to accommodate subplot titles AND metadata
margin_subplot_top_inches = 1. # margin to accommodate just subplot titles
nplots = len(qv)
nrow = np.int(np.ceil(nplots/3))

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

# Generate the actual figure of the correct dimensions
fig = plt.figure(figsize=(fig_width_inches, fig_height_inches))

for iplot in range(nplots):
    ax_left = margin_x + (iplot%ncol)*(subplot_width + margin_x)
    ax_bottom = 1 - margin_top - subplot_height - \
            (iplot//ncol)*(subplot_height + margin_subplot_top)
    ax = fig.add_axes((ax_left, ax_bottom, subplot_width, subplot_height))
    iq = qv[iplot]
    field = vals[:, :, lut[iq]]
    if not minmax is None:
        this_minmax = minmax[2*iplot], minmax[2*iplot + 1]
    else:
        this_minmax = None
    plot_azav (field, rr, cost, sint, fig=fig, ax=ax, minmax=this_minmax,\
            plotcontours=plotcontours, plotlatlines=plotlatlines,\
            rvals=rvals)
    ax.set_title('iq = %i' %iq, verticalalignment='bottom', **csfont)

# Put some metadata in upper left
fsize = 12
fig.text(margin_x, 1. - 0.1*margin_top, dirname_stripped,\
         ha='left', va='top', fontsize=fsize, **csfont)
fig.text(margin_x, 1 - 0.4*margin_top,\
         str(iter1).zfill(8) + ' to ' + str(iter2).zfill(8),\
         ha='left', va='top', fontsize=fsize, **csfont)

if saveplot:
    savefile = plotdir + savename 
    print ('Saving plot at ' + savefile + ' ...')
    plt.savefig(savefile, dpi=300)
if showplot:
    plt.show()
plt.close()
