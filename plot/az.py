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
rbcz = None
plotlatlines = False
rvals = [] # radii to mark on the meridional plane
minmax = None # if specified, must give minmax pair for each quantity 
posdef = None # 1 for each quantity
logscale = None # 1 for each quantity
symlog = None # 1 for each quantity
qv = [1, 2, 3] # by default plot the velocity components
qv_sub = []
ncol = 3 # in the figure, put three plots per row
AZ_Avgs_file = get_widest_range_file(datadir, 'AZ_Avgs')

args = sys.argv[2:]
nargs = len(args)
for i in range(nargs):
    arg = args[i]
    if arg == '-qv':
        qv_str = args[i+1].split()
        qv = []
        for qv_str_elem in qv_str:
            qv.append(int(qv_str_elem))
    elif arg == '-sub':
        strings = args[i+1].split()
        qv_sub = []
        for st in strings:
            qv_sub.append(int(st))
    elif arg == '-rbcz':
        rbcz = float(args[i+1])
    elif arg == '-minmax':
        minmax_str = args[i+1].split()
        minmax = []
        for st in minmax_str:
            minmax.append(float(st))
    elif arg == '-log':
        logscale_str = args[i+1].split()
        if logscale_str == ['all']:
            logscale = 'all'
        else:
            logscale = []
            for st in logscale_str:
                logscale.append(my_bool(st))
    elif arg == '-posdef':
        posdef_str = args[i+1].split()
        if posdef_str == ['all']:
            posdef = 'all'
        else:
            posdef = []
            for st in posdef_str:
                posdef.append(my_bool(st))
    elif arg == '-symlog':
        symlog_str = args[i+1].split()
        if symlog_str == ['all']:
            symlog = 'all'
        else:
            symlog = []
            for st in symlog_str:
                symlog.append(my_bool(st))
    elif arg == '-noshow':
        showplot = False
    elif arg == '-save':
        saveplot = True
        savename = args[i+1]
    elif arg == '-nocontour':
        plotcontours = False
    elif arg == '-plotlat':
        plotlatlines = True
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
    elif arg == '-ncol':
        ncol = int(args[i+1])
    elif arg == '-usefile':
        AZ_Avgs_file = args[i+1]
        AZ_Avgs_file = AZ_Avgs_file.split('/')[-1]

# Get the AZ_Avg data
print ('Getting AZ_Avgs data from ' + datadir + AZ_Avgs_file + ' ...')
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

# Get necessary grid info
rr = di['rr']
cost = di['cost']
sint = di['sint']

# Set up the actual figure from scratch
fig_width_inches = 7. # TOTAL figure width, in inches
    # (i.e., 8x11.5 paper with 1/2-inch margins)
margin_inches = 1./8. # margin width in inches (for both x and y) and 
    # horizontally in between figures
margin_top_inches = 1 # wider top margin to accommodate subplot titles AND metadata
margin_bottom_inches = 0.75*(2 - (rbcz is None)) 
    # larger bottom margin to make room for colorbar(s)
margin_subplot_top_inches = 1/4 # margin to accommodate just subplot titles
nplots = len(qv) + len(qv_sub)//2
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
    ax_bottom = 1. - margin_top - margin_subplot_top -\
            subplot_height - (iplot//ncol)*(subplot_height +\
            margin_subplot_top + margin_bottom)
    ax = fig.add_axes((ax_left, ax_bottom, subplot_width, subplot_height))
    if iplot < len(qv):
        iq = qv[iplot]
        field = vals[:, :, lut[iq]]
        title = 'iq = %i' %iq
    else:
        iq1 = qv_sub[0]
        iq2 = qv_sub[1]
        field = vals[:, :, lut[iq1]] - vals[:, :, lut[iq2]]
        title = 'iq = %i - %i' %(iq1, iq2)
    if not minmax is None:
        this_minmax = minmax[2*iplot], minmax[2*iplot + 1]
    else:
        this_minmax = None

    if posdef is None:
        this_posdef = False
    elif posdef == 'all':
        this_posdef = True
    else:
        this_posdef = posdef[iplot]

    if logscale is None:
        this_logscale = False
    elif logscale == 'all':
        this_logscale = True
    else:
        this_logscale = logscale[iplot]

    if symlog is None:
        this_symlog = False
    elif symlog == 'all':
        this_symlog = True
    else:
        this_symlog = symlog[iplot]

    plot_azav (field, rr, cost, fig=fig, ax=ax, minmax=this_minmax,\
            plotcontours=plotcontours, plotlatlines=plotlatlines,\
            rvals=rvals, posdef=this_posdef, logscale=this_logscale,\
            symlog=this_symlog, rbcz=rbcz)
    ax.set_title(title, verticalalignment='bottom', **csfont)

# Label averaging interval
if rotation:
    time_string = ('t = %.1f to %.1f ' %(t1/time_unit, t2/time_unit))\
            + time_label + '\n' + (r'$\ (\Delta t = %.1f\ $'\
            %((t2 - t1)/time_unit)) + time_label + ')'
else:
    time_string = ('t = %.3f to %.3f ' %(t1/time_unit, t2/time_unit))\
            + time_label + (r'$\ (\Delta t = %.3f\ $'\
            %((t2 - t1)/time_unit)) + time_label + ')'

# Put some metadata in upper left
fsize = 12
line_height = 1./4./fig_height_inches
fig.text(margin_x, 1. - margin_y, dirname_stripped,\
         ha='left', va='top', fontsize=fsize, **csfont)
fig.text(margin_x, 1 - margin_y - line_height,\
        time_string, ha='left', va='top', fontsize=fsize, **csfont)

if saveplot:
    savefile = plotdir + savename + '.png'
    print ('Saving plot at ' + savefile + ' ...')
    plt.savefig(savefile, dpi=300)
if showplot:
    plt.show()
plt.close()
