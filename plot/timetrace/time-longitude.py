# Author: Loren Matilsky
# Date created: 08/18/2019
#import matplotlib as mpl
#mpl.use('TkAgg')
import matplotlib.pyplot as plt
plt.rcParams['mathtext.fontset'] = 'dejavuserif'
csfont = {'fontname':'DejaVu Serif'}
import numpy as np
import pickle
import sys, os
sys.path.append(os.environ['raco'])
sys.path.append(os.environ['rapl'])
from common import *
from plotcommon import *

# Get the run directory on which to perform the analysis
dirname = sys.argv[1]

# Data and plot directories
datadir = dirname + '/data/'
dirname_stripped = strip_dirname(dirname)

# domain bounds
ncheby, domain_bounds = get_domain_bounds(dirname)
ri = np.min(domain_bounds)
ro = np.max(domain_bounds)
d = ro - ri

# Get global rotation rate; this script fails for non-rotating models
angular_velocity = get_parameter(dirname, 'angular_velocity')
Prot = 2*np.pi/angular_velocity
Om0_nhz = angular_velocity/(2*np.pi)*1e9
print('Om0: ', Om0_nhz, ' nHz')

# more defaults
minmax = None # this may be true at first
minmax_wasnone = True # This always will be true unless user specifies
            # values through -minmax
# By default (if tminmax is None) will plot over whole time trace
tminmax = None
saveplot = None # turned off by default if saving one figure, can change
# with -save option
showplot = False # will only show if plotting one figure
labelbytime = False # by default label by first/last iteration number

rvals = 'all'  # by default, plot all available time-lat levels
    # user specifies another choice via -rvals '[val1] [val2] ... [valn]'
    # where 'vals' have dimensional units in cm: 4.8e10, 5e10, etc.
irvals = None # user can also specify -irvals '2 3 9', etc.

# By default, plot longitudinal B field
qval = 803

Om_subtract = None # by default, do not subtract any differential rotation
# if nonzero, subtract a CONSTANT DR

# Default no. rotations to plot per inch (vertically)
rpi = 100.

# Get command-line arguments
plotdir = None

args = sys.argv[2:]
nargs = len(args)
tag = ''
plottimes = None
clat = 10. 
dlat = 0.
the_file = None
for i in range(nargs):
    arg = args[i]
    if arg == '-plotdir':
        plotdir = args[i+1]
    if arg == '-minmax':
        minmax = float(args[i+1]), float(args[i+2])
        minmax_wasnone = False
    elif arg == '-usefile':
        the_file = args[i+1]
        the_file = the_file.split('/')[-1]
    elif arg == '-rvals':
        strings = args[i+1].split()
        rvals = []
        for j in range(len(strings)):
            rvals.append(float(strings[j]))
    elif arg == '-irvals':
        irvals = []
        strings = args[i+1].split()
        for j in range(len(strings)):
            irvals.append(int(strings[j]))
    elif arg == '-tminmax':
        tminmax = float(args[i+1]), float(args[i+2])
    elif arg == '-save':
        saveplot = True
    elif arg == '-tlabel':
        labelbytime = True
    elif arg == '-rpi': # rpi = rotations per inch
        rpi = float(args[i+1])
    elif arg == '-qval':
        qval = int(args[i+1])
    elif arg == '-tag':
        tag = args[i+1] + '_'
    elif arg == '-clat':
        clat = float(args[i+1])
    elif arg == '-dlat':
        dlat = float(args[i+1])
    elif arg == '-om': # enter this value in nHz, in the LAB frame
        Om_subtract = float(args[i+1]) - Om0_nhz
        Om_subtract *= (2*np.pi/1e9) # convert nHz --> rad s^-1
    elif arg == '-times':
        strings = args[i+1].split()
        plottimes = []
        for string in strings:
            plottimes.append(float(string))

# Get plot directory and create if not already there
if plotdir is None:
    plotdir = dirname + '/plots/time-lon' + tag + '/'
    if labelbytime:
        plotdir = dirname + '/plots/time-lon' + tag + '_tlabel/'
    if not os.path.isdir(plotdir):
        os.makedirs(plotdir)

# Find the time-longitude file(s) the data directory, for clat and dlat. 
# If there are 
# multiple, by default choose the one with widest range in the trace.
if the_file is None:
    if clat >= 0.:
        hemisphere = 'N'
    else:
        hemisphere = 'S'
    the_file = get_widest_range_file(datadir, 'time_lon_clat' + hemisphere +\
            '%02.0f_dlat%03.0f' %(np.abs(clat), dlat))

# Read in the time-longitude data (dictionary form)
print ('Getting time-longitude trace from ' + the_file)
di = get_dict(the_file)

# this is clat + hemisphere we ended up with
clat = di['clat']
if clat >= 0.:
    hemisphere = 'N'
else:
    hemisphere = 'S'

vals = di['vals']
times = di['times']
iters = di['iters']
qvals = di['qvals']
rvals_avail = di['rvals']
clat = di['clat']
dlat = di['dlat']

# get grid info
di_grid = get_grid_info(dirname)

lons = di_grid['lons']
nphi = di_grid['nphi']

# baseline time unit
iter1, iter2 = get_iters_from_file(the_file)
time_unit, time_label, rotation = get_time_unit(dirname)

# determine desired levels to plot
if irvals is None:
    if rvals == 'all':
        irvals = np.arange(len(rvals_avail))
    else:
        irvals = []
        for rval in rvals:
            ir = np.argmin(np.abs(rvals_avail - rval))
            irvals.append(ir)

if saveplot is None:
    if len(irvals) == 1:
        saveplot = False
    else:
        saveplot = True
if len(irvals) == 1:
    showplot = True

# Get raw trace of "qval"
if tminmax is None:
    it1, it2 = 0, len(times) - 1
else:
    it1 = np.argmin(np.abs(times/time_unit - tminmax[0]))
    it2 = np.argmin(np.abs(times/time_unit - tminmax[1]))

if qval == 4001: # shorthand for b_r * b_phi
    qind1 = np.argmin(np.abs(qvals - 801))
    qind2 = np.argmin(np.abs(qvals - 803))
    quant = vals[it1:it2+1, :, :, qind1]*vals[it1:it2+1, :, :, qind2]
elif qval == 4002: # shorthand for b_theta * b_phi
    qind1 = np.argmin(np.abs(qvals - 802))
    qind2 = np.argmin(np.abs(qvals - 803))
    quant = vals[it1:it2+1, :, :, qind1]*vals[it1:it2+1, :, :, qind2]
else:
    qind = np.argmin(np.abs(qvals - qval))
    quant = vals[it1:it2+1, :, :, qind]

times = times[it1:it2+1]/Prot
t1, t2 = times[0], times[-1] # These begin times and end times
        # will be used for labeling the plots

# Make meshgrid of time/radius
lons_2d, times_2d = np.meshgrid(lons, times, indexing='ij')

# Subtract DR, if desired
if not Om_subtract is None:
    phi_deflections = (times*Prot*Om_subtract*180./np.pi) % 360.
    for it in range(len(times)):
        phi_deflection = phi_deflections[it]
        nroll = int(phi_deflection/360*nphi)
        quant[it, :, :] = np.roll(quant[it, :, :], -nroll, axis=0)

# Create figure dimensions
# Create figure as a vertical strip, taking into account the desired
# rpi
ncol = 1
nplots = ncol
nrow = 1
#nrow = np.int(np.ceil(nplots/ncol))
fig_width_inches = 5.
margin_inches = 1./16. # margin width in inches (for both x and y) and 
    # horizontally in between figures
margin_right_inches = 1.
margin_left_inches = 1.
margin_bottom_inches = 0.5
# larger bottom margin to make room for colorbar(s)
margin_top_inches = 1. # wider top margin to accommodate subplot titles AND metadata
subplot_height_inches = (times[-1] - times[0])/rpi
fig_height_inches = margin_top_inches + nrow*(subplot_height_inches +\
        margin_bottom_inches)
print('figsize = ', fig_width_inches, fig_height_inches)

subplot_width_inches = (fig_width_inches - margin_left_inches -\
        margin_right_inches - (ncol - 1)*margin_inches)/ncol
    # Make the subplot width so that ncol subplots fit together side-by-side
    # with margins in between them and at the left and right.

# "Margin" in "figure units"; figure units extend from 0 to 1 in BOTH 
# directions, so unitless dimensions of margin will be different in x and y
# to force an equal physical margin
margin_x = margin_inches/fig_width_inches
margin_y = margin_inches/fig_height_inches
margin_top = margin_top_inches/fig_height_inches
margin_bottom = margin_bottom_inches/fig_height_inches
margin_right = margin_right_inches/fig_width_inches
margin_left = margin_left_inches/fig_width_inches

# Subplot dimensions in figure units
subplot_width = subplot_width_inches/fig_width_inches
subplot_height = subplot_height_inches/fig_height_inches

# Loop over the desired radii and save plots
for i in range(len(irvals)):
    ir = irvals[i]
    rval = rvals_avail[ir]
    print('plotting r/rsun = %0.3f (ir = %02i)' %(rval, ir))
    
    quant_loc = quant[:, :, ir]

    # Make appropriate file name to save
    if labelbytime:
        savename = dirname_stripped + '_time-lon_' +\
                ('Prot%05.0f-to-%05.0f_clat%s%02.0f_dlat%02.0f_'\
                %(t1, t2, hemisphere, np.abs(clat), dlat)) +\
                ('qv%04i_' %qval) + ('rval%0.3f' %rval) + '.png'
    else:
        savename = dirname_stripped + '_time-lon_' +\
                ('%08i_%08i_clat%s%02.0f_dlat%02.0f_'\
                %(iter1, iter2, hemisphere, np.abs(clat), dlat)) +\
                ('qv%04i_' %qval) + ('rval%0.3f' %rval) + '.png'

    if minmax is None or minmax_wasnone:
        std_quant = np.std(quant_loc)
        minmax = -3.*std_quant, 3.*std_quant
        minmax_wasnone = True

    # Generate the actual figure of the correct dimensions
    fig = plt.figure(figsize=(fig_width_inches, fig_height_inches))

    for iplot in range(nplots):
        ax_left = margin_left + (iplot%ncol)*(subplot_width + margin_x)
        ax_bottom = 1. - margin_top - subplot_height
        ax = fig.add_axes((ax_left, ax_bottom, subplot_width, subplot_height))
     
        im = ax.pcolormesh(lons_2d, times_2d, quant_loc.T,\
                vmin=minmax[0], vmax=minmax[1], cmap='RdYlBu_r')
        ax.invert_yaxis()

        # Put colorbar next to plot 
        # First make room and then find location of subplot

        ax_xmin, ax_xmax, ax_ymin, ax_ymax = axis_range(ax)
        ax_delta_x = ax_xmax - ax_xmin
        ax_delta_y = ax_ymax - ax_ymin
        ax_center_x = ax_xmin + 0.5*ax_delta_x

        cbar_width = 0.15*(1 - ax_xmax)
        cbar_height = 0.4*ax_delta_y

        cbar_left = ax_xmax + 0.1*(1 - ax_xmax)
        cbar_bottom = ax_ymin + 0.5*ax_delta_y - 0.5*cbar_height
        cax = fig.add_axes((cbar_left, cbar_bottom, cbar_width, cbar_height))
        cax.set_title('cgs', **csfont)
        plt.colorbar(im, cax=cax)

        # Set x label
        ax.set_xlabel('longitude (deg)', fontsize=12)

        # Label y (time) axis
        timeunit = r'$$'
        ylabel = r'$\rm{time}\ (P_{\rm{rot}})}$'
    #    ylabel = r'$\rm{time}\ (\longleftarrow)}$'
        ax.set_ylabel(ylabel, **csfont, fontsize=12)

        ax.set_xlim((lons[0], lons[-1]))

        # Label the plots by B_r, B_theta, B_phi
        ax_xmin, ax_xmax = ax.get_xlim()
        ax_Dx = ax_xmax - ax_xmin
        ax.text(ax_xmin + 1.01*ax_Dx, 0.,  'qval = %i' %qval, **csfont)

        # Put some useful information on the title
        title = dirname_stripped + (', iq = %i' %qval) + '\n' +\
                (r'$r/R_\odot\ =\ %0.3f$' %rval) + '\n' +\
                (r'$\lambda_c=%02.0f^\circ\ \ \ \ \ \Delta\lambda=%02.0f^\circ$'\
                %(clat, dlat)) 
        ax.set_title(title, **csfont)

        # Get ticks everywhere
        plt.sca(ax)
        plt.minorticks_on()
        plt.tick_params(top=True, right=True, direction='in', which='both')

        # Put grid of white lines on plot
        ax.grid(color='k', linestyle='-', linewidth=1., alpha=0.3)
        ax.grid(which='minor', color='k', linestyle='-', linewidth=0.5, alpha=0.3)

    #Save the plot
    if saveplot:
        print ('Saving the time-latitude plot at ' + plotdir +\
                savename)
        plt.savefig(plotdir + savename, dpi=300)

    # Show the plot if only plotting at one latitude
    if showplot:
        plt.show()
    plt.close()