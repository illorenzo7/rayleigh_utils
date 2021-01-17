# Author: Loren Matilsky
# Date created: 03/02/2019
import matplotlib as mpl
mpl.use('TkAgg')
import matplotlib.pyplot as plt
plt.rcParams['mathtext.fontset'] = 'dejavuserif'
csfont = {'fontname':'DejaVu Serif'}
import numpy as np
import sys, os
sys.path.append(os.environ['raco'])
from common import *
from plotcommon import axis_range
from tl_util import plot_tl

# Get the run and data directories
dirname = sys.argv[1]
dirname_stripped = strip_dirname(dirname)

# domain bounds
ncheby, domain_bounds = get_domain_bounds(dirname)
ri = np.min(domain_bounds)
ro = np.max(domain_bounds)
d = ro - ri
datadir = dirname + '/data/'

# Find the time/latitude file(s) the data directory. If there are 
# multiple, by default choose the one with widest range in the trace.
the_file = get_widest_range_file(datadir, 'time-latitude')

# more defaults
minmax = None
xminmax = None
xmin = None
xmax = None
saveplot = None # turned off by default if saving one figure, can change
# with -save option
showplot = False # only show if plotting one figure
labelbytime = False # by default label by first/last iteration number
# not first/last time

rvals = 'all'  # by default, plot all available time-lat levels
    # user specifies another choice via -rvals '[val1] [val2] ... [valn]'
    # where 'vals' have dimensional units in cm: 4.8e10, 5e10, etc.
irvals = None # user can also specify -irvals '2 3 9', etc.
navg = 1 # by default average over 1 AZ_Avgs instance (no average)
# for navg > 1, a "sliding average" will be used.
tag = '' # optional way to tag save directory
lats = [0.]
plottimes = None

# Get command-line arguments
args = sys.argv[2:]
nargs = len(args)
for i in range(nargs):
    arg = args[i]
    if arg == '-minmax':
        strings = args[i+1].split()
        minmax = []
        for st in strings:
            minmax.append(float(st))
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
    elif arg == '-navg':
        navg = int(args[i+1])
        if navg % 2 == 0:
            print ("Please don't enter even values for navg!")
            print ("Replacing navg = %i with navg = %i" %(navg, navg + 1))
            navg += 1
    elif arg == '-xminmax':
        xminmax = float(args[i+1]), float(args[i+2])
    elif arg == '-xmin':
        xmin = float(args[i+1])
    elif arg == '-xmax':
        xmax = float(args[i+1])
    elif arg == '-save':
        saveplot = True
    elif arg == '-tlabel':
        labelbytime = True
    elif arg == '-tag':
        tag = '_' + args[i+1]
    elif arg == '-lats':
        strings = args[i+1].split()
        for string in strings:
            lats.append(float(string))
    elif arg == '-times':
        strings = args[i+1].split()
        plottimes = []
        for string in strings:
            plottimes.append(float(string))

# Get plot directory and create if not already there
plotdir = dirname + '/plots/time-lat' + tag + '_v/'
if labelbytime:
    plotdir = dirname + '/plots/time-lat' + tag + '_v_tlabel/'
if not os.path.isdir(plotdir):
    os.makedirs(plotdir)

# Read in the time-latitude data (dictionary form)
print ('Getting time-latitude trace from ' + datadir + the_file)
di = get_dict(datadir + the_file)
vals = di['vals']
times = di['times']
iters = di['iters']
rr = di['rr']
irvals_avail = di['rinds']
rvals_avail = rr[irvals_avail]
ri = di['ri']; ro = di['ro']; shell_depth = ro - ri
tt_lat = di['tt_lat']
ntheta = di['ntheta']

qvals = np.array(di['qvals'])

niter = di['niter']
nr = di['nr']
nq = di['nq']

iter1 = di['iter1']
iter2 = di['iter2']

# Get the baseline time unit
rotation = get_parameter(dirname, 'rotation')
if rotation:
    time_unit = compute_Prot(dirname)
    time_label = r'$\rm{P_{rot}}$'
else:
    time_unit = compute_tdt(dirname)
    time_label = r'$\rm{TDT}$'

vr_index = np.argmin(np.abs(qvals - 1))
vt_index = np.argmin(np.abs(qvals - 2))
vp_index = np.argmin(np.abs(qvals - 3))

# determine desired levels to plot
if irvals is None:
    if rvals == 'all':
        irvals = np.arange(len(irvals_avail))
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

# Get raw traces of vr, vtheta, vphi (in m/s)
vr = vals[:, :, :, vr_index]/100.
vt = vals[:, :, :, vt_index]/100.
vp = vals[:, :, :, vp_index]/100.

# Normalize the time 
times /= time_unit

# Make meshgrid of time/latitude
# Take into account if user specified xmin, xmax
if xminmax is None:
    xminmax = np.min(times), np.max(times)
# Change JUST xmin or xmax, if desired
if not xmin is None:
    xminmax = xmin, xminmax[1]
if not xmax is None:
    xminmax = xminmax[0], xmax

it1 = np.argmin(np.abs(times - xminmax[0]))
it2 = np.argmin(np.abs(times - xminmax[1]))
t1, t2 = times[it1], times[it2] # These begin times and end times
        # will be used for labeling the plots

# set figure dimensions
suvplot_width_inches = 6.5
margin_inches = 1./4.
margin_bottom_inches = 1./2. # space for x-axis label
margin_top_inches = 1./2.
margin_left_inches = 5./8. # space for latitude label
margin_right_inches = 0.9

fig_width_inches = suvplot_width_inches + margin_right_inches +\
        margin_left_inches
suvplot_height_inches = 2.0

nrow = 3
fig_height_inches = nrow*suvplot_height_inches +\
        (nrow - 1)*margin_inches + margin_bottom_inches +\
        margin_top_inches

margin_x = margin_inches/fig_width_inches
margin_y = margin_inches/fig_height_inches
suvplot_width = suvplot_width_inches/fig_width_inches
suvplot_height = suvplot_height_inches/fig_height_inches
margin_left = margin_left_inches/fig_width_inches
margin_bottom = margin_bottom_inches/fig_height_inches
margin_top = margin_top_inches/fig_height_inches

# field units and labels
units = r'$\rm{m\ s^{-1}}$'
labels = [r'$\langle v_r\rangle$', r'$\langle v_\theta\rangle$',\
        r'$\langle v_\phi\rangle$']

# Loop over the desired radii and save plots
for i in range(len(irvals)):
    ir = irvals[i]
    rval = rvals_avail[ir]/rsun 
    print('plotting r/rsun = %0.3f (ir = %02i)' %(rval, ir))
    vr_loc = vr[:, :, ir]
    vt_loc = vt[:, :, ir]
    vp_loc = vp[:, :, ir]
    
    # Make appropriate file name to save
    if labelbytime:
        savename = dirname_stripped + '_time-lat_v_' +\
                ('Prot%05.0f-to-%05.0f_' %(t1, t2)) +\
            ('rval%0.3f' %rval) + '.png'
    else:
        savename = dirname_stripped + '_time-lat_v_' +\
                ('%08i_%08i_' %(iter1, iter2)) +\
            ('rval%0.3f' %rval) + '.png'

    if minmax is None:
        minmax_vr = None
        minmax_vt = None
        minmax_vp = None
    else:
        if len(minmax) == 2:
            minmax_vr = minmax
            minmax_vt = minmax
            minmax_vp = minmax
        elif len(minmax) == 6:
            minmax_vr = minmax[0], minmax[1]
            minmax_vt = minmax[2], minmax[3]
            minmax_vp = minmax[4], minmax[5]
        else:
            print ("error: minmax must have length 2 or 6, as in")
            print ("-minmax '-5e3 5e3 -1e4 1e4 -5e4 5e4'")
            print ("exiting")
            sys.exit()

    # Create figure with  3 panels in a row (time-radius plots of
    #       vr, vtheta, and vphi)
    fig = plt.figure(figsize=(fig_width_inches, fig_height_inches))
    ax1 = fig.add_axes((margin_left, 1. - margin_top -\
            suvplot_height - 0*(suvplot_height + margin_y),\
            suvplot_width, suvplot_height))
    ax2 = fig.add_axes((margin_left, 1. - margin_top -\
            suvplot_height - 1*(suvplot_height + margin_y),\
            suvplot_width, suvplot_height))
    ax3 = fig.add_axes((margin_left, 1. - margin_top -\
            suvplot_height - 2*(suvplot_height + margin_y),\
            suvplot_width, suvplot_height))

    # Plot evolution of each (zonally averaged) field component
    plot_tl(vr_loc, times, tt_lat, fig=fig, ax=ax1, navg=navg,\
            minmax=minmax_vr, units=units, xminmax=xminmax, yvals=lats,\
            plottimes=plottimes)
    plot_tl(vt_loc, times, tt_lat, fig=fig, ax=ax2, navg=navg,\
            minmax=minmax_vt, units=units, xminmax=xminmax, yvals=lats,\
            plottimes=plottimes)
    plot_tl(vp_loc, times, tt_lat, fig=fig, ax=ax3, navg=navg,\
            minmax=minmax_vp, units=units, xminmax=xminmax, yvals=lats,\
            plottimes=plottimes)

    # Label with the field components
    for irow in range(nrow):
        label = labels[irow]
        fig.text(margin_left + 0.5*margin_x, 1. - margin_top - \
                0.5*margin_y - irow*(suvplot_height + margin_y), label,\
                va='top', ha='left', fontsize=14,\
                bbox=dict(facecolor='white'))

    # Turn the x tick labels off for the top two strips
    ax1.set_xticklabels([])
    ax2.set_xticklabels([])

    # Label x (time) axis
    ax3.set_xlabel('time (' + time_label + ')', **csfont)
    # Label y-axis (latitude in degrees)
    ax2.set_ylabel('latitude (deg)', **csfont)

    # Put some useful information on the title
    averaging_time = (times[-1] - times[0])/niter*navg
    title = dirname_stripped + '     ' + (r'$r/R_\odot\ =\ %0.3f$' %rval)
    if navg > 1:
        title += '     ' + ('t_avg = %.1f Prot' %averaging_time)
    else:
        title += '     t_avg = none'
    ax1.set_title(title, **csfont)

    # Save the plot
    if saveplot:
        print ('Saving the time-latitude plot at ')
        print (plotdir + savename)
        print ("=======================================")
        plt.savefig(plotdir + savename, dpi=200)

    # Show the plot if only plotting at one latitude
    if showplot:
        plt.show()
    plt.close()
