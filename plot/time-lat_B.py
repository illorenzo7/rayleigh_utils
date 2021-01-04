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
from common import get_widest_range_file, strip_dirname, rsun, get_dict
from plotcommon import axis_range
from get_parameter import get_parameter
from time_scales import compute_Prot, compute_tdt
from tl_util import plot_tl

# Get the run and data directories
dirname = sys.argv[1]
dirname_stripped = strip_dirname(dirname)
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
lats = None

# Get command-line arguments
args = sys.argv[2:]
nargs = len(args)
for i in range(nargs):
    arg = args[i]
    if arg == '-minmax':
        try: # See if user wants to set ranges for B_r, B_theta, and B_phi
            minmax = float(args[i+1]), float(args[i+2]), float(args[i+3]),\
                    float(args[i+4]), float(args[i+5]), float(args[i+6])
        except:
            minmax = float(args[i+1]), float(args[i+2])
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
        lats_str = args[i+1].split()
        lats = []
        for lat_str in lats_str:
            lats.append(float(lat_str))

# Get plot directory and create if not already there
plotdir = dirname + '/plots/time-lat_B' + tag + '/'
if labelbytime:
    plotdir = dirname + '/plots/time-lat_B_tlabel' + tag + '/'
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

br_index = np.argmin(np.abs(qvals - 801))
bt_index = np.argmin(np.abs(qvals - 802))
bp_index = np.argmin(np.abs(qvals - 803))

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

# Get raw traces of br, btheta, bphi
br = vals[:, :, :, br_index]
bt = vals[:, :, :, bt_index]
bp = vals[:, :, :, bp_index]

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
subplot_width_inches = 6.5
margin_inches = 1./4.
margin_bottom_inches = 1./2. # space for x-axis label
margin_top_inches = 1./2.
margin_left_inches = 5./8. # space for latitude label
margin_right_inches = 0.9

fig_width_inches = subplot_width_inches + margin_right_inches +\
        margin_left_inches
subplot_height_inches = 2.0

nrow = 3
fig_height_inches = nrow*subplot_height_inches +\
        (nrow - 1)*margin_inches + margin_bottom_inches +\
        margin_top_inches

margin_x = margin_inches/fig_width_inches
margin_y = margin_inches/fig_height_inches
subplot_width = subplot_width_inches/fig_width_inches
subplot_height = subplot_height_inches/fig_height_inches
margin_left = margin_left_inches/fig_width_inches
margin_bottom = margin_bottom_inches/fig_height_inches
margin_top = margin_top_inches/fig_height_inches

# field units and labels
units = r'$\rm{G}$'
labels = [r'$\langle B_r\rangle$', r'$\langle B_\theta\rangle$',\
        r'$\langle B_\phi\rangle$']

# Loop over the desired radii and save plots
for i in range(len(irvals)):
    ir = irvals[i]
    rval = rvals_avail[ir]/rsun 
    print('plotting r/rsun = %0.3f (ir = %02i)' %(rval, ir))
    br_loc = br[:, :, ir]
    bt_loc = bt[:, :, ir]
    bp_loc = bp[:, :, ir]
    
    # Make appropriate file name to save
    if labelbytime:
        savename = dirname_stripped + '_time-lat_B_' +\
                ('Prot%05.0f-to-%05.0f_' %(t1, t2)) +\
            ('rval%0.3f' %rval) + '.png'
    else:
        savename = dirname_stripped + '_time-lat_B_' +\
                ('%08i_%08i_' %(iter1, iter2)) +\
            ('rval%0.3f' %rval) + '.png'

    if minmax is None:
        minmax_br = None
        minmax_bt = None
        minmax_bp = None
    else:
        if len(minmax) == 2:
            minmax_br = minmax
            minmax_bt = minmax
            minmax_bp = minmax
        elif len(minmax) == 6:
            minmax_br = minmax[0], minmax[1]
            minmax_bt = minmax[2], minmax[3]
            minmax_bp = minmax[4], minmax[5]

    # Create figure with  3 panels in a row (time-radius plots of
    #       br, btheta, and bphi)
    fig = plt.figure(figsize=(fig_width_inches, fig_height_inches))
    ax1 = fig.add_axes((margin_left, 1. - margin_top -\
            subplot_height - 0*(subplot_height + margin_y),\
            subplot_width, subplot_height))
    ax2 = fig.add_axes((margin_left, 1. - margin_top -\
            subplot_height - 1*(subplot_height + margin_y),\
            subplot_width, subplot_height))
    ax3 = fig.add_axes((margin_left, 1. - margin_top -\
            subplot_height - 2*(subplot_height + margin_y),\
            subplot_width, subplot_height))

    # Plot evolution of each (zonally averaged) field component
    plot_tl(br_loc, times, tt_lat, fig=fig, ax=ax1, navg=navg,\
            minmax=minmax_br, units=units, xminmax=xminmax, yvals=lats)
    plot_tl(bt_loc, times, tt_lat, fig=fig, ax=ax2, navg=navg,\
            minmax=minmax_bt, units=units, xminmax=xminmax, yvals=lats)
    plot_tl(bp_loc, times, tt_lat, fig=fig, ax=ax3, navg=navg,\
            minmax=minmax_bp, units=units, xminmax=xminmax, yvals=lats)

    # Label with the field components
    for irow in range(nrow):
        label = labels[irow]
        fig.text(margin_left + 0.5*margin_x, 1. - margin_top - \
                0.5*margin_y - irow*(subplot_height + margin_y), label,\
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
