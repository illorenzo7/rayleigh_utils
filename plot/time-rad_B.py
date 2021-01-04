# Author: Loren Matilsky
# Date created: 03/02/2019
import matplotlib as mpl
mpl.use('TkAgg')
import matplotlib.pyplot as plt
plt.rcParams['mathtext.fontset'] = 'dejavuserif'
csfont = {'fontname':'DejaVu Serif'}
import numpy as np
import pickle
import sys, os
sys.path.append(os.environ['raco'])
from common import get_file_lists, get_widest_range_file, strip_dirname,\
        rsun, get_dict
from get_parameter import get_parameter
from time_scales import compute_Prot, compute_tdt
from tl_util import plot_tl

# Get the run and directories
dirname = sys.argv[1]
dirname_stripped = strip_dirname(dirname)
datadir = dirname + '/data/'

# defaults, then CLAs

# Find the time/latitude file(s) the data directory. If there are 
# multiple, by default choose the one with widest range in the trace.
the_file = get_widest_range_file(datadir, 'time-radius')
minmax = None
xminmax = None
xmin = None
xmax = None
saveplot = True
showplot = True # will only show if plotting one figure

labelbytime = False # by default label by first/last iteration number
# not first/last time

# more defaults
desired_lats = [0.] # by default, plot time-radius diagram for fields 
    # at the equator
rbcz = None
navg = 1 # by default don't average in time
tag = '' # optional way to tag save directory
rvals = []

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
    elif arg == '-xminmax':
        xminmax = float(args[i+1]), float(args[i+2])
    elif arg == '-xmin':
        xmin = float(args[i+1])
    elif arg == '-xmax':
        xmax = float(args[i+1])
    elif arg == '-usefile':
        the_file = args[i+1]
        the_file = the_file.split('/')[-1]
    elif arg == '-lats':
        string_lats = args[i+1].split()
        if string_lats == ['all']:
            desired_lats = 'all'
        else:
            desired_lats = []
            for j in range(len(string_lats)):
                desired_lats.append(float(string_lats[j]))
    elif arg == '-navg':
        navg = int(args[i+1])
        if navg % 2 == 0:
            print ("Please don't enter even values for navg!")
            print ("Replacing navg = %i with navg = %i" %(navg, navg + 1))
            navg += 1
    elif arg == '-nosave':
        saveplot = False
    elif arg == '-noshow':
        showplot = False
    elif arg == '-tlabel':
        labelbytime = True
    elif arg == '-rbcz':
        rbcz = float(args[i+1])/rsun
    elif arg == '-tag':
        tag = '_' + args[i+1]
    elif arg == '-rvals':
        rvals_str = args[i+1].split()
        rvals = []
        for rval_str in rvals_str:
            rvals.append(float(rval_str)/rsun)

# Get plot directory and create if not already there
plotdir = dirname + '/plots/time-rad' + tag + '_B/'
if labelbytime:
    plotdir = dirname + '/plots/time-rad' + tag + '_B_tlabel/'
if not os.path.isdir(plotdir):
    os.makedirs(plotdir)

# Read in the time-latitude data (dictionary form)
print ('Getting time-radius trace from ' + datadir + the_file)
di = get_dict(datadir + the_file)

vals = di['vals']

times = di['times']
iters = di['iters']
rr = di['rr']
tt_lat = di['tt_lat']
lats_sampled = np.array(di['lats']) 
    # specific latitudes sampled for the trace
theta_inds = np.array(di['theta_inds'])
qvals = np.array(di['qvals'])

niter = di['niter']
nr = di['nr']
nlats = di['nlats']
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

if desired_lats == 'all':
    i_desiredlats = np.arange(len(theta_inds))
    lats_to_plot = lats_sampled
else:
    i_desiredlats = []
    lats_to_plot = []
    for desired_lat in desired_lats:
        i_desiredlat = np.argmin(np.abs(lats_sampled - desired_lat))
        i_desiredlats.append(i_desiredlat)
        lats_to_plot.append(lats_sampled[i_desiredlat])

# Get raw traces of br, btheta, bphi
br = vals[:, :, :, br_index]
bt = vals[:, :, :, bt_index]
bp = vals[:, :, :, bp_index]

# Normalize the time and radius
times /= time_unit
rr /= rsun

# Make meshgrid of time/radius
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
t1, t2 = times[0], times[-1] # These begin times and end times
        # will be used for labeling the plots

# set figure dimensions
subplot_width_inches = 6.5
margin_inches = 1./4.
margin_bottom_inches = 1./2. # space for x-axis label
margin_top_inches = 1./2.
margin_left_inches = 5./8. # space for latitude label
margin_right_inches = 0.9*(2 - (rbcz is None))

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

# Loop over the desired latitudes and save plots
alphabet =\
        ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm',\
        'n', 'o', 'p', 'q', 'r', 's', 't', 'u', 'v', 'w', 'x', 'y', 'z']
        # tags for latitudes so they appear in order. Don't plot more than 26 
# latitudes!

for i in range(len(lats_to_plot)):
    i_desiredlat = i_desiredlats[i]
    lat_to_plot = lats_to_plot[i] 
    br_loc = br[:, i_desiredlat, :]
    bt_loc = bt[:, i_desiredlat, :]
    bp_loc = bp[:, i_desiredlat, :]

    # Make appropriate file name to save
    if lat_to_plot >= 0.0:
        hemisphere = "N"
    else:
        hemisphere = "S"

    if labelbytime:
        savename = dirname_stripped + '_time-rad_B_' +\
                ('Prot%05.0f-to-%05.0f_' %(t1, t2)) +\
                '_latval_' + alphabet[i] +\
                ('%s%2.1f' %(hemisphere,np.abs(lat_to_plot))) + '.png'
    else:
        savename = dirname_stripped + '_time-rad_B_' +\
                str(iter1).zfill(8) + '_' + str(iter2).zfill(8) +\
                '_latval_' + alphabet[i] + '_' +\
                ('%s%2.1f' %(hemisphere,np.abs(lat_to_plot))) + '.png'

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

    plot_tl(br_loc, times, rr, fig=fig, ax=ax1, navg=navg, yvals=rvals,\
            minmax=minmax_br, units=units, xminmax=xminmax, rbcz=rbcz)
    plot_tl(bt_loc, times, rr, fig=fig, ax=ax2, navg=navg, yvals=rvals,\
            minmax=minmax_bt, units=units, xminmax=xminmax, rbcz=rbcz)
    plot_tl(bp_loc, times, rr, fig=fig, ax=ax3, navg=navg, yvals=rvals,\
            minmax=minmax_bp, units=units, xminmax=xminmax, rbcz=rbcz)

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
    # Label y-axis (radius in units of rsun)
    ax2.set_ylabel(r'$\rm{radius}\ (R_\odot)$', **csfont)

    # Put some useful information on the title
    averaging_time = (times[-1] - times[0])/niter*navg
    title = dirname_stripped + '     ' +\
            (r'${\rm{latitude}}\ =\ %0.1f\ {\rm{deg}}$' %lat_to_plot)
    if navg > 1:
        title += '     ' + (('t_avg = %.2f ' + time_label) %averaging_time)
    else:
        title += '     t_avg = none'
    ax1.set_title(title, **csfont)

    # Save the plot
    if saveplot:
        print ('Saving the time-radius plot at ' + plotdir + savename)
        plt.savefig(plotdir + savename, dpi=200)

    # Show the plot if only plotting at one latitude
    if len(lats_to_plot) == 1 and showplot:
        plt.show()
    plt.close()
