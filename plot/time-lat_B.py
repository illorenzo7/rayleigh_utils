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
from plotcommon import axis_range
from get_parameter import get_parameter

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
saveplot = True
showplot = True # will only show if plotting one figure
labelbytime = False # by default label by first/last iteration number
# not first/last time

desired_rvals = [0.83] # by default, plot time-radius diagram for fields 
    # mid-CZ (units of solar radius)
navg = 1 # by default average over 1 AZ_Avgs instance (no average)
# for navg > 1, a "sliding average" will be used.

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
        string_desired_rvals = args[i+1].split()
        if string_desired_rvals == ['all']:
            desired_rvals = 'all'
        else:
            desired_rvals = []
            for j in range(len(string_desired_rvals)):
                desired_rvals.append(float(string_desired_rvals[j]))
    elif arg == '-navg':
        navg = int(args[i+1])
        if navg % 2 == 0:
            print ("Please don't enter even values for navg!")
            print ("Replacing navg = %i with navg = %i" %(navg, navg + 1))
            navg += 1
    elif arg == '-xminmax':
        xminmax = float(args[i+1]), float(args[i+2])
    elif arg == '-nosave':
        saveplot = False
    elif arg == '-noshow':
        showplot = False
    elif arg == '-tlabel':
        labelbytime = True

# Get plot directory and create if not already there
plotdir = dirname + '/plots/time-lat/'
if labelbytime:
    plotdir = dirname + '/plots/time-lat_tlabel/'
if (not os.path.isdir(plotdir)):
    os.makedirs(plotdir)

# Read in the time-latitude data (dictionary form)
print ('Getting time-latitude trace from ' + datadir + the_file)
di = get_dict(datadir + the_file)

vals = di['vals']

times = di['times']
iters = di['iters']
rr = di['rr']
ri = di['ri']; ro = di['ro']; shell_depth = ro - ri
tt_lat = di['tt_lat']
rinds = di['rinds'] # radial locations sampled for the trace
ntheta = di['ntheta']
rvals_sampled = rr[rinds]/rsun

qvals = np.array(di['qvals'])

niter = di['niter']
nr = di['nr']
nrvals = di['ndepths']
nq = di['nq']

iter1 = di['iter1']
iter2 = di['iter2']

# Get global rotation rate; this script fails for non-rotating models
angular_velocity = get_parameter(dirname, 'angular_velocity')
Prot = 2*np.pi/angular_velocity

br_index = np.argmin(np.abs(qvals - 801))
bt_index = np.argmin(np.abs(qvals - 802))
bp_index = np.argmin(np.abs(qvals - 803))

i_desiredrvals = []
rvals_to_plot = []
if desired_rvals == 'all':
    i_desiredrvals = np.arange(len(rinds))
    rvals_to_plot = rvals_sampled
else:
    i_desiredrvals = []
    rvals_to_plot = []
    for desired_rval in desired_rvals:
        i_desiredrval = np.argmin(np.abs(rvals_sampled - desired_rval))
        i_desiredrvals.append(i_desiredrval)
        rvals_to_plot.append(rvals_sampled[i_desiredrval])

# Get raw traces of br, btheta, bphi
br = vals[:, :, :, br_index]
bt = vals[:, :, :, bt_index]
bp = vals[:, :, :, bp_index]

# Average these traces in time (if navg = 1, [...]_trace_av = [...]
over2 = navg//2
br_trace_av = np.zeros((niter - navg + 1, ntheta, nrvals))
bt_trace_av = np.zeros((niter - navg + 1, ntheta, nrvals))
bp_trace_av = np.zeros((niter - navg + 1, ntheta, nrvals))
for i in range(navg):
    br_trace_av += br[i:niter - navg + 1 + i]
    bt_trace_av += bt[i:niter - navg + 1 + i]
    bp_trace_av += bp[i:niter - navg + 1 + i]
br_trace_av /= navg
bt_trace_av /= navg
bp_trace_av /= navg

times_trace = times[over2:niter - over2]/Prot # time_trace is in units of
    # Prot

# Make meshgrid of time/latitude
# Take into account if user specified xmin, xmax
if not xminmax is None:
    it1 = np.argmin(np.abs(times_trace - xminmax[0]))
    it2 = np.argmin(np.abs(times_trace - xminmax[1]))
    times_trace = times_trace[it1:it2+1]
    br_trace_av = br_trace_av[it1:it2+1]
    bt_trace_av = bt_trace_av[it1:it2+1]
    bp_trace_av = bp_trace_av[it1:it2+1]
t1, t2 = times_trace[0], times_trace[-1] # These begin times and end times
        # will be used for labeling the plots
times2, tt_lat2 = np.meshgrid(times_trace, tt_lat, indexing='ij')

# Loop over the desired radii and save plots
for i in range(len(i_desiredrvals)):
    i_desiredrval = i_desiredrvals[i]
    rval_to_plot = rvals_to_plot[i]
    br_trace = br_trace_av[:, :, i_desiredrval]
    bt_trace = bt_trace_av[:, :, i_desiredrval]
    bp_trace = bp_trace_av[:, :, i_desiredrval]
    
    # Make appropriate file name to save
    if labelbytime:
        savename = dirname_stripped + '_time-lat_B_' +\
                ('Prot%05.0f-to-%05.0f_' %(t1, t2)) +\
            ('rval%0.3f' %rval_to_plot) + '.png'
    else:
        savename = dirname_stripped + '_time-lat_B_' +\
                ('%08i_%08i_' %(iter1, iter2)) +\
            ('rval%0.3f' %rval_to_plot) + '.png'

    if minmax is None:
        std_br = np.std(br_trace)
        std_bt = np.std(bt_trace)
        std_bp = np.std(bp_trace)
        minmax_br = -3.*std_br, 3.*std_br
        minmax_bt = -3.*std_bt, 3.*std_bt
        minmax_bp = -3.*std_bp, 3.*std_bp
    else:
        if len(minmax) == 2:
            minmax_br = minmax
            minmax_bt = minmax
            minmax_bp = minmax
        elif len(minmax) == 6:
            minmax_br = minmax[0], minmax[1]
            minmax_bt = minmax[2], minmax[3]
            minmax_bp = minmax[4], minmax[5]
     
    # Create figure with  3 panels in a row (time-latitude plots of
    #       br, btheta, and bphi)
    fig, axs = plt.subplots(3, 1, figsize=(12, 8), sharex=True, sharey=True)
    ax1 = axs[0]; ax2 = axs[1]; ax3 = axs[2]

    # first plot: evolution of B_r
    im1 = ax1.pcolormesh(times2, tt_lat2, br_trace,\
            vmin=minmax_br[0], vmax=minmax_br[1], cmap='RdYlBu_r')
    im2 = ax2.pcolormesh(times2, tt_lat2, bt_trace,\
            vmin=minmax_bt[0], vmax=minmax_bt[1], cmap='RdYlBu_r')
    im3 = ax3.pcolormesh(times2, tt_lat2, bp_trace,\
            vmin=minmax_bp[0], vmax=minmax_bp[1], cmap='RdYlBu_r')

    # Put colorbar next to all plots (possibly normalized separately)
    # First make room and then find location of subplots
    plt.subplots_adjust(left=0.1, right=0.85, wspace=0.03, top=0.9)

    # First, B_r:
    ax_xmin, ax_xmax, ax_ymin, ax_ymax = axis_range(ax1)
    ax_delta_x = ax_xmax - ax_xmin
    ax_delta_y = ax_ymax - ax_ymin
    ax_center_x = ax_xmin + 0.5*ax_delta_x

    cbar_left = ax_xmax + 0.3*(1 - ax_xmax)
    cbar_bottom = ax_ymin
    cbar_width = 0.07*(1 - ax_xmax)
    cbar_height = ax_delta_y
    cax = fig.add_axes((cbar_left, cbar_bottom, cbar_width, cbar_height))
    cax.set_title('G', **csfont)
    plt.colorbar(im1, cax=cax)

    # Next, B_theta:
    ax_xmin, ax_xmax, ax_ymin, ax_ymax = axis_range(ax2)
    ax_delta_x = ax_xmax - ax_xmin
    ax_delta_y = ax_ymax - ax_ymin
    ax_center_x = ax_xmin + 0.5*ax_delta_x

    cbar_left = ax_xmax + 0.3*(1 - ax_xmax)
    cbar_bottom = ax_ymin
    cbar_width = 0.07*(1 - ax_xmax)
    cbar_height = ax_delta_y
    cax = fig.add_axes((cbar_left, cbar_bottom, cbar_width, cbar_height))
    cax.set_title('G', **csfont)
    plt.colorbar(im2, cax=cax)

    # Finally, B_phi:
    ax_xmin, ax_xmax, ax_ymin, ax_ymax = axis_range(ax3)
    ax_delta_x = ax_xmax - ax_xmin
    ax_delta_y = ax_ymax - ax_ymin
    ax_center_x = ax_xmin + 0.5*ax_delta_x

    cbar_left = ax_xmax + 0.3*(1 - ax_xmax)
    cbar_bottom = ax_ymin
    cbar_width = 0.07*(1 - ax_xmax)
    cbar_height = ax_delta_y
    cax = fig.add_axes((cbar_left, cbar_bottom, cbar_width, cbar_height))
    cax.set_title('G', **csfont)
    plt.colorbar(im3, cax=cax)

    # Label x (time) axis
    timeunit = r'$P_{\rm{rot}}$'
    xlabel = 'time (' + timeunit + ')'
    ax3.set_xlabel(xlabel, **csfont)

    ax3.set_xlim((t1, t2))

    # Label y-axis (radius in units of rsun)
    ax2.set_ylabel('latitude (deg.)', **csfont)
    ax2.set_yticks(np.arange(-90, 90, 30))

    # Label the plots by B_r, B_theta, B_phi
    ax_xmin, ax_xmax = ax1.get_xlim()
    ax_Dx = ax_xmax - ax_xmin
    ax1.text(ax_xmin + 1.01*ax_Dx, 0.,  r'$B_r$')
    ax2.text(ax_xmin + 1.01*ax_Dx, 0.,  r'$B_\theta$')
    ax3.text(ax_xmin + 1.01*ax_Dx, 0.,  r'$B_\phi$')

    # Put some useful information on the title
    averaging_time = (times[-1] - times[0])/niter*navg/Prot
    title = dirname_stripped + '     ' +\
            (r'$r/R_\odot\ =\ %0.3f$' %rval_to_plot) +\
            '     ' + ('P_rot = %.1f days' %(Prot/86400.))
    if navg > 1:
        title += '     ' + ('t_avg = %.1f Prot' %averaging_time)
    else:
        title += '     t_avg = none'

    ax1.set_title(title, **csfont)
    # Get ticks everywhere
    plt.sca(ax1)
    plt.minorticks_on()
    plt.tick_params(top=True, right=True, direction='in', which='both')

    plt.sca(ax2)
    plt.minorticks_on()
    plt.tick_params(top=True, right=True, direction='in', which='both')

    plt.sca(ax3)
    plt.minorticks_on()
    plt.tick_params(top=True, right=True, direction='in', which='both')

    # Save the plot
    if saveplot:
        print ('Saving the time-latitude plot at ' + plotdir + savename)
        plt.savefig(plotdir + savename, dpi=200)

    # Show the plot if only plotting at one latitude
    if len(rvals_to_plot) == 1 and showplot:
        plt.show()
    plt.close()
