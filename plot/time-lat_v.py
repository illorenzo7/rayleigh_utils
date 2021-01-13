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
from common import *
        rsun, get_dict
from plotcommon import axis_range
from get_parameter import get_parameter

# Get the run directory on which to perform the analysis
dirname = sys.argv[1]

# Data and plot directories
datadir = dirname + '/data/'
plotdir = dirname + '/plots/'
if (not os.path.isdir(plotdir)):
    os.makedirs(plotdir)
dirname_stripped = strip_dirname(dirname)

# Find the time/latitude file(s) the data directory. If there are 
# multiple, by default choose the one with widest range in the trace.
time_latitude_file = get_widest_range_file(datadir, 'time-latitude')

# more defaults
minmax = None
xminmax = None
saveplot = True
showplot = True # will only show if plotting one figure

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
        time_latitude_file = args[i+1]
        time_latitude_file = time_latitude_file.split('/')[-1]
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

# Read in the time-latitude data (dictionary form)
print ('Getting time-latitude trace from ' + datadir +\
       time_latitude_file + ' ...')
di = get_dict(datadir + time_latitude_file)

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

vr_index = np.argmin(np.abs(qvals - 1))
vt_index = np.argmin(np.abs(qvals - 2))
vp_index = np.argmin(np.abs(qvals - 3))

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

# Get raw traces of vr, vtheta, vphi
vr = vals[:, :, :, vr_index]/100. # convert cm/s --> m/s
vt = vals[:, :, :, vt_index]/100.
vp = vals[:, :, :, vp_index]/100.

# Average these traces in time (if navg = 1, [...]_trace_av = [...]
over2 = navg//2
vr_trace_av = np.zeros((niter - navg + 1, ntheta, nrvals))
vt_trace_av = np.zeros((niter - navg + 1, ntheta, nrvals))
vp_trace_av = np.zeros((niter - navg + 1, ntheta, nrvals))
for i in range(navg):
    vr_trace_av += vr[i:niter - navg + 1 + i]
    vt_trace_av += vt[i:niter - navg + 1 + i]
    vp_trace_av += vp[i:niter - navg + 1 + i]
vr_trace_av /= navg
vt_trace_av /= navg
vp_trace_av /= navg

times_trace = times[over2:niter - over2]/Prot # time_trace is in units of
    # Prot

# Make meshgrid of time/radius
# Take into account if user specified xmin, xmax
if not xminmax is None:
    it1 = np.argmin(np.abs(times_trace - xminmax[0]))
    it2 = np.argmin(np.abs(times_trace - xminmax[1]))
    times_trace = times_trace[it1:it2+1]
    vr_trace_av = vr_trace_av[it1:it2+1]
    vt_trace_av = vt_trace_av[it1:it2+1]
    vp_trace_av = vp_trace_av[it1:it2+1]
t1, t2 = times_trace[0], times_trace[-1] # These begin times and end times
        # will be used for labeling the plots
times2, tt_lat2 = np.meshgrid(times_trace, tt_lat, indexing='ij')

# Loop over the desired radii and save plots
for i in range(len(i_desiredrvals)):
    i_desiredrval = i_desiredrvals[i]
    rval_to_plot = rvals_to_plot[i]
    vr_trace = vr_trace_av[:, :, i_desiredrval]
    vt_trace = vt_trace_av[:, :, i_desiredrval]
    vp_trace = vp_trace_av[:, :, i_desiredrval]
    
    # Make appropriate file name to save
    savename = dirname_stripped + '_time-latitude_v_' +\
            ('Prot%05.0f-to-%05.0f_' %(t1, t2)) +\
        ('rval%0.3f' %rval_to_plot) + '.png'

    if minmax is None:
        # Exclude high latitude regions for vr and and vtheta
        max_lat = 60.
        ith1, ith2 = np.argmin(np.abs(tt_lat + max_lat)),\
                np.argmin(np.abs(tt_lat - max_lat)), 
        std_vr = np.std(vr_trace[:, ith1:ith2+1])
        std_vt = np.std(vt_trace[:, ith1:ith2+1])
        std_vp = np.std(vp_trace[:, ith1:ith2+1])
        minmax_vr = -3.*std_vr, 3.*std_vr
        minmax_vt = -3.*std_vt, 3.*std_vt
        minmax_vp = -3.*std_vp, 3.*std_vp
    else:
        if len(minmax) == 2:
            minmax_vr = minmax
            minmax_vt = minmax
            minmax_vp = minmax
        elif len(minmax) == 6:
            minmax_vr = minmax[0], minmax[1]
            minmax_vt = minmax[2], minmax[3]
            minmax_vp = minmax[4], minmax[5]
     
    # Create figure with  3 panels in a row (time-latitude plots of
    #       vr, vtheta, and vphi)
    fig, axs = plt.subplots(3, 1, figsize=(12, 8), sharex=True, sharey=True)
    ax1 = axs[0]; ax2 = axs[1]; ax3 = axs[2]

    # first plot: evolution of B_r
    im1 = ax1.pcolormesh(times2, tt_lat2, vr_trace,\
            vmin=minmax_vr[0], vmax=minmax_vr[1], cmap='RdYlBu_r')
    im2 = ax2.pcolormesh(times2, tt_lat2, vt_trace,\
            vmin=minmax_vt[0], vmax=minmax_vt[1], cmap='RdYlBu_r')
    im3 = ax3.pcolormesh(times2, tt_lat2, vp_trace,\
            vmin=minmax_vp[0], vmax=minmax_vp[1], cmap='RdYlBu_r')

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
    cax.set_title('m/s', **csfont)
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
    cax.set_title('m/s', **csfont)
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
    cax.set_title('m/s', **csfont)
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
    ax1.text(ax_xmin + 1.01*ax_Dx, 0.,  r'$v_r$')
    ax2.text(ax_xmin + 1.01*ax_Dx, 0.,  r'$v_\theta$')
    ax3.text(ax_xmin + 1.01*ax_Dx, 0.,  r'$v_\phi$')

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
        print ('Saving the time-latitude plot at ' + plotdir +\
                savename + ' ...')
        plt.savefig(plotdir + savename, dpi=200)

    # Show the plot if only plotting at one latitude
    if len(rvals_to_plot) == 1 and showplot:
        plt.show()
    plt.close()
