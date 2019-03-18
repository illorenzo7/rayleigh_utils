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
sys.path.append(os.environ['co'])
from common import get_file_lists, get_widest_range_file, strip_dirname,\
        rsun
from get_parameter import get_parameter

def axis_range(ax): # gets subplot coordinates on a figure in "normalized"
        # coordinates
    pos = plt.get(ax, 'position')
    bottom_left = pos.p0
    top_right = pos.p1
    xmin, xmax = bottom_left[0], top_right[0]
    ymin, ymax = bottom_left[1], top_right[1]
    
    return xmin, xmax, ymin, ymax

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
time_radius_file = get_widest_range_file(datadir, 'time-radius')

# more defaults
user_specified_minmax = False
user_specified_xminmax = False
tag = ''

desired_lats = [0.] # by default, plot time-radius diagram for fields 
    # at the equator
navg = 11 # by default average over 11 AZ_Avgs files for each time

# Get command-line arguments
args = sys.argv[2:]
nargs = len(args)
for i in range(nargs):
    arg = args[i]
    if (arg == '-minmax'):
        user_specified_minmax = True
        my_min = float(args[i+1])
        my_max = float(args[i+2])
    elif (arg == '-usefile'):
        time_radius_file = args[i+1]
        time_radius_file = time_radius_file.split('/')[-1]
    elif (arg == '-lats'):
        string_lats = args[i+1].split()
        if string_lats == ['all']:
            desired_lats = 'all'
        else:
            desired_lats = []
            for j in range(len(string_lats)):
                desired_lats.append(float(string_lats[j]))
    elif (arg == '-navg'):
        navg = int(args[i+1])
        if (navg % 2 == 0):
            print ("Please don't enter even values for navg!")
            print ("Replacing navg = %i with navg = %i" %(navg, navg + 1))
            navg += 1
    elif (arg == '-xminmax'):
        user_specified_xminmax = True
        xmin = float(args[i+1])
        xmax = float(args[i+2])
    elif (arg == '-tag'):
        tag = '_' + args[i+1]

# Read in the time-latitude data (dictionary form)
print ('Getting time-radius trace from ' + datadir +\
       time_radius_file + ' ...')
try:
    di = np.load(datadir + time_radius_file, encoding='latin1').item()
except:
    f = open(datadir + time_radius_file, 'rb')
    di = pickle.load(f)
    f.close()

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

# Get global rotation rate, if present
rotation = get_parameter(dirname, 'rotation')
if rotation:
    angular_velocity = get_parameter(dirname, 'angular_velocity')
    Prot = 2*np.pi/angular_velocity
    tnorm = Prot # normalize the time axis by rotation period if applicable
else:
    tnorm = 86400. # normalize by "days"

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

# Average these traces in time
over2 = navg//2
br_trace_all = np.zeros((niter - navg + 1, nlats, nr))
bt_trace_all = np.zeros((niter - navg + 1, nlats, nr))
bp_trace_all = np.zeros((niter - navg + 1, nlats, nr))
for i in range(navg):
    br_trace_all += br[i:niter - navg + 1 + i]
    bt_trace_all += bt[i:niter - navg + 1 + i]
    bp_trace_all += bp[i:niter - navg + 1 + i]
br_trace_all /= navg
bt_trace_all /= navg
bp_trace_all /= navg

times_trace = times[over2:niter - over2]/tnorm
# Make meshgrid of time/radius
# Take into account if user specified xmin, xmax
if user_specified_xminmax:
    it1 = np.argmin(np.abs(times_trace - xmin))
    it2 = np.argmin(np.abs(times_trace - xmax))
    times_trace = times_trace[it1:it2+1]
    br_trace_all = br_trace_all[it1:it2+1]
    bt_trace_all = bt_trace_all[it1:it2+1]
    bp_trace_all = bp_trace_all[it1:it2+1]

times2, rr2 = np.meshgrid(times_trace, rr/rsun, indexing='ij')

# Loop over the desired latitudes and save plots
for i in range(len(lats_to_plot)):
    i_desiredlat = i_desiredlats[i]
    lat_to_plot = lats_to_plot[i]
    br_trace = br_trace_all[:, i_desiredlat, :]
    bt_trace = bt_trace_all[:, i_desiredlat, :]
    bp_trace = bp_trace_all[:, i_desiredlat, :]
    
    # Make appropriate file name to save
    savename = dirname_stripped + '_time-radius_B_' +\
        ('lat%+2.1f_' %lat_to_plot) + str(iter1).zfill(8) + '_' +\
        str(iter2).zfill(8) + tag + '.png'

    if not user_specified_minmax:
        std_br = np.std(br_trace)
        std_bt = np.std(bt_trace)
        std_bp = np.std(bp_trace)
        my_min_br, my_max_br = -3.*std_br, 3.*std_br
        my_min_bt, my_max_bt = -3.*std_bt, 3.*std_bt
        my_min_bp, my_max_bp = -3.*std_bp, 3.*std_bp
    else:
        my_min_br, my_max_br = my_min, my_max
        my_min_bt, my_max_bt = my_min, my_max
        my_min_bp, my_max_bp = my_min, my_max
            
    # Create figure with  3 panels in a row (time-latitude plots of
    #       br, btheta, and bphi)
    fig, axs = plt.subplots(3, 1, figsize=(12, 8), sharex=True, sharey=True)
    ax1 = axs[0]; ax2 = axs[1]; ax3 = axs[2]

    # first plot: evolution of B_r
    im1 = ax1.pcolormesh(times2, rr2, br_trace[:, :],\
            vmin=my_min_br, vmax=my_max_br, cmap='RdYlBu_r')
    im2 = ax2.pcolormesh(times2, rr2, bt_trace[:, :],\
            vmin=my_min_bt, vmax=my_max_bt, cmap='RdYlBu_r')
    im3 = ax3.pcolormesh(times2, rr2, bp_trace[:, :],\
            vmin=my_min_bp, vmax=my_max_bp, cmap='RdYlBu_r')

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
    cax.set_title('Gauss', **csfont)
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
    cax.set_title('Gauss', **csfont)
    plt.colorbar(im2, cax=cax)

    # Next, B_phi:
    ax_xmin, ax_xmax, ax_ymin, ax_ymax = axis_range(ax3)
    ax_delta_x = ax_xmax - ax_xmin
    ax_delta_y = ax_ymax - ax_ymin
    ax_center_x = ax_xmin + 0.5*ax_delta_x

    cbar_left = ax_xmax + 0.3*(1 - ax_xmax)
    cbar_bottom = ax_ymin
    cbar_width = 0.07*(1 - ax_xmax)
    cbar_height = ax_delta_y
    cax = fig.add_axes((cbar_left, cbar_bottom, cbar_width, cbar_height))
    cax.set_title('Gauss', **csfont)
    plt.colorbar(im3, cax=cax)

    # Label x (time) axis
    if rotation:
        timeunit = r'$P_{\rm{rot}}$'
    else:
        timeunit = 'days'

    xlabel = 'time (' + timeunit + ')'
    ax3.set_xlabel(xlabel, **csfont)

    # Label y-axis (radius in units of rsun)
    ax2.set_ylabel(r'$\rm{radius}\ (R_\odot)$', **csfont)


    # Label the plots by B_r, B_theta, B_phi
    min_r = np.min(rr)/rsun
    shell_depth = (np.max(rr) - np.min(rr))/rsun
    mid_rval = min_r + 0.5*shell_depth
    ax_xmin, ax_xmax = ax1.get_xlim()
    ax_Dx = ax_xmax - ax_xmin
    ax1.text(ax_xmin + 1.01*ax_Dx, mid_rval,  r'$B_r$')
    ax2.text(ax_xmin + 1.01*ax_Dx, mid_rval,  r'$B_\theta$')
    ax3.text(ax_xmin + 1.01*ax_Dx, mid_rval,  r'$B_\phi$')

    # Put some useful information on the title
    averaging_time = (times[-1] - times[0])/niter*navg/86400.
    title = dirname_stripped + '     ' +\
            ('lat = %+2.1f deg' %lat_to_plot) + '     ' +\
            ('t_avg = %.1f days' %averaging_time) +\
            '     ' + ('time_unit = %.1f days' %(tnorm/86400.))

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
    print ('Saving the time-latitude plot at ' + plotdir +\
            savename + ' ...')
    plt.savefig(plotdir + savename, dpi=300)

    # Show the plot if only plotting at one latitude
    if len(lats_to_plot) == 1:
        plt.show()
