# Author: Loren Matilsky
# Date created: 03/02/2019
import matplotlib as mpl
mpl.use('TkAgg')
import matplotlib.pyplot as plt
from matplotlib import colors
plt.rcParams['mathtext.fontset'] = 'dejavuserif'
csfont = {'fontname':'DejaVu Serif'}
import numpy as np
import pickle
import sys, os
sys.path.append(os.environ['raco'])
from common import get_file_lists, get_widest_range_file, strip_dirname,\
        rsun, get_dict, get_satvals
from get_parameter import get_parameter
from plotcommon import axis_range, xy_grid, integerticks

# Get the run directory on which to perform the analysis
dirname = sys.argv[1]

# Data and plot directories
datadir = dirname + '/data/'
plotdir = dirname + '/plots/'
if (not os.path.isdir(plotdir)):
    os.makedirs(plotdir)
dirname_stripped = strip_dirname(dirname)

# Find the lvals_trace file(s) in the data directory. If there are 
# multiple, by default choose the one with widest range in the trace.
trace_lvals_file = get_widest_range_file(datadir, 'trace_lvals')

# more defaults
user_specified_minmax = False
user_specified_xminmax = False
user_specified_lmax = False
separate = False # don't separate even and odd values
nosave = False
tag = ''

desired_rvals = [0.83] # by default, plot time-radius diagram for fields 
    # mid-CZ (units of solar radius)
navg = 11 # by default average over 11 Shell_Spectra files for each time
power_type = 'mean' # can be tot, mean, conv
itype = 0
power_types_dict = {'tot': 0, 'mean': 1, 'conv': 2}

# Get command-line arguments
args = sys.argv[2:]
nargs = len(args)
for i in range(nargs):
    arg = args[i]
    if (arg == '-minmax'):
        user_specified_minmax = True
        my_min_br = float(args[i+1])
        my_max_br = float(args[i+2])
        try:
            my_min_bt = float(args[i+3])
            my_max_bt = float(args[i+4])
            my_min_bp = float(args[i+5])
            my_max_bp = float(args[i+6])
        except: # do the same values for each component
            my_min_bt = my_min_br
            my_max_bt = my_max_br
            my_min_bp = my_min_br
            my_max_bp = my_max_br
    elif (arg == '-nosave'):
        nosave = True
    elif (arg == '-separate'):
        separate = True
    elif (arg == '-usefile'):
        trace_lvals_file = args[i+1]
        trace_lvals_file = trace_lvals_file.split('/')[-1]
    elif (arg == '-rvals'):
        string_desired_rvals = args[i+1].split()
        if string_desired_rvals == ['all']:
            desired_rvals = 'all'
        else:
            desired_rvals = []
            for j in range(len(string_desired_rvals)):
                desired_rvals.append(float(string_desired_rvals[j]))
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
    elif (arg == '-lmax'):
        user_specified_lmax = True
        my_lmax = float(args[i+1])
    elif (arg == '-tag'):
        tag = '_' + args[i+1]
    elif (arg == '-ptype'):
        power_type = args[i+1]
        if not (power_type == 'mean' or power_type == 'conv'\
                or power_type == 'tot'):
            power_type = input('Please enter "tot," "mean," or "conv" for power_type:\n')
        itype = power_types_dict[power_type]

# Read in the time-lvals data (dictionary form)
print ('Reading in lvals trace from ' + datadir +\
       trace_lvals_file + ' ...')
di = get_dict(datadir + trace_lvals_file)

vals = di['vals']
times = di['times']
iters = di['iters']

lvals = di['lvals']
nell = di['nell']
lmax = di['lmax']

rvals = di['rvals']
rinds = di['rinds'] # radial locations sampled for the trace
rvals_sampled = rvals/rsun

qv = np.array(di['qv'])
lut = di['lut']

niter = di['niter']
nr = di['nr']
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

br_index = lut[801]
bt_index = lut[802]
bp_index = lut[803]

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
br = vals[:, :, :, br_index, :]
bt = vals[:, :, :, bt_index, :]
bp = vals[:, :, :, bp_index, :]

# Average these traces in time
over2 = navg//2
br_trace_all = np.zeros((niter - navg + 1, nell, nr, 3))
bt_trace_all = np.zeros((niter - navg + 1, nell, nr, 3))
bp_trace_all = np.zeros((niter - navg + 1, nell, nr, 3))

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

# If user specified lmax, change the default lvals array
il_max = nell - 1
if user_specified_lmax:
    il_max = np.argmin(np.abs(lvals - my_lmax))
    # it is easier when said and done if lvals ends on an odd value
    # (So there are equal numbers even and odd l)
#    if lvals[il_max] % 2 != 1:
#        il_max -= 1
    lvals = lvals[:il_max+1]
    br_trace_all = br_trace_all[:, :il_max+1]
    bt_trace_all = bt_trace_all[:, :il_max+1]
    bp_trace_all = bp_trace_all[:, :il_max+1]

if separate:
    # Break up quantities into l-even/l-odd
    lvals_int = lvals.astype(int)
    il_even = np.where(lvals_int % 2 == 0)[0]
    il_odd = np.where(lvals_int % 2 == 1)[0]
    lvals_even = lvals[il_even]
    lvals_odd = lvals[il_odd]

# Make grid of times and l-values
times2, lvals2 = np.meshgrid(times_trace, lvals, indexing='ij')
times2, lvals2 = xy_grid(times2, lvals2)

# Loop over the desired radii and save plots
for i in range(len(i_desiredrvals)):
    i_desiredrval = i_desiredrvals[i]
    rval_to_plot = rvals_to_plot[i]
    br_trace = br_trace_all[:, :, i_desiredrval, itype] 
    bt_trace = bt_trace_all[:, :, i_desiredrval, itype] 
    bp_trace = bp_trace_all[:, :, i_desiredrval, itype]
    
    # Make appropriate file name to save
    savename = dirname_stripped + '_time-lval_B_' + \
            str(itype).zfill(2) + '-' + power_type + '_' + \
        ('rval%0.3f_' %rval_to_plot) + str(iter1).zfill(8) + '_' +\
        str(iter2).zfill(8) + tag + '.png'

    if not user_specified_minmax:
        my_mins = []
        my_maxes = []
        cut = len(times_trace)//2

        for field in [br_trace[cut:, 1:il_max+1],\
                bt_trace[cut:, 1:il_max+1], bp_trace[cut:, 1:il_max+1]]:
            my_min, my_max = get_satvals(field+1e-22, logscale=True) 
            #+1e-22 to guard against exactly 0 values
            my_mins.append(my_min)
            my_maxes.append(my_max)
    else:
        my_mins = [my_min_br, my_min_bt, my_min_bp]
        my_maxes = [my_max_br, my_max_bt, my_max_bp]
    
    # Create figure with  3 panels in a row (time-lval plots of
    #       br, btheta, and bphi)
    fig, axs = plt.subplots(3, 1, figsize=(12, 8), sharex=True, sharey=True)
    ax1 = axs[0]; ax2 = axs[1]; ax3 = axs[2]

    # first plot: evolution of B_r spectrum
    if separate:
        times2_even, lvals2_even = np.meshgrid(times_trace, lvals_even,\
                indexing='ij')
        times2_even, lvals2_even = xy_grid(times2_even, lvals2_even)

        times2_odd, lvals2_odd = np.meshgrid(times_trace, lvals_odd,\
                indexing='ij')
        times2_odd, lvals2_odd = xy_grid(times2_odd, lvals2_odd)

        # First B_r power
        im1 = ax1.pcolormesh(times2_even, lvals2_even,\
                br_trace[:, il_even], cmap='Greys',\
                norm=colors.LogNorm(vmin=my_mins[0], vmax=my_maxes[0]))
        ax1.pcolormesh(times2_odd, -lvals2_odd, br_trace[:, il_odd],\
                cmap='Greys',\
                norm=colors.LogNorm(vmin=my_mins[0], vmax=my_maxes[0]))

        # Next, B_theta power
        im2 = ax2.pcolormesh(times2_even, lvals2_even,\
                bt_trace[:, il_even], cmap='Greys',\
                norm=colors.LogNorm(vmin=my_mins[1], vmax=my_maxes[1]))
        ax2.pcolormesh(times2_odd, -lvals2_odd, bt_trace[:, il_odd],\
                cmap='Greys',\
                norm=colors.LogNorm(vmin=my_mins[1], vmax=my_maxes[1]))

        # Finally, B_phi power
        im3 = ax3.pcolormesh(times2_even, lvals2_even,\
                bp_trace[:, il_even], cmap='Greys',\
                norm=colors.LogNorm(vmin=my_mins[2], vmax=my_maxes[2]))
        ax3.pcolormesh(times2_odd, -lvals2_odd, bp_trace[:, il_odd],\
                cmap='Greys',\
                norm=colors.LogNorm(vmin=my_mins[2], vmax=my_maxes[2]))

        # Make tick values something sensible
        tickvals, ticklabels = integerticks(np.max(lvals))
        ax1.set_yticks(tickvals)
        ax1.set_yticklabels(ticklabels)
    else:
        im1 = ax1.pcolormesh(times2, lvals2, br_trace, cmap='Greys',\
                norm=colors.LogNorm(vmin=my_mins[0], vmax=my_maxes[0]))
        im2 = ax2.pcolormesh(times2, lvals2, bt_trace, cmap='Greys',\
                norm=colors.LogNorm(vmin=my_mins[1], vmax=my_maxes[1]))
        im3 = ax3.pcolormesh(times2, lvals2, bp_trace, cmap='Greys',\
                norm=colors.LogNorm(vmin=my_mins[2], vmax=my_maxes[2]))

    # Put colorbar next to all plots (possibly normalized separately)
    # First make room and then find location of subplots
    plt.subplots_adjust(left=0.1, right=0.85, wspace=0.03, top=0.9)

    # First, B_r:
    ax_xmin, ax_xmax, ax_ymin, ax_ymax = axis_range(ax1)
    ax_delta_x = ax_xmax - ax_xmin
    ax_delta_y = ax_ymax - ax_ymin
    ax_center_x = ax_xmin + 0.5*ax_delta_x

    # Make the y label
    if separate:
        fig.text(0.5*ax_xmin, ax_ymin + 0.25*ax_delta_y, 'l odd',\
                rotation=90, va='center', ha='right')
        fig.text(0.5*ax_xmin, ax_ymin + 0.75*ax_delta_y, 'l even',\
                rotation=90, va='center', ha='right')
    else:
        ax1.set_ylabel('l-degree')

    cbar_left = ax_xmax + 0.3*(1 - ax_xmax)
    cbar_bottom = ax_ymin
    cbar_width = 0.07*(1 - ax_xmax)
    cbar_height = ax_delta_y
    cax = fig.add_axes((cbar_left, cbar_bottom, cbar_width, cbar_height))
    cax.set_title(r'$G^2$', **csfont)
    plt.colorbar(im1, cax=cax)

    # Next, B_theta:
    ax_xmin, ax_xmax, ax_ymin, ax_ymax = axis_range(ax2)
    ax_delta_x = ax_xmax - ax_xmin
    ax_delta_y = ax_ymax - ax_ymin
    ax_center_x = ax_xmin + 0.5*ax_delta_x

    # Make the y label
    if separate:
        fig.text(0.5*ax_xmin, ax_ymin + 0.25*ax_delta_y, 'l odd',\
                rotation=90, va='center', ha='right')
        fig.text(0.5*ax_xmin, ax_ymin + 0.75*ax_delta_y, 'l even',\
                rotation=90, va='center', ha='right')
    else:
        ax2.set_ylabel('l-degree')

    cbar_left = ax_xmax + 0.3*(1 - ax_xmax)
    cbar_bottom = ax_ymin
    cbar_width = 0.07*(1 - ax_xmax)
    cbar_height = ax_delta_y
    cax = fig.add_axes((cbar_left, cbar_bottom, cbar_width, cbar_height))
    cax.set_title(r'$G^2$', **csfont)
    plt.colorbar(im2, cax=cax)

    # Next, B_phi:
    ax_xmin, ax_xmax, ax_ymin, ax_ymax = axis_range(ax3)
    ax_delta_x = ax_xmax - ax_xmin
    ax_delta_y = ax_ymax - ax_ymin
    ax_center_x = ax_xmin + 0.5*ax_delta_x

    # Make the y label
    if separate:
        fig.text(0.5*ax_xmin, ax_ymin + 0.25*ax_delta_y, 'l odd',\
                rotation=90, va='center', ha='right')
        fig.text(0.5*ax_xmin, ax_ymin + 0.75*ax_delta_y, 'l even',\
                rotation=90, va='center', ha='right')
    else:
        ax3.set_ylabel('l-degree')

    cbar_left = ax_xmax + 0.3*(1 - ax_xmax)
    cbar_bottom = ax_ymin
    cbar_width = 0.07*(1 - ax_xmax)
    cbar_height = ax_delta_y
    cax = fig.add_axes((cbar_left, cbar_bottom, cbar_width, cbar_height))
    cax.set_title(r'$G^2$', **csfont)
    plt.colorbar(im3, cax=cax)

    # Label x (time) axis
    if rotation:
        timeunit = r'$P_{\rm{rot}}$'
    else:
        timeunit = 'days'

    xlabel = 'time (' + timeunit + ')'
    ax3.set_xlabel(xlabel, **csfont)

#    if user_specified_xminmax:
#        ax3.set_xlim((xmin, xmax))
#    if user_specified_lminmax:
#        ax3.set_ylim((ymin, ymax))

    # Label y-axis (radius in units of rsun)
#    ax2.set_ylabel('total power: l-value', **csfont)

    # Label the plots by B_r, B_theta, B_phi
    ax_xmin, ax_xmax = ax1.get_xlim()
    ax_Dx = ax_xmax - ax_xmin
    ax1.text(ax_xmin + 1.01*ax_Dx, 0.,  r'$B_r$')
    ax2.text(ax_xmin + 1.01*ax_Dx, 0.,  r'$B_\theta$')
    ax3.text(ax_xmin + 1.01*ax_Dx, 0.,  r'$B_\phi$')

    # Put some useful information on the title
    averaging_time = (times[-1] - times[0])/niter*navg/86400.
    title = power_type + ' l-power, Magnetic field\n' +\
            dirname_stripped + '     ' +\
            (r'$r/R_\odot\ =\ %0.3f$' %rval_to_plot) + '     ' +\
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

    if not nosave:
        # Save the plot
        print ('Saving the time-lval plot at ' + plotdir +\
            savename + ' ...')
        plt.savefig(plotdir + savename, dpi=300)

    # Show the plot if only plotting at one latitude
    if len(rvals_to_plot) == 1:
        plt.show()
