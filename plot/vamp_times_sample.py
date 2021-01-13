# Author: Loren Matilsky
# Created: 11/03/2019
# This script generates the spherically averaged convective velocity
# amplitudes for many instantaneous shell slices
# (use a range specified at the command line)
# to investigate why stuff blows up
# Saves plots as
# plots/vamp_times_sample/[iter].png

# Import relevant modules
import numpy as np
import matplotlib as mpl
mpl.use('TkAgg')
import matplotlib.pyplot as plt
plt.rcParams['mathtext.fontset'] = 'dejavuserif'
csfont = {'fontname':'DejaVu Serif'}
import sys, os
sys.path.append(os.environ['raco'])
sys.path.append(os.environ['rapp'])
from rayleigh_diagnostics import Shell_Avgs
from common import *

# Get directory name and stripped_dirname for plotting purposes
dirname = sys.argv[1]
dirname_stripped = strip_dirname(dirname)

# domain bounds
ncheby, domain_bounds = get_domain_bounds(dirname)
ri = np.min(domain_bounds)
ro = np.max(domain_bounds)
d = ro - ri

# Directory with Shell_Avgs
radatadir = dirname + '/Shell_Avgs/'

# Set defaults
rnorm = None
minmax = None
logscale = False
rvals = [] # user can specify radii to mark by vertical lines

# Get all the file names in datadir and their integer counterparts
file_list, int_file_list, nfiles = get_file_lists(radatadir)

# Read in CLAs
args = sys.argv[2:]
nargs = len(args)

the_tuple = get_desired_range(int_file_list, args)
if the_tuple is None: # By default average over the last 10 files
    index_first, index_last = nfiles - 11, nfiles - 1  
else:
    index_first, index_last = the_tuple

# Change defaults
for i in range(nargs):
    arg = args[i]
    if arg == '-rnorm':
        rnorm = float(args[i+1])
    elif arg == '-minmax':
        minmax = float(args[i+1]), float(args[i+2])
    elif arg == '-log':
        logscale = True
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

# Get the baseline time unit
rotation = get_parameter(dirname, 'rotation')
if rotation:
    time_unit = compute_Prot(dirname)
    time_label = r'$\rm{P_{rot}}$'
else:
    time_unit = compute_tdt(dirname)
    time_label = r'$\rm{TDT}$'

# Make plot directory
if logscale:
    plotdir = dirname + '/plots/vamp_times_sample_log/'
else:
    plotdir = dirname + '/plots/vamp_times_sample/'
if not os.path.isdir(plotdir):
    os.makedirs(plotdir)

# Read in vavg data from Shell_Avgs
sh0 = Shell_Avgs(radatadir + file_list[index_first], '')

# Grid info
rr = sh0.radius
nr = sh0.nr
ri, ro = np.min(rr), np.max(rr)
d = ro - ri
rr_height = (rr - ri)/d
rr_depth = (ro - rr)/d

# User can specify what to normalize the radius by
# By default, normalize by the solar radius
if rnorm is None:
    rr_n = rr/rsun
else:
    rr_n = rr/rnorm                                           

# Plotting loop
print ('Plotting Shell_Avgs files %s through %s'\
       %(file_list[index_first], file_list[index_last]))
print('Saving plots to '  + plotdir)
for i in range(index_first, index_last + 1):
    if i == index_first:
        sh = sh0
    else:   
        sh = Shell_Avgs(radatadir + file_list[i], '')

    local_ntimes = sh.niter
    vals = sh.vals
    lut = sh.lut
    for j in range(local_ntimes):
        iter_loc = sh.iters[j]
        t_loc = sh.time[j]

        # Make the savename like for Mollweide times sample
        savename = 'vamp_iter' + str(iter_loc).zfill(8) + '.png'
        print('Plotting: ' + savename)

        # Convective velocity amplitudes...
        vsq_r, vsq_t, vsq_p = vals[:, 0, lut[422], j],\
                vals[:, 0, lut[423], j],\
                vals[:, 0, lut[424], j] 
        vsq = vsq_r + vsq_t + vsq_p

        amp_v = np.sqrt(vsq)/100.
        amp_vr = np.sqrt(vsq_r)/100.
        amp_vt = np.sqrt(vsq_t)/100.
        amp_vp = np.sqrt(vsq_p)/100.

        # Create the plot
        fig = plt.figure()
        ax = fig.add_subplot(111)

        # Get extrema values for diff. rot.
        maxes = [] # Get the max-value of Omega for plotting purposes
        mins = []  # ditto for the min-value
                                                       

        # Plot Re vs radius
        ax.plot(rr_n, amp_v, label= r'$(v^\prime)_{\rm{rms}}$')
        ax.plot(rr_n, amp_vr, label= r'$(v^\prime_r)_{\rm{rms}}$')
        ax.plot(rr_n, amp_vt, label = r'$(v^\prime_\theta)_{\rm{rms}}$')
        ax.plot(rr_n, amp_vp, label = r'$(v^\prime_\phi)_{\rm{rms}}$')

        # Label the axes
        if rnorm is None:
            plt.xlabel(r'$r/R_\odot$',fontsize=12, **csfont)
        else:
            plt.xlabel(r'r/(%.1e cm)' %rnorm, fontsize=12, **csfont)

        plt.ylabel('velocity (m/s)',fontsize=12,\
                **csfont)

        # Set the axis limits
        xmin, xmax = np.min(rr_n), np.max(rr_n)
        plt.xlim((xmin, xmax))

        # Compute maximum/minimum Reynolds numbers 
        # (ignore the upper/lower 5%
        # of the shell to avoid extreme values associated with boundaries
        ir1, ir2 = np.argmin(np.abs(rr_depth - 0.05)),\
                np.argmin(np.abs(rr_depth - 0.95))
        amp_min = min(np.min(amp_vr[ir1:ir2]), np.min(amp_vt[ir1:ir2]),\
                np.min(amp_vp[ir1:ir2]))
        amp_max = np.max(amp_v[ir1:ir2])

        if minmax is None:
            if logscale:
                ratio = amp_max/amp_min
                ybuffer = ratio**0.2
                ymin = amp_min/ybuffer
                ymax = amp_max*ybuffer
            else:
                difference = amp_max - amp_min
                ybuffer = 0.2*difference
                ymin, ymax = 0., amp_max + ybuffer
        else:
            ymin, ymax = minmax

        if logscale:
            plt.yscale('log')

        plt.ylim((ymin, ymax))

        xvals = np.linspace(xmin, xmax, 100)
        yvals = np.linspace(ymin, ymax, 100)

        # Mark radii if desired
        if not rvals is None:
            for rval in rvals:
                if rnorm is None:
                    rval_n = rval/rsun
                else:
                    rval_n = rval/rnorm
                plt.plot(rval_n + np.zeros(100), yvals, 'k--')

        # Make title
        if rotation:
            time_string = ('t = %.1f ' %(t_loc/time_unit)) + time_label +\
                    ' (1 ' + time_label + (' = %.2f days)'\
                    %(time_unit/86400.))
        else:
            time_string = ('t = %.3f ' %(t_loc/time_unit)) + time_label +\
                    ' (1 ' + time_label + (' = %.1f days)'\
                    %(time_unit/86400.))

        # Create a title    
        title = dirname_stripped + '\n' +'Velocity Amplitudes' +\
                '     ' + time_string

        plt.title(title, **csfont)
        plt.legend()

        # Get ticks everywhere
        plt.minorticks_on()
        plt.tick_params(top=True, right=True, direction='in', which='both')
        plt.tight_layout()

        savefile = plotdir + savename
        plt.savefig(savefile, dpi=300)
        plt.close()
