# Author: Loren Matilsky
# Created: 11/03/2019
# This script generates the spherically averaged convective velocity
# amplitudes for many instantaneous shell slices
# (use a range specified at the command line)
# to investigate why stuff blows up
# Saves plots as
# plots/vamp_vs_time/[iter].png

# Import relevant modules
import numpy as np
import matplotlib as mpl
mpl.use('TkAgg')
import matplotlib.pyplot as plt
plt.rcParams['mathtext.fontset'] = 'dejavuserif'
csfont = {'fontname':'DejaVu Serif'}
import sys, os
sys.path.append(os.environ['rapp'])
sys.path.append(os.environ['raco'])
from rayleigh_diagnostics import Shell_Avgs
from common import strip_dirname, get_widest_range_file,\
        get_iters_from_file, get_dict, rsun, get_file_lists,\
        get_desired_range
from get_eq import get_eq
from get_parameter import get_parameter
from time_scales import compute_Prot, compute_tdt

# Get directory name and stripped_dirname for plotting purposes
dirname = sys.argv[1]
dirname_stripped = strip_dirname(dirname)

# Get the baseline time unit
rotation = get_parameter(dirname, 'rotation')
if rotation:
    time_unit = compute_Prot(dirname)
    time_label = r'$\rm{P_{rot}}$'
else:
    time_unit = compute_tdt(dirname)
    time_label = r'$\rm{TDT}$'

# Directory with data and plots, make the plotting directory if it doesn't
# already exist    
plotdir = dirname + '/plots/thermo_shav_vs_time/'
if not os.path.isdir(plotdir):
    os.makedirs(plotdir)
radatadir = dirname + '/Shell_Avgs/'

# Set defaults
rnorm = None
minmax = None
minmax_was_None = True
logscale = False
rvals = None # user can specify radii to mark by vertical lines
nrec = 1 # by default only plot 1 record from each Shell_Avgs file
nskip = 1 # by default don't skip any Shell_Avgs files in the range
    # for nskip = 3, only read every third Shell_Avgs file, etc.
ntot = None # user can specify a total number of plots they want to see
    # in the desired range

# Get all the file names in datadir and their integer counterparts
file_list, int_file_list, nfiles = get_file_lists(radatadir)

# Read in CLAs
args = sys.argv[2:]
nargs = len(args)

index_first, index_last = get_desired_range(int_file_list, args)

# Change defaults
for i in range(nargs):
    arg = args[i]
    if arg == '-rnorm':
        rnorm = float(args[i+1])
    elif arg == '-minmax':
        minmax = float(args[i+1]), float(args[i+2])
        minmax_was_None = False
    elif arg == '-log':
        logscale = True
    elif arg == '-rvals':
        rvals_str = args[i+1].split()
        rvals = []
        for rval_str in rvals_str:
            rvals.append(float(rval_str))
    elif arg == '-nrec':
        nrec = int(args[i+1])
    elif arg == '-nskip':
        nskip = int(args[i+1])
    elif arg == '-ntot':
        ntot = int(args[i+1])

if not ntot is None: # This overrides nskip even if user specified it
    nskip = (index_last - index_first)//ntot

# Read in vavg data from Shell_Avgs
sh0 = Shell_Avgs(radatadir + file_list[index_first], '')

# Grid info
rr = sh0.radius
nr = sh0.nr
ri, ro = np.min(rr), np.max(rr)
shell_depth = ro - ri
d = ro - ri
rr_height = (rr - ri)/d
rr_depth = (ro - rr)/d

# Read reference state
eq = get_eq(dirname)
prs_spec_heat = 3.5e8
ref_rho = eq.density
ref_prs = eq.pressure
ref_temp = eq.temperature

# User can specify what to normalize the radius by
# By default, normalize by the solar radius
if rnorm is None:
    rr_n = rr/rsun
else:
    rr_n = rr/rnorm                                           

# Plotting loop
print ('Plotting Shell_Avgs files %s through %s ...'\
       %(file_list[index_first], file_list[index_last]))

for i in range(index_first, index_last + 1, nskip):
    if i == index_first:
        sh = sh0
    else:   
        sh = Shell_Avgs(radatadir + file_list[i], '')

    vals = sh.vals
    lut = sh.lut

    nrec_tot = sh.niter
    nstep = nrec_tot // nrec
    if nstep == 0:
        nstep = 1
    for j in range(0, nrec_tot, nstep):
        iter_loc = sh.iters[j]
        t_loc = sh.time[j]

        # Thermal deviations
        entropy = vals[:, 0, lut[501], j]/prs_spec_heat
        prs = vals[:, 0, lut[502], j]/ref_prs
        # Calculate temp. from EOS
        poly_n = 1.5
        temp = prs/(poly_n + 1.) + entropy

        # Calculate  density from Ideal Gas Law
        rho = prs - temp

        # Create the plot
        fig = plt.figure()
        ax = fig.add_subplot(111)

        # Plot thermal deviations vs radius
        ax.plot(rr_n, entropy,\
                label=r'$\langle S\rangle_{\rm{sph}}/c_{\rm{p}}$')
        ax.plot(rr_n, prs,\
                label=r'$\langle P\rangle_{\rm{sph}}/\overline{P}$')
        ax.plot(rr_n, temp,\
                label=r'$\langle T\rangle_{\rm{sph}}/\overline{T}$')
        ax.plot(rr_n, rho,\
                label=r'$\langle \rho\rangle_{\rm{sph}}/\overline{\rho}$')

        # Label the axes
        if rnorm is None:
            plt.xlabel(r'$r/R_\odot$',fontsize=12, **csfont)
        else:
            plt.xlabel(r'r/(%.1e cm)' %rnorm, fontsize=12, **csfont)

        plt.ylabel('thermal deviation',fontsize=12,\
                **csfont)

        # Set the axis limits
        xmin, xmax = np.min(rr_n), np.max(rr_n)
        plt.xlim((xmin, xmax))
        if minmax_was_None:
            # Compute maximum/minimum thermo. pert. (ignore the upper/lower 5%
            # of the shell to avoid extreme values associated with boundary conditions
            rr_depth = (ro - rr)/shell_depth
            mindepth, maxdepth = 0.0, 1.0
            ir1, ir2 = np.argmin(np.abs(rr_depth - mindepth)),\
                    np.argmin(np.abs(rr_depth - maxdepth))

            mmin = min(np.min(entropy[ir1:ir2+1]), np.min(prs[ir1:ir2+1]),\
                    np.min(temp[ir1:ir2+1]))
            mmax = np.max(rho[ir1:ir2+1])
            difference = mmax - mmin
            ybuffer = 0.2*difference
            ymin, ymax = mmin - ybuffer, mmax + ybuffer
        else:
            ymin, ymax = minmax

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

        # Create a title    
        if rotation:
            time_string = ('t = %.1f ' %(t_loc/time_unit)) + time_label +\
                    ' (1 ' + time_label + (' = %.2f days)'\
                    %(time_unit/86400.))
        else:
            time_string = ('t = %.3f ' %(t_loc/time_unit)) + time_label +\
                    ' (1 ' + time_label + (' = %.1f days)'\
            %(time_unit/86400.))

        plt.title(dirname_stripped + '\n' + time_string)
        plt.legend()

        # Get ticks everywhere
        plt.minorticks_on()
        plt.tick_params(top=True, right=True, direction='in', which='both')
        plt.tight_layout()

        savefile = plotdir + str(iter_loc).zfill(8) + '.png'
        print('Saving plot at ' + savefile)
        plt.savefig(savefile, dpi=300)
        plt.close()
