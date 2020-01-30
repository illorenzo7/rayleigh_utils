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
from rayleigh_diagnostics import Shell_Avgs
from common import strip_dirname, get_widest_range_file,\
        get_iters_from_file, get_dict, rsun, get_file_lists,\
        get_desired_range

# Get directory name and stripped_dirname for plotting purposes
dirname = sys.argv[1]
dirname_stripped = strip_dirname(dirname)

# Directory with data and plots, make the plotting directory if it doesn't
# already exist    
plotdir = dirname + '/plots/thermo_shav_vs_time/'
if not os.path.isdir(plotdir):
    os.makedirs(plotdir)
radatadir = dirname + '/Shell_Avgs/'

# Set defaults
rnorm = None
minmax = None
logscale = False
rvals = None # user can specify radii to mark by vertical lines

# Get all the file names in datadir and their integer counterparts
file_list, int_file_list, nfiles = get_file_lists(radatadir)

# Read in CLAs
args = sys.argv[2:]
nargs = len(args)

if nargs == 0:
    index_first, index_last = nfiles - 11, nfiles - 1  
    # By default average over the last 10 files
else:
    index_first, index_last = get_desired_range(int_file_list, args)

# Change defaults
for i in range(nargs):
    arg = args[i]
    if arg == '-rnorm':
        rnorm = float(args[i+1])
    elif arg == '-minmax':
        minmax = float(args[i+1]), float(args[i+2])
    elif arg == '-log':
        logscale = True
    elif arg == '-rvals':
        rvals_str = args[i+1].split()
        rvals = []
        for rval_str in rvals_str:
            rvals.append(float(rval_str))

# Read in vavg data from Shell_Avgs
sh0 = Shell_Avgs(radatadir + file_list[index_first], '')

# Grid info
rr = sh0.radius
nr = sh0.nr
ri, ro = np.min(rr), np.max(rr)
d = ro - ri
rr_height = (rr - ri)/d
rr_depth = (ro - rr)/d

# Read reference state
prs_spec_heat = 3.5e8
try:
    ref = ReferenceState(dirname + '/reference', '')
    ref_rho = ref.density
    ref_prs = ref.pressure
    ref_temp = ref.temperature
except:
    eq = equation_coefficients()
    eq.read(dirname + '/equation_coefficients')
    ref_rho = eq.functions[0]
    ref_temp = eq.functions[3]
    gam = 5.0/3.0
    gas_R = (gam - 1.0)/gam*prs_spec_heat
    ref_prs = gas_R*ref_rho*ref_temp

# User can specify what to normalize the radius by
# By default, normalize by the solar radius
if rnorm is None:
    rr_n = rr/rsun
else:
    rr_n = rr/rnorm                                           

# Plotting loop
print ('Plotting Shell_Avgs files %s through %s ...'\
       %(file_list[index_first], file_list[index_last]))

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

        # Compute maximum/minimum Reynolds numbers (ignore the upper/lower 5%
        # of the shell to avoid extreme values associated with boundaries
#        ir1, ir2 = np.argmin(np.abs(rr_depth - 0.05)),\
#                np.argmin(np.abs(rr_depth - 0.95))
#        amp_min = min(np.min(amp_vr[ir1:ir2]), np.min(amp_vt[ir1:ir2]),\
#                np.min(amp_vp[ir1:ir2]))
#        amp_max = np.max(amp_v[ir1:ir2])

#        if minmax is None:
#            if logscale:
#                ratio = amp_max/amp_min
#                ybuffer = 0.2*ratio
#                ymin = amp_min/ybuffer
#                ymax = amp_max*ybuffer
#            else:
#                difference = amp_max - amp_min
#                ybuffer = 0.2*difference
#                ymin, ymax = amp_min - ybuffer, amp_max + ybuffer
#        else:
#            ymin, ymax = minmax

#        plt.ylim((ymin, ymax))
#        if logscale:
#            plt.yscale('log')
        ymin, :
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
        plt.title(dirname_stripped + '\n' +'Velocity Amplitudes, ' +\
                  str(iter_loc).zfill(8))
        plt.legend()

        # Get ticks everywhere
        plt.minorticks_on()
        plt.tick_params(top=True, right=True, direction='in', which='both')
        plt.tight_layout()

        savefile = plotdir + str(iter_loc).zfill(8) + '.png'
        print('Saving plot at ' + savefile + ' ...')
        plt.savefig(savefile, dpi=300)
        plt.close()
