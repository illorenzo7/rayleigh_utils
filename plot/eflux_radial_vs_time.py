###############################################
# Author: Loren Matilsky
# Date created: 05/06/2020
#
# This script plots the radial energy fluxes as functions of
# For more plots, specify -nrec

import matplotlib as mpl
mpl.use('TkAgg')
import matplotlib.pyplot as plt
plt.rcParams['mathtext.fontset'] = 'dejavuserif'
csfont = {'fontname':'DejaVu Serif'}
import numpy as np
import sys, os
sys.path.append(os.environ['rapp'])
sys.path.append(os.environ['raco'])
from rayleigh_diagnostics import Shell_Avgs
from common import *
from compute_grid_info import compute_theta_grid

# Get the run directory on which to perform the analysis
dirname = sys.argv[1]
dirname_stripped = strip_dirname(dirname)

# domain bounds
ncheby, domain_bounds = get_domain_bounds(dirname)
ri = np.min(domain_bounds)
ro = np.max(domain_bounds)
d = ro - ri

# Directory with data and plots, make the plotting directory if it doesn't
# already exist    
datadir = dirname + '/data/'
plotdir = dirname + '/plots/eflux_radial_vs_time/'
if (not os.path.isdir(plotdir)):
    os.makedirs(plotdir)

# Get list of shell slices at all possible times to plot
radatadir = dirname + '/Shell_Avgs/'
file_list, int_file_list, nfiles = get_file_lists(radatadir)

# Get command-line arguments to adjust the interval of averaging files
minmax = None
xminmax = None
minmax_was_None = True
xminmax_was_None = True
rnorm = None
rvals = []
plot_enth_fluc = False
mark_bcz = False
lw = 1.0 # regular width lines
dpi = 300.
rmaxwindow = None # If present, set the bounds to zoom in on the boundaries
        # with a window of a particular size
rminwindow = None

args = sys.argv[2:]
nargs = len(args)

# get the time range to make plots
the_tuple = get_desired_range(int_file_list, args)
if the_tuple is None: # By default plot the last 10 Shell_Slices
    index_first, index_last = nfiles - 11, nfiles - 1  
else:
    index_first, index_last = the_tuple

nrec = 1 # by default only plot 1 record from each Shell_Avgs file
nskip = 1 # by default don't skip any Shell_Avgs files in the range
    # for nskip = 3, only read every third Shell_Avgs file, etc.
ntot = None # user can specify a total number of plots they want to see
    # in the desired range

for i in range(nargs):
    arg = args[i]
    if arg == '-minmax':
        minmax = float(args[i+1]), float(args[i+2])
        minmax_was_None = False
    elif arg == '-xminmax':
        xminmax = float(args[i+1]), float(args[i+2])
        xminmax_was_None = False
    elif arg == '-rnorm':
        rnorm = float(args[i+1])
    elif arg == '-fluc':
        plot_enth_fluc = True
    elif arg == '-depths':
        strings = args[i+1].split()
        for st in strings:
            rval = ro - float(st)*d
            rvals.append(rval)
    elif arg == '-depthscz':
        rm = domain_bounds[1]
        dcz = ro - rm
        strings = args[i+1].split()
        for st in strings:
            rval = ro - float(st)*dcz
            rvals.append(rval)
    elif arg == '-depthsrz':
        rm = domain_bounds[1]
        drz = rm - ri
        strings = args[i+1].split()
        for st in strings:
            rval = rm - float(st)*drz
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
    elif arg == '-bcz': # try to estimate the real BCZ (and mark it)
                        # from where the enthalpy flux first goes negative
        mark_bcz = True
    elif arg == '-lw':
        lw = float(args[i+1])
    elif arg == '-dpi':
        dpi = float(args[i+1])
    elif arg == '-rminw':
        rminwindow = float(args[i+1])
    elif arg == '-rmaxw':
        rmaxwindow = float(args[i+1])
    elif arg == '-nrec':
        nrec = int(args[i+1])
    elif arg == '-nskip':
        nskip = int(args[i+1])
    elif arg == '-ntot':
        ntot = int(args[i+1])

if not ntot is None: # This overrides nskip even if user specified it
    nskip = (index_last - index_first)//ntot

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
plotdir = dirname + '/plots/eflux_radial_vs_time/'
if not os.path.isdir(plotdir):
    os.makedirs(plotdir)

# Determine the simulation is magnetic
magnetism = get_parameter(dirname, 'magnetism')

# Plotting loop
print ('Plotting Shell_Avgs files %s through %s ...'\
       %(file_list[index_first], file_list[index_last]))

# Read in grid info
print ("Getting grid info from Shell_Avgs/" + file_list[index_first])
sh0 = Shell_Avgs(radatadir + file_list[index_first], '')
rr = sh0.radius
nr = sh0.nr
ri, ro = np.min(rr), np.max(rr)
# Make lstar = lsun unless otherwise specified in main_input
try:
    # First see if we can get c_10 from equation_coefficients:
    try:
        eq = equation_coefficients()
        eq.read(dirname + '/equation_coefficients')
        lstar = eq.constants[9]
        print("Got luminosity from 'equation_coefficients' file")
    except: # otherwise get "luminosity" from main_input
        lstar = get_parameter(dirname, 'luminosity')
        print ("Got luminosity from 'main_input' file")
except:
    lstar = lsun
    print ("Cannot find luminosity in either 'equation_coefficients'")
    print("or 'main_input' files. Setting luminosity to lsun.")

for i in range(index_first, index_last + 1, nskip):
    if i == index_first:
        sh = sh0
    else:   
        sh = Shell_Avgs(radatadir + file_list[i], '')
    sh = Shell_Avgs(radatadir + file_list[i], '')

    nrec_tot = sh.niter
    nstep = nrec_tot // nrec
    if nstep == 0:
        nstep = 1
    for j in range(0, nrec_tot, nstep):
        vals = sh.vals[:,0,:,j]
        lut = sh.lut

        iter_loc = sh.iters[j]
        t_loc = sh.time[j]

        qindex_hflux = lut[1433]
        qindex_cflux = lut[1470]
        qindex_kflux = lut[1923]
        qindex_vflux = lut[1935]
        qindex_eflux = lut[1455]

        hflux = vals[:, lut[1433]]
        eflux = vals[:, lut[1455]]
        cflux = vals[:, lut[1470]]
        kflux = vals[:, lut[1923]]
        vflux = -vals[:, lut[1935]]
        tflux = hflux + eflux + cflux + kflux + vflux # compute the total flux

        if magnetism:
            qindex_mflux = lut[2001] # this is actually (-4*pi) TIMES 
                                # the correct Poynting flux
            mflux = -vals[:, qindex_mflux]/(4*np.pi)
            tflux += mflux

        fpr = 4*np.pi*rr**2

        # Compute the integrated fluxes
        hflux_int = hflux*fpr
        eflux_int = eflux*fpr
        cflux_int = cflux*fpr
        kflux_int = kflux*fpr
        vflux_int = vflux*fpr
        tflux_int = tflux*fpr

        if magnetism:
            mflux_int = mflux*fpr

        # Create the plot; start with plotting all the energy fluxes

        # User can specify what to normalize the radius by
        # By default, normalize by the solar radius
        if rnorm is None:
            rr_n = rr/rsun
        else:
            rr_n = rr/rnorm                                           

        plt.plot(rr_n, hflux_int/lstar, label=r'$\rm{F}_{heat}$', linewidth=lw)
        plt.plot(rr_n, eflux_int/lstar, 'm', label = r'$\rm{F}_{enth}$',\
                linewidth=lw)
        plt.plot(rr_n, cflux_int/lstar, label = r'$\rm{F}_{cond}$', linewidth=lw)
        plt.plot(rr_n, kflux_int/lstar, label = r'$\rm{F}_{KE}$', linewidth=lw)
        plt.plot(rr_n, vflux_int/lstar, label = r'$\rm{F}_{visc}$', linewidth=lw)
        plt.plot(rr_n, tflux_int/lstar, label= r'$\rm{F}_{total}$',\
                linewidth=lw, color='black')
        if magnetism:
            plt.plot(rr_n, mflux_int/lstar, label=r'$\rm{F}_{Poynting}$',\
                linewidth=lw)

        # Mark zero line
        plt.plot(rr_n, np.zeros_like(rr_n), 'k--', linewidth=lw)

        # Get the y-axis in scientific notation
        plt.ticklabel_format(useMathText=True, axis='y', scilimits=(0,0))

        # Get ticks everywhere
        plt.minorticks_on()
        plt.tick_params(top=True, right=True, direction='in', which='both')

        # Set the x limits
        if xminmax_was_None:
            xminmax = np.min(rr_n), np.max(rr_n)
        # If user set -rminwindow or -rmaxwindow, this trumps the bounds
        if not rminwindow is None:
            rmin = np.min(rr_n)
            rmax = np.max(rr_n)
            Delta_r = rmax - rmin
            xminmax = rmin - rminwindow*Delta_r, rmin + rminwindow*Delta_r
        if not rmaxwindow is None:
            rmin = np.min(rr_n)
            rmax = np.max(rr_n)
            Delta_r = rmax - rmin
            xminmax = rmax - rmaxwindow*Delta_r, rmax + rmaxwindow*Delta_r
        plt.xlim(xminmax[0], xminmax[1])

        # Set the y-limits (the following values seem to "work well" for my models
        # so far...perhaps adjust this in the future. 

        if minmax_was_None:
            minmax = -0.7, 1.3
        if not rminwindow is None:
            minmax = -rminwindow, rminwindow
        if not rmaxwindow is None:
            minmax = -rmaxwindow, rmaxwindow
        plt.ylim(minmax[0], minmax[1])

        # Label the axes
        if rnorm is None:
            plt.xlabel(r'$r/R_\odot$',fontsize=12, **csfont)
        else:
            plt.xlabel(r'r/(%.1e cm)' %rnorm, fontsize=12, **csfont)

        # Try to find the BCZ from where enthalpy flux goes negative,
        # if desired avoid the outer boundary
        yvals = np.linspace(minmax[0], minmax[1], 100)
        if mark_bcz:
            irneg = np.argmin(eflux[20:] > 0) + 20
            irpos = np.argmin(eflux[20:] > 0) + 19
            rrneg = rr_n[irneg]
            rrpos = rr_n[irpos]
            efluxneg = eflux[irneg]
            efluxpos = eflux[irpos]
            slope =  (efluxpos - efluxneg)/(rrpos - rrneg)
            rbcz = rrpos - efluxpos/slope
            plt.plot(rbcz + np.zeros(100), yvals, 'k--', linewidth=lw)

        # Mark radii if desired
        if not rvals is None:
            for rval in rvals:
                if rnorm is None:
                    rval_n = rval/rsun
                else:
                    rval_n = rval/rnorm
        #        plt.ylim(ymin, ymax)
                plt.plot(rval_n + np.zeros(100), yvals, 'k--', linewidth=lw)


        plt.ylabel(r'$4\pi r^2\ \rm{\times \ (energy \ flux)}\ /\ L_*$',\
                fontsize=12, **csfont)

        # Make title
        # Create a title    
        if rotation:
            time_string = ('t = %.1f ' %(t_loc/time_unit)) + time_label +\
                    ' (1 ' + time_label + (' = %.2f days)'\
                    %(time_unit/86400.))
        else:
            time_string = ('t = %.3f ' %(t_loc/time_unit)) + time_label +\
                    ' (1 ' + time_label + (' = %.1f days)'\
            %(time_unit/86400.))

        the_title = dirname_stripped + '\n' + time_string
        if mark_bcz:
            the_title += ('\n' + r'$r_{BCZ}/\rm{rnorm} = %.3f$' %rbcz)

        plt.title(the_title, **csfont)

        # Create a see-through legend
        plt.legend(loc='lower left', shadow=True, ncol=3, fontsize=10)

        # Last command
        plt.tight_layout()

        # Save the plot
        savefile = plotdir + str(iter_loc).zfill(8) + '.png'
        print('Saving plot at ' + savefile)
        plt.savefig(savefile, dpi=dpi)
        plt.close()
