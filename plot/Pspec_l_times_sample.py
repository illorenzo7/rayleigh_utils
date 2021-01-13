# Author: Loren Matilsky
# Created: 01/30/2020
# This script generates powerspectra (amplitude-squared power vs spherical-
# harmonic-degree l) for m = 0, m != 0, and total for a given variable at a
# a given radius, from the Shell_Spectra data
# Plots the last 10 Shell_Spectra (by default) at an instant, saving the 
# plots in directory
# plots/Pspec_l_vs_time/

import matplotlib as mpl
mpl.use('TkAgg')
import matplotlib.pyplot as plt
plt.rcParams['mathtext.fontset'] = 'dejavuserif'
csfont = {'fontname':'DejaVu Serif'}
import numpy as np
import sys, os
sys.path.append(os.environ['raco'])
sys.path.append(os.environ['rapp'])
from varprops import texunits, texlabels, var_indices
from common import *
from rayleigh_diagnostics import Shell_Spectra
from time_scales import compute_Prot, compute_tdt
from translate_times import translate_times

# Get directory name and stripped_dirname for plotting purposes
dirname = sys.argv[1]
dirname_stripped = strip_dirname(dirname)

# Directory with Shell_Spectra
radatadir = dirname + '/Shell_Spectra/'

# Set defaults
ir = 0 # by default plot just below the surface
rval = None # can also find ir by finding the closest point
            # to a local radius divided by rsun
varname = 'vr'
logscale = True
minmax = None

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

# Change the defaults using the CLAs
for i in range(nargs):
    arg = args[i]
    if arg == '-var':
        varname = args[i+1]
    elif arg == '-ir':
        ir = int(args[i+1])
    elif arg == '-rval':
        rval = float(args[i+1])
    elif arg == '-minmax':
        minmax = float(args[i+1]), float(args[i+2])
    elif arg == '-nolog':
        logscale = False
    elif arg == '-show':
        showplot = True

# Get the baseline time unit
rotation = get_parameter(dirname, 'rotation')
if rotation:
    time_unit = compute_Prot(dirname)
    time_label = r'$\rm{P_{rot}}$'
else:
    time_unit = compute_tdt(dirname)
    time_label = r'$\rm{TDT}$'

# directory for plots (depends on whether logscale = True, so do this after
# command-line-arguments are read
if logscale:
    plotdir = dirname + '/plots/Pspec_l/times_sample/'
else:
    plotdir = dirname + '/plots/Pspec_l/times_sample_notlog/'

if not os.path.isdir(plotdir):
    os.makedirs(plotdir)

# Read in Pspec data, use first file for grid stuff
spec0 = Shell_Spectra(radatadir + file_list[index_first], '')
rvals = spec0.radius/rsun
qv = spec0.qv
nr = spec0.nr
nell = spec0.nell
lvals = np.arange(nell + 0.)

# Plotting loop
print ('Plotting Shell_Spectra files %s through %s ...'\
       %(file_list[index_first], file_list[index_last]))
print('Saving plots to '  + plotdir)

for i in range(index_first, index_last + 1):
    fname = file_list[i]
    if i == index_first:
        spec = spec0
    else:   
        spec = Shell_Spectra(radatadir + fname, '')

    for j in range(spec.niter):
        # Get local time (in seconds)
        t_loc = spec.time[j]
        iter_loc = spec.iters[j]

        # Get the l power (various units)
        # take square root, so has same unit as variable itself
        lpower = np.sqrt(spec.lpower[:, :, :, j, :])
        if varname in ['vr', 'vt', 'vp', 'vtot']: # convert to m/s
            lpower /= 100.

        # Find desired radius (by default ir=0--near outer surface)
        if not rval is None:
            ir = np.argmin(np.abs(spec.radius/rsun - rval))
        rval = spec.radius[ir]/rsun 
        # in any case, this is the actual rvalue we get

        # Get power associated with desired quantity 
        # (should really have a contigency
        # where the routine exits if the desired quantity isn't present)
        if (varname == 'vtot' or varname == 'btot' or varname == 'omtot'):
            if varname == 'vtot':
                qv_vals = [1, 2, 3]
                varlabel = r'$|\mathbf{v}|$'
                units = r'$\rm{m}\ \rm{s}^{-1}$'
            elif varname == 'omtot':
                qv_vals = [301, 302, 303]
                varlabel = r'$|\mathbf{\omega}|$'
                units = r'$\rm{s^{-1}}$'
            elif varname == 'btot':
                qv_vals = [801, 802, 803]
                varlabel = r'$|\mathbf{B}|$'
                units = r'$\rm{G}$'
            field_lpower = np.sqrt(lpower[:, ir,\
                    spec.lut[qv_vals[0]], :]**2. +\
                    lpower[:, ir, spec.lut[qv_vals[1]], :]**2. +\
                    lpower[:, ir, spec.lut[qv_vals[2]], :]**2.)
        else:
            field_lpower = lpower[:, ir, spec.lut[var_indices[varname]], :]
            varlabel = texlabels[varname]
            units = texunits[varname] 

        # Get saturation values if not specified--avoid the two extreme 
        # l-values in calculation
        if minmax is None:
            lpower_cut = field_lpower[1:-1, :]
            amp_min, amp_max = np.min(lpower_cut), np.max(lpower_cut)
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

        # Make the savename like for Mollweide times sample
        savename = 'Pspec_l_' + varname + ('_rval%0.3f' %rval) + '_iter' +\
                str(iter_loc).zfill(8) + '.png'
        print('Plotting: ' + savename)

        # Get the tot, m=0, and m!=0 power
        lpower_tot = field_lpower[:, 0]
        lpower_m0 = field_lpower[:, 1]
        lpower_mnot0 = field_lpower[:, 2]

        if logscale:
            plt.loglog(lvals, lpower_tot, label=r'$\rm{tot}$')
            plt.loglog(lvals, lpower_m0, label=r'$\rm{m = 0}$')
            plt.loglog(lvals, lpower_mnot0, label=r'$\rm{|m| > 0}$')
        else:
            plt.plot(lvals, lpower_tot, label=r'$\rm{tot}$')
            plt.plot(lvals, lpower_m0, label=r'$\rm{m = 0}$')
            plt.plot(lvals, lpower_mnot0, label=r'$\rm{|m| > 0}$')

        # set bounds
        plt.xlim(1, nell - 1)
        plt.ylim(ymin, ymax)

        # make legend
        plt.legend()

        # label axes
        plt.xlabel(r'$\ell$')
        plt.ylabel(r'$\sqrt{P_\ell}$' + ' (' + units + ')')

        # Get ticks everywhere
        plt.minorticks_on()
        plt.tick_params(top=True, right=True, direction='in', which='both')

        # Make title
        # Compute l_rms
        l_rms = np.sum(lpower_mnot0**2*lvals)/np.sum(lpower_mnot0**2)

        if rotation:
            time_string = ('t = %.1f ' %(t_loc/time_unit)) + time_label +\
                    ' (1 ' + time_label + (' = %.2f days)'\
                    %(time_unit/86400.))
        else:
            time_string = ('t = %.3f ' %(t_loc/time_unit)) + time_label +\
                    ' (1 ' + time_label + (' = %.1f days)'\
                    %(time_unit/86400.))

        title = dirname_stripped +\
            '\n' + r'$\rm{Pspec\ vs.\ \ell}$' + '     '  + time_string +\
            '\n' + varlabel + '     ' + (r'$r/R_\odot\ =\ %0.3f$' %rval) +\
                '     ' + (r'$\ell_{\rm{rms},\ m\neq0} = %.1f$' %l_rms)
        plt.title(title, **csfont)

        # Final commands
        plt.tight_layout()
        plt.savefig(plotdir + savename, dpi=300)
        plt.close()
