# Author: Loren Matilsky
# Created: 09/12/2019
# This script generates powerspectra (amplitude-squared power) vs spherical-
# harmonic modes l,m for all l, m at 
# a given radius, from the Shell_Spectra data in
# the Rayleigh run directory indicated by [dirname]. To use  time-averaged 
# Shell_Spectra file different than the one associated with the 
# longest averaging range, use -usefile [complete path-name of desired vavg 
# file]
# Saves plot in
# [dirname]/plots/specav_lm/rvals_sample/[varname]_[first iter]_[last iter]_[rval].png

import matplotlib as mpl
mpl.use('TkAgg')
from matplotlib import ticker, colors
import matplotlib.pyplot as plt
plt.rcParams['mathtext.fontset'] = 'dejavuserif'
csfont = {'fontname':'DejaVu Serif'}
import numpy as np
import sys, os
sys.path.append(os.environ['raco'])
sys.path.append(os.environ['rapp'])
from varprops import texunits, texlabels, var_indices
from common import *
        get_iters_from_file, get_dict, rsun, sci_format, get_satvals, rms
from plotcommon import axis_range, xy_grid
from time_scales import compute_Prot, compute_tdt
from translate_times import translate_times
from rayleigh_diagnostics import Shell_Spectra

# Get directory name and stripped_dirname for plotting purposes
dirname = sys.argv[1]
dirname_stripped = strip_dirname(dirname)

# Directory with data and plots, make the plotting directory if it doesn't
# already exist    
datadir = dirname + '/data/'

# Data with Shell_Slices
radatadir = dirname + '/Shell_Spectra/'
file_list, int_file_list, nfiles = get_file_lists(radatadir)

# Set defaults
ir = 0 # by default plot just below the surface
rval = None # can also find ir by finding the closest point
            # to a local radius divided by rsun
iiter = 0 # by default get info from first Shell_Spectra
varname = 'vr'
minmax = None
logscale = True
fs = 12.
tag = ''

# Get specific range desired for plotting
args = sys.argv[2:]
nargs = len(args)
the_tuple = get_desired_range(int_file_list, args)
if the_tuple is None: # By default plot the last 10 Shell_Slices
    index_first, index_last = nfiles - 11, nfiles - 1  
else:
    index_first, index_last = the_tuple

# Read other command-line arguments (CLAs)
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
    elif arg == '-tag':
        tag = '_' + args[i+1]
    elif arg == '-nolog':
        logscale = False

# Get the baseline time unit
rotation = get_parameter(dirname, 'rotation')
if rotation:
    time_unit = compute_Prot(dirname)
    time_label = r'$\rm{P_{rot}}$'
else:
    time_unit = compute_tdt(dirname)
    time_label = r'$\rm{TDT}$'

# Get general info from first file
fname = file_list[iiter]

# Read in last shell spectrum
spec = Shell_Spectra(radatadir + fname, '')
rvals = spec.radius/rsun
qv = spec.qv
lut = spec.lut
nr = spec.nr
nell = spec.nell
nm = spec.nm
lvals = np.arange(nell, dtype='float')

# Get local time (in seconds)
t_loc = spec.time[0]

# directory for plots 
plotdir = dirname + '/plots/spec_l/times_sample' + tag + '/'
if not os.path.isdir(plotdir):
    os.makedirs(plotdir)

# Now make plots
for i in range(index_first, index_last + 1):
    # Read in desired shell spectrum
    fname = file_list[i]
    spec = Shell_Spectra(radatadir + fname, '')

    # Find desired radius (by default ir=0--near outer surface)
    if not rval is None:
        ir = np.argmin(np.abs(spec.radius/rsun - rval))
    rval = spec.radius[ir]/rsun 
    # in any case, this is the actual rvalue we get

    # Loop over the slices in the file 
    for j in range(spec.niter):
        # Get the l-power
        lpower = spec.lpower[:, :, :, j, :] 

        # Get local time (in seconds)
        t_loc = spec.time[j]
        iter_loc = spec.iters[j]

        # Savename
        savename = 'spec_l_' + varname + ('_rval%0.3f' %rval) +\
                '_iter' + str(iter_loc).zfill(8) + '.png'
        print('Plotting: ' + savename)
        if varname == 'vtot' or varname == 'btot' or varname == 'omtot':
            if varname == 'vtot':
                desired_qv_vals = [1, 2, 3]
                varlabel = r'$|\mathbf{v}|$'
                units = r'$\rm{m}\ \rm{s}^{-1}$'
            elif varname == 'btot':
                desired_qv_vals = [801, 802, 803]
                varlabel = r'$|\mathbf{B}|$'
                units = r'$\rm{G}$'
            elif varname == 'omtot':
                desired_qv_vals = [301, 302, 303]
                varlabel = r'$|\mathbf{\omega}|$'
                units = r'$\rm{s}^{-1}$'
            iq_vals = []
            for desired_qv in desired_qv_vals:
                iq_vals.append(np.argmin(np.abs(qv - desired_qv)))
            lpower_tot = lpower[:, ir, iq_vals[0], 0] +\
                    lpower[:, ir, iq_vals[1], 0] +\
                    lpower[:, ir, iq_vals[2], 0]
            lpower_m0 = lpower[:, ir, iq_vals[0], 1] +\
                    lpower[:, ir, iq_vals[1], 1] +\
                    lpower[:, ir, iq_vals[2], 1]
            lpower_mnot0 = lpower[:, ir, iq_vals[0], 2] +\
                    lpower[:, ir, iq_vals[1], 2] +\
                    lpower[:, ir, iq_vals[2], 2]
        else:
            desired_qv = var_indices[varname]
            iq = np.argmin(np.abs(qv - desired_qv))
            varlabel = texlabels[varname]
            units = texunits[varname] 
            lpower_tot = lpower[:, ir, iq, 0]
            lpower_m0 = lpower[:, ir, iq, 1]
            lpower_mnot0 = lpower[:, ir, iq, 2]

        if logscale:
            plt.loglog(lvals, lpower_tot, label='tot')
            plt.loglog(lvals, lpower_m0, label='m = 0')
            plt.loglog(lvals, lpower_mnot0, label='|m| > 0')
        else:
            plt.plot(lvals, lpower_tot, label='tot')
            plt.plot(lvals, lpower_m0, label='m = 0')
            plt.plot(lvals, lpower_mnot0, label='|m| > 0')

        # set bounds
        plt.xlim(1, nell - 1)
        # Get saturation values if not specified--avoid the two extreme
        # l values in calculation
        if minmax is None:
            # ignore the first and last l-values
            if varname == 'vtot' or varname == 'btot' or varname == 'omtot':
                lpower_cut = lpower[1:-1, ir, iq_vals[0], :] +\
                    lpower[1:-1, ir, iq_vals[1], :] +\
                    lpower[1:-1, ir, iq_vals[2], :]
            else:
                lpower_cut = lpower[1:-1, ir, iq, :]
            if logscale:
                minmax_loc = np.min(lpower_cut)/2., np.max(lpower_cut)*2.
            else:
                minmax_loc = 0., np.max(lpower_cut)*1.1
        else:
            minmax_loc = minmax
        plt.ylim(minmax_loc[0], minmax_loc[1])

        # make legend
        plt.legend()

        # label axes
        plt.xlabel(r'$\ell$')
        plt.ylabel(r'$P_\ell\equiv\Sigma_m|\hat{\psi}_{\ell m}|^2}$' +\
                ' (' + units + r'$)^2$')

        # Get ticks everywhere
        plt.minorticks_on()
        plt.tick_params(top=True, right=True, direction='in', which='both')

        # Make title
        # Label the current time in the title
        if rotation:
            time_string = ('t = %.1f ' %(t_loc/time_unit)) + time_label +\
                    ' (1 ' + time_label + (' = %.2f days)'\
                    %(time_unit/86400.))
        else:
            time_string = ('t = %.3f ' %(t_loc/time_unit)) + time_label +\
                    ' (1 ' + time_label + (' = %.1f days)'\
                    %(time_unit/86400.))

        # Make title
        # Compute l_rms
        l_rms = np.sum(lpower_mnot0**2*lvals)/np.sum(lpower_mnot0**2)

        title = varlabel + '     ' + (r'$r/R_\odot\ =\ %0.3f$' %rval) +\
                '     ' + (r'$l_{\rm{rms},\ m\neq0} = %.1f$' %l_rms) +\
                '\n' + time_string
        plt.title(title)

        # Final command
        plt.tight_layout()

        print ('Saving ' + plotdir + savename)

        plt.savefig(plotdir + savename, dpi=300)
        plt.close()
