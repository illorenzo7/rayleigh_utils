# Author: Loren Matilsky
# Created: 09/12/2019
# This script generates powerspectra (amplitude-squared power vs spherical-
# harmonic-degree l) for m = 0, m != 0, and total for a given variable at a
# a given radius, from the Shell_Spectra data
# the Rayleigh run directory indicated by [dirname]. To use  time-averaged 
# Shell_Spectra file different than the one associated with the 
# longest averaging range, use -usefile [complete path-name of desired vavg 
# file]
# Saves plot in
# [dirname]_Pspec_l_[varname]_[rval]_[first iter]_[last iter].npy

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
        get_iters_from_file, get_dict, rsun, sci_format

from plotcommon import default_axes_1by1, axis_range
from time_scales import compute_Prot, compute_tdt
from translate_times import translate_times

# Get directory name and stripped_dirname for plotting purposes
dirname = sys.argv[1]
dirname_stripped = strip_dirname(dirname)

# Directory with data and plots, make the plotting directory if it doesn't
# already exist    
datadir = dirname + '/data/'

# Set defaults
varname = 'vr'
desired_rvals = ['all']
logscale = True
Shell_Spectra_file = get_widest_range_file(datadir, 'Shell_Spectra')
tag = ''
ir = 0 # by default plot just below the surface
rval = None # can also find ir by finding the closest point
varlist = None

# Read command-line arguments (CLAs)
args = sys.argv[2:]
nargs = len(args)
for i in range(nargs):
    arg = args[i]
    if arg == '-ir':
        ir = int(args[i+1])
    elif arg == '-rval':
        rval = float(args[i+1])
    elif arg == '-usefile':
        Shell_Spectra_file = args[i+1]
        Shell_Spectra_file = Shell_Spectra_file.split('/')[-1]
    elif arg == '-tag':
        tag = '_' + args[i+1]
    elif arg == '-nolog':
        logscale = False
    elif arg == '-vars':
        varlist = args[i+1].split()

# directory for plots (depends on whether logscale = True, so do this after
# command-line-arguments are read
plotdir = dirname + '/plots/specav_l/vars_sample' + tag + '/'
if not logscale:
    plotdir = plotdir + 'notlog/'
if not os.path.isdir(plotdir):
    os.makedirs(plotdir)

# Read in Pspec data
print ('Reading Shell_Spectra data from ' + datadir + Shell_Spectra_file)
di = get_dict(datadir + Shell_Spectra_file)
rvals = di['rvals']/rsun
qv = di['qv']
nr = di['nr']
lvals = di['lvals']
nell = di['nell']
iter1, iter2 = di['iter1'], di['iter2']

# Find desired radius (by default ir=0--near outer surface)
if not rval is None:
    ir = np.argmin(np.abs(rvals - rval))
rval = rvals[ir] # in any case, this is the actual rvalue we get

# Get the time range in sec
t1 = translate_times(iter1, dirname, translate_from='iter')['val_sec']
t2 = translate_times(iter2, dirname, translate_from='iter')['val_sec']

# Get the baseline time unit
rotation = get_parameter(dirname, 'rotation')
if rotation:
    time_unit = compute_Prot(dirname)
    time_label = r'$\rm{P_{rot}}$'
else:
    time_unit = compute_tdt(dirname)
    time_label = r'$\rm{TDT}$'

# Get squared power vs. l
lpower = di['lpower']

if varname in ['vr', 'vt', 'vp', 'vtot']:
    lpower /= 1.0e4 # change (cm/s)^2 --> (m/s)^2

# Loop over vars and make plots
if varlist is None:
    varlist = ['vr', 'vt', 'vp', 'vtot', 'omr', 'omt', 'omp', 'omtot',\
            's', 'p']
    magnetism = get_parameter(dirname, 'magnetism')
    if magnetism:
        varlist.extend(['br', 'bt', 'bp', 'btot'])

# Now make plots
for varname in varlist:
    print('Plotting specav_l: ' + varname +\
            (', r/rsun = %0.3f (ir = %02i), ' %(rval, ir)))

    # Get power associated with desired quantity (should really have a
    # contigency where the routine exits if the desired quantity isn't
    # present
    if not (varname == 'vtot' or varname == 'btot' or varname == 'omtot'):
        desired_qv = var_indices[varname]
        iq = np.argmin(np.abs(qv - desired_qv))
        varlabel = texlabels[varname]
        units = texunits[varname] 
    else:
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

    savename = 'specav_l_' + str(iter1).zfill(8) + '_' +\
            str(iter2).zfill(8) +  ('_rval%0.3f_' %rval) +\
            varname + '.png'
    if varname == 'vtot' or varname == 'btot' or varname=='omtot':
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
    # Get saturation values--avoid the two extreme l values in calculation
    if varname == 'vtot' or varname == 'btot' or varname == 'omtot':
        lpower_cut = lpower[1:-1, ir, iq_vals[0], :] +\
            lpower[1:-1, ir, iq_vals[1], :] +\
            lpower[1:-1, ir, iq_vals[2], :]
    else:
        lpower_cut = lpower[1:-1, ir, iq, :]
    if logscale:
        minmax = np.min(lpower_cut)/2., np.max(lpower_cut)*2.
    else:
        minmax = 0., np.max(lpower_cut)*1.1
    plt.ylim(minmax[0], minmax[1])

    # make legend
    plt.legend()

    # label axes
    plt.xlabel(r'$\ell$')
    plt.ylabel(r'$P_\ell\equiv\Sigma_m|\hat{\psi}_{\ell m}|^2}$' +\
            ' (' + units + r'$)^2$')

    # Get ticks everywhere
    plt.minorticks_on()
    plt.tick_params(top=True, right=True, direction='in', which='both')

    # Label averaging interval
    if rotation:
        time_string = ('t = %.1f to %.1f ' %(t1/time_unit, t2/time_unit))\
                + time_label + ' ' + (r'$\ (\Delta t = %.1f\ $'\
                %((t2 - t1)/time_unit)) + time_label + ')'
    else:
        time_string = ('t = %.3f to %.3f ' %(t1/time_unit, t2/time_unit))\
                + time_label + (r'$\ (\Delta t = %.3f\ $'\
                %((t2 - t1)/time_unit)) + time_label + ')'

    # Make title
    # Compute l_rms
    l_rms = np.sum(lpower_mnot0**2*lvals)/np.sum(lpower_mnot0**2)

    title = varlabel + '     ' + (r'$r/R_\odot\ =\ %0.3f$' %rval) +\
            '     ' + (r'$l_{\rm{rms},\ m\neq0} = %.1f$' %l_rms) + '\n' +\
            time_string
    plt.title(title)

    # Final command
    plt.tight_layout()

    print ('Saving ' + plotdir + savename)

    plt.savefig(plotdir + savename, dpi=300)
    plt.close()
