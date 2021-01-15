# Author: Loren Matilsky
# Created: 09/12/2019
# This script generates parity traces (amplitude-squared parity power (even
# or odd) vs time for a particular m-value and for a given variable at a
# a given radius, from the Shell_Spectra data
# the Rayleigh run directory indicated by [dirname]. To use  time-trace
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
        get_dict, rsun, sci_format, allthrees_start, get_satvals

from plotcommon import default_axes_1by1, axis_range
from get_parameter import get_parameter

# Get directory name and stripped_dirname for plotting purposes
dirname = sys.argv[1]
dirname_stripped = strip_dirname(dirname)

# Directory with data and plots, make the plotting directory if it doesn't
# already exist    
datadir = dirname + '/data/'

# Set defaults
varname = 'bp'
desired_rvals = ['all']
desired_mval = 0.
ir_vals = None
showplot = False
rnorm = None
minmax = None
the_file = get_widest_range_file(datadir, 'parity_vs_m')
include_tl = False

# Read command-line arguments (CLAs)
args = sys.argv[2:]
nargs = len(args)
for i in range(nargs):
    arg = args[i]
    if arg == '-var':
        varname = args[i+1]
    elif arg == '-rvals':
        my_str = args[i+1].split()
        desired_rvals = []
        for j in range(len(my_str)):
            desired_rvals.append(float(my_str[j]))
    elif arg == '-ir':
        my_str = args[i+1].split()
        ir_vals = []
        for j in range(len(my_str)):
            ir_vals.append(int(my_str[j]))
    elif arg == '-usefile':
        the_file = args[i+1]
        the_file = Shell_Spectra_file.split('/')[-1]
    elif arg == '-rnorm':
        rnorm = float(args[i+1])
    elif arg == '-minmax':
        minmax = float(args[i+1]), float(args[i+2])
    elif arg == '-show':
        showplot = True
    elif arg == '-mval':
        desired_mval = float(args[i+1])
    elif arg == '-tl':
        include_tl = True


# directory for plots (depends on whether logscale = True, so do this after
# command-line-arguments are read
plotdir = dirname + '/plots/Pspec/parity_vs_m/'
if (not os.path.isdir(plotdir)):
    os.makedirs(plotdir)

# Read in Pspec data
print ('Reading parity trace data from ' + datadir + the_file +\
        ' ...')
di = get_dict(datadir + the_file)
rvals = di['rvals']/rsun
mvals = di['mvals']
qv = di['qv']
nr = di['nr']
lvals = di['lvals']
nell = di['nell']
iter1, iter2 = di['iter1'], di['iter2']

# Get trace values
vals = di['vals']
Om0 = get_parameter(dirname, 'angular_velocity')
Prot = 2*np.pi/Om0
times = di['times']/Prot - allthrees_start

# Get time-latitude data
if include_tl:
    tl_file = get_widest_range_file(datadir, 'time-latitude')
    di = get_dict(datadir + tl_file)
    vals_tl = di['vals']

    times_tl = di['times']/Prot - allthrees_start
    tt_lat = di['tt_lat']
    times2, tt_lat2 = np.meshgrid(times_tl, tt_lat, indexing='ij')
    ntheta = di['ntheta']
    qvals_tl = np.array(di['qvals'])

if ir_vals is None:
    ir_vals = []
    if desired_rvals == ['all']:
        ir_vals = np.arange(nr)
    else:
        for i in range(len(desired_rvals)):
            rval = desired_rvals[i]
            ir_vals.append(np.argmin(np.abs(rvals - rval)))

# Get the m-index
im = np.argmin(np.abs(mvals - desired_mval))

# What is given is now what is desired
desired_rvals = rvals[ir_vals]

# Please please PLEASE, always sample power spectra and time-latitude data
# at the same depths ...

# Get power associated with desired quantity (should really have a contigency
# where the routine exits if the desired quantity isn't present
if not (varname == 'vtot' or varname == 'btot'):
    desired_qv = var_indices[varname]
    iq = np.argmin(np.abs(qv - desired_qv))
    if include_tl:
        iq_tl = np.argmin(np.abs(qvals_tl - desired_qv))
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
    iq_vals = []
    for desired_qv in desired_qv_vals:
        iq_vals.append(np.argmin(np.abs(qv - desired_qv)))

for ir in range(len(ir_vals)):
    rval =  desired_rvals[ir]
    savedir = plotdir + varname + ('/rval%0.3f/' %rval)
    if (not os.path.isdir(savedir)):
        os.makedirs(savedir)

    savename = savedir + ('m%04d_' %im) +\
            str(iter1).zfill(8) + '_' + str(iter2).zfill(8) + '.png'
    print('Plotting parity_vs_time: ' + varname +\
            (', r/rsun = %0.3f (ir = %02i), ' %(rval, ir_vals[ir])) + ' ...')
    if varname == 'vtot' or varname == 'btot':
        parity_even = vals[:, im, ir, iq_vals[0], 0] +\
                vals[:, im, ir, iq_vals[1], 0] + vals[:, im, ir, iq_vals[2], 0]
        parity_odd = vals[:, im, ir, iq_vals[0], 1] +\
                vals[:, im, ir, iq_vals[1], 1] + vals[:, im, ir, iq_vals[2], 1]
        parity_even_allm = np.sum(vals[:, :, ir, iq_vals[0], 0] +\
                vals[:, :, ir, iq_vals[1], 0] +\
                vals[:, :, ir, iq_vals[2], 0], axis=1)
        parity_odd_allm = np.sum(vals[:, :, ir, iq_vals[0], 1] +\
                vals[:, :, ir, iq_vals[1], 1] +\
                vals[:, :, ir, iq_vals[2], 1], axis=1)

        # Time-latitude stuff...don't mess with vtot/btot for now...

    else:
        parity_even = vals[:, im, ir, iq, 0]
        parity_odd = vals[:, im, ir, iq, 1]
        parity_even_allm = np.sum(vals[:, :, ir, iq, 0], axis=1)
        parity_odd_allm = np.sum(vals[:, :, ir, iq, 1], axis=1)
        if include_tl:
            quant_tl = vals_tl[:, :, ir, iq_tl]

    power_allm = parity_even_allm + parity_odd_allm
    # compute the actual parity
    parity = (parity_even - parity_odd)/(parity_even + parity_odd)

    # Now make the plot, finally
    plt.figure(figsize=(7.25, 3))
    if include_tl:
        fig, axs = plt.subplots(2, 1, figsize=(7.25, 4), sharex=True)
        plt.sca(axs[0])
    lw = 0.3
    fs = 7
    plt.plot(times, parity_even/(parity_even + parity_odd), label='even',\
            linewidth=lw)
    plt.plot(times, parity_odd/(parity_even + parity_odd), label='odd',\
            linewidth=lw)
    plt.plot(times, parity, label='parity', linewidth=lw)

    # set bounds
    plt.xlim(times[0], times[-1])
    plt.ylim(-1, 1)

    # make legend
    plt.legend(fontsize=fs)

    # label axes
    plt.xlabel('time (rotations)', fontsize=fs)

    # Get ticks everywhere
    plt.minorticks_on()
    plt.tick_params(top=True, right=True, direction='in', which='both')

    # Make title
    title = varlabel + '     ' + (r'$r/R_\odot\ =\ %0.3f$' %rval) +\
            '     ' + (r'$m = %i$' %im) + '\n' +\
            ('%08i to %08i' %(iter1, iter2))    
    plt.title(title, fontsize=fs)
    # set tick fontsize
    plt.tick_params(labelsize=fs)

    if include_tl:
        min_tl, max_tl = get_satvals(quant_tl)
        axs[1].pcolormesh(times2, tt_lat2, quant_tl, vmin=min_tl,\
                vmax=max_tl, cmap='RdYlBu_r')

    # Final command
    plt.tight_layout()

    print ('Saving ' + plotdir + savename + ' ...')

    plt.savefig(savename, dpi=300)
    if showplot:
        plt.show()
    plt.close()
