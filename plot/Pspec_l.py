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

from common import strip_dirname, get_widest_range_file, get_iters_from_file,\
        get_dict, rsun, sci_format

from plotcommon import default_axes_1by1, axis_range
from get_parameter import get_parameter

# Get directory name and stripped_dirname for plotting purposes
dirname = sys.argv[1]
dirname_stripped = strip_dirname(dirname)

# Directory with data and plots, make the plotting directory if it doesn't
# already exist    
datadir = dirname + '/data/'

# Set defaults
varname = 'vr'
desired_rvals = ['all']
ir_vals = None
logscale = True
showplot = False
rnorm = None
minmax = None
Shell_Spectra_file = get_widest_range_file(datadir, 'Shell_Spectra')

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
        Shell_Spectra_file = args[i+1]
        Shell_Spectra_file = Shell_Spectra_file.split('/')[-1]
    elif arg == '-rnorm':
        rnorm = float(args[i+1])
    elif arg == '-minmax':
        minmax = float(args[i+1]), float(args[i+2])
    elif arg == '-nolog':
        logscale = False
    elif arg == '-show':
        showplot = True

# directory for plots (depends on whether logscale = True, so do this after
# command-line-arguments are read
plotdir = dirname + '/plots/Pspec/Pspec_l/'
if not logscale:
    plotdir = plotdir + 'notlog/'
if (not os.path.isdir(plotdir)):
    os.makedirs(plotdir)

# Read in Pspec data
print ('Reading Shell_Spectra data from ' + datadir + Shell_Spectra_file +\
        ' ...')
di = get_dict(datadir + Shell_Spectra_file)
rvals = di['rvals']/rsun
qv = di['qv']
nr = di['nr']
lvals = di['lvals']
nell = di['nell']
iter1, iter2 = di['iter1'], di['iter2']

# Scale so power corresponds to representative field strength
lpower = di['lpower']

fourpi = 4*np.pi
lpower = np.sqrt(nell*lpower/fourpi)

if varname in ['vr', 'vt', 'vp', 'vtot']:
    lpower /= 100 # change cm/s --> m/s

if ir_vals is None:
    ir_vals = []
    if desired_rvals == ['all']:
        ir_vals = np.arange(nr)
    else:
        for i in range(len(desired_rvals)):
            rval = desired_rvals[i]
            ir_vals.append(np.argmin(np.abs(rvals - rval)))

# What is given is now what is desired
desired_rvals = rvals[ir_vals]

# Get power associated with desired quantity (should really have a contigency
# where the routine exits if the desired quantity isn't present
if not (varname == 'vtot' or varname == 'btot'):
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
    iq_vals = []
    for desired_qv in desired_qv_vals:
        iq_vals.append(np.argmin(np.abs(qv - desired_qv)))

# Get saturation values if not specified--avoid the two extreme l values in 
# calculation
if minmax is None:
    # ignore the first and last l-values
    if varname == 'vtot' or varname == 'btot':
        lpower_cut = np.sqrt(lpower[1:-1, ir_vals, iq_vals[0], :]**2 +\
            lpower[1:-1, ir_vals, iq_vals[1], :]**2 +\
            lpower[1:-1, ir_vals, iq_vals[2], :]**2)
    else:
        lpower_cut = lpower[1:-1, ir_vals, iq, :]
    if logscale:
        minmax = np.min(lpower_cut)/2., np.max(lpower_cut)*2.
    else:
        minmax = 0., np.max(lpower_cut)*1.1

for ir in range(len(ir_vals)):
    rval =  desired_rvals[ir]
    savename = varname + ('_rval%0.3f' %rval) + '_' +\
          str(iter1).zfill(8) + '_' + str(iter2).zfill(8) + '.png'
    print('Plotting Pspec_l: ' + varname + (', r/rsun = %0.3f (ir = %02i), '\
            %(rval, ir_vals[ir])) + ' ...')
    if varname == 'vtot' or varname == 'btot':
        lpower_tot = np.sqrt(lpower[:, ir, iq_vals[0], 0]**2 +\
                lpower[:, ir, iq_vals[1], 0]**2 +\
                lpower[:, ir, iq_vals[2], 0]**2)
        lpower_m0 = np.sqrt(lpower[:, ir, iq_vals[0], 1]**2 +\
                lpower[:, ir, iq_vals[1], 1]**2 +\
                lpower[:, ir, iq_vals[2], 1]**2)
        lpower_mnot0 = np.sqrt(lpower[:, ir, iq_vals[0], 2]**2 +\
                lpower[:, ir, iq_vals[1], 2]**2 +\
                lpower[:, ir, iq_vals[2], 2]**2)
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
    plt.ylim(minmax[0], minmax[1])

    # make legend
    plt.legend()

    # label axes
    plt.xlabel('l-value')
    plt.ylabel(r'$\sqrt{(l_{\rm{max}} + 1)P_l/4\pi}$' + ' (' + units + ')')

    # Get ticks everywhere
    plt.minorticks_on()
    plt.tick_params(top=True, right=True, direction='in', which='both')

    # Make title
    # Compute l_rms
    l_rms = np.sum(lpower_mnot0**2*lvals)/np.sum(lpower_mnot0**2)

    title = varlabel + '     ' + (r'$r/R_\odot\ =\ %0.3f$' %rval) +\
            '     ' + (r'$l_{\rm{rms},\ m\neq0} = %.1f$' %l_rms) + '\n' +\
            ('%08i to %08i' %(iter1, iter2))    
    plt.title(title)

    # Final command
    plt.tight_layout()

    print ('Saving ' + plotdir + savename + ' ...')

    plt.savefig(plotdir + savename, dpi=300)
    if showplot:
        plt.show()
    plt.close()