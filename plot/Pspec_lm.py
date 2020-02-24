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
# [dirname]/plots/Pspec/Pspec_lm/[varname]_[rval]_[first iter]_[last iter].png

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
        get_dict, rsun, sci_format, get_satvals, rms

from plotcommon import default_axes_1by1, axis_range, xy_grid
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
the_file = get_widest_range_file(datadir, 'Shell_Spectra')
lminmax = None

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
    elif arg == '-nolog':
        logscale = False
    elif arg == '-show':
        showplot = True
    elif arg == '-lminmax':
        lminmax = float(args[i+1]), float(args[i+2])

# directory for plots (depends on whether logscale = True, so do this after
# command-line-arguments are read
plotdir = dirname + '/plots/Pspec/Pspec_lm/'
if not logscale:
    plotdir = plotdir + 'notlog/'
if (not os.path.isdir(plotdir)):
    os.makedirs(plotdir)

# Read in Pspec data
print ('Reading Shell_Spectra data from ' + datadir + the_file +\
        ' ...')
di = get_dict(datadir + the_file)
rvals = di['rvals']/rsun
qv = di['qv']
lut = di['lut']
nr = di['nr']
lvals = di['lvals']
mvals = di['mvals']
nell = di['nell']
nm = di['nm']

if not lminmax is None:
    il1 = np.argmin(np.abs(lvals - lminmax[0]))
    il2 = np.argmin(np.abs(lvals - lminmax[1]))
else:
    il1, il2 = 0, nell - 1

lvals = lvals[il1:il2+1]
mvals = mvals[il1:il2+1]

lvals_2d, mvals_2d = np.meshgrid(lvals, mvals, indexing='ij')
lvals_2d_new, mvals_2d_new = xy_grid(lvals_2d, mvals_2d)
iter1, iter2 = di['iter1'], di['iter2']

# Get full power
fullpower = di['fullpower']

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

# Get power in l,m space for desired variable
if varname == 'vtot' or varname == 'btot':
    power = np.abs(fullpower[:, :, :, iq_vals[0]])**2 +\
            np.abs(fullpower[:, :, :, iq_vals[1]])**2 +\
            np.abs(fullpower[:, :, :, iq_vals[2]])**2 

else:
    power = np.abs(fullpower[:, :, :, iq])**2

# Now make plots
for ir in range(len(ir_vals)):
    rval =  desired_rvals[ir]
    savename = varname + ('_rval%0.3f' %rval) + '_' +\
          str(iter1).zfill(8) + '_' + str(iter2).zfill(8) + '.png'
    print('Plotting Pspec_lm: ' + varname + (', r/rsun = %0.3f (ir = %02i), '\
            %(rval, ir_vals[ir])) + ' ...')
    power_loc = power[il1:il2+1, il1:il2+1, ir]

    print ('max power: %1.1e' %np.max(power_loc))
    # Get minmax, if not specified
    if minmax is None:
        minmax_loc = 0, 5*rms(power_loc)
    else:
        minmax_loc = minmax
    # Make plot
    plt.pcolormesh(lvals_2d_new, mvals_2d_new, power_loc, cmap='Greys',\
            vmin=minmax_loc[0], vmax=minmax_loc[1])

    # set bounds
#    plt.xlim(0.5, nell - 0.5)
#    plt.ylim(0.5, nm - 0.5)

    # label axes
    plt.xlabel('l-value')
    plt.ylabel('m-value')

    # Get ticks everywhere
    plt.minorticks_on()
    plt.tick_params(top=True, right=True, direction='in', which='both')

    # Get colorbar
    cbar = plt.colorbar()

    # Make title
    # Compute l_rms and m_rms
    l_rms = np.sum(power_loc*lvals_2d)/np.sum(power_loc)
    m_rms = np.sum(power_loc*mvals_2d)/np.sum(power_loc)

    title = varlabel + '     ' + (r'$r/R_\odot\ =\ %0.3f$' %rval) + '\n' +\
            (r'$l_{\rm{rms}} = %.1f$' %l_rms) + '     ' +\
            (r'$m_{\rm{rms}} = %.1f$' %m_rms)+ '\n' +\
            ('%08i to %08i' %(iter1, iter2))    
    plt.title(title)

    # Final command
    plt.tight_layout()

    print ('Saving ' + plotdir + savename + ' ...')

    plt.savefig(plotdir + savename, dpi=300)
    if showplot:
        plt.show()
    plt.close()
