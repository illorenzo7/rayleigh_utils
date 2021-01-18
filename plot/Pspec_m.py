# Author: Loren Matilsky
# Created: 10/10/2019
# This script generates powerspectra (amplitude-squared power) vs spherical-
# harmonic modes m (summed over all l)  for
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

from common import *

from plotcommon import default_axes_1by1, axis_range, xy_grid

# Get directory name and stripped_dirname for plotting purposes
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

# Set defaults
varname = 'vr'
desired_rvals = ['all']
ir_vals = None
logscale = True
showplot = False
ylog = False
linear = False
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
    elif arg == '-lin':
        linear = True
        logscale = False
    elif arg == '-ylog':
        logscale = False
        ylog = True
    elif arg == '-show':
        showplot = True

# directory for plots (depends on whether logscale = True, so do this after
# command-line-arguments are read
plotdir = dirname + '/plots/Pspec/Pspec_m/'
if linear:
    plotdir = plotdir + 'notlog/'
elif ylog:
    plotdir = plotdir + 'ylog/'
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

iter1, iter2 = di['iter1'], di['iter2']

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
    varlabel = texlabels.get(varname, 'qval = ' + varname)
    units = texunits.get(varname, 'cgs') 
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

# Now sum power over all al modes
power_vs_m = np.sum(power, axis=0)


# Now make plots
for ir in range(len(ir_vals)):
    rval =  desired_rvals[ir]
    savename = varname + ('_rval%0.3f' %rval) + '_' +\
          str(iter1).zfill(8) + '_' + str(iter2).zfill(8) + '.png'
    print('Plotting Pspec_lm: ' + varname + (', r/rsun = %0.3f (ir = %02i), '\
            %(rval, ir_vals[ir])) + ' ...')
   
    power_loc = power_vs_m[:, ir]/np.sum(power_vs_m[:, ir], axis=0)
    # Get saturation values if not specified--avoid the two extreme m values in 
    # calculation
    if minmax is None:
        if logscale or ylog:
                minmax = np.min(power_loc[1:-1])/2., np.max(power_loc[1:-1])*2.
        elif linear:
           minmax = 0., np.max(power_loc[1:-1])*1.1

    # Make the plot
    if logscale:
        plt.loglog(mvals, power_loc, label='tot')
    elif linear:
        plt.plot(mvals, power_loc, label='tot')
    elif ylog:
        plt.plot(mvals, power_loc, label='tot')
        plt.yscale('log')

    # set bounds
    plt.xlim(1, nm - 1)
    plt.ylim(minmax[0], minmax[1])

    # label axes
    plt.xlabel('m-value')
    plt.ylabel('fractional power')

    # Get ticks everywhere
    plt.minorticks_on()
    plt.tick_params(top=True, right=True, direction='in', which='both')

    # Make title
    # Label averaging interval
    if rotation:
        time_string = ('t = %.1f to %.1f ' %(t1/time_unit, t2/time_unit))\
                + time_label + ' ' + (r'$\ (\Delta t = %.1f\ $'\
                %((t2 - t1)/time_unit)) + time_label + ')'
    else:
        time_string = ('t = %.3f to %.3f ' %(t1/time_unit, t2/time_unit))\
                + time_label + (r'$\ (\Delta t = %.3f\ $'\
                %((t2 - t1)/time_unit)) + time_label + ')'
    # Compute l_rms
    m_rms = np.sum(mvals*power_loc)/np.sum(power_loc)

    title = varlabel + '     ' + (r'$r/R_\odot\ =\ %0.3f$' %rval) +\
            '     ' + (r'$m_{\rm{rms}} = %.1f$' %m_rms) + '\n' +\
            time_string
    plt.title(title)

    # Final command
    plt.tight_layout()

    print ('Saving ' + plotdir + savename + ' ...')
    plt.savefig(plotdir + savename, dpi=300)
    if showplot:
        plt.show()
    plt.close()
