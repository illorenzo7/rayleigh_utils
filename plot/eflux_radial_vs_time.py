###############################################
# Author: Loren Matilsky
# Date created: 05/06/2020
#
# This script plots the radial energy fluxes as functions of
# radius using from the Shell_Avgs data, one for EACH TIME in the record
# (Note that this will typically generate a LOT of plots)

import matplotlib as mpl
mpl.use('TkAgg')
import matplotlib.pyplot as plt
plt.rcParams['mathtext.fontset'] = 'dejavuserif'
csfont = {'fontname':'DejaVu Serif'}
import numpy as np
import sys, os
sys.path.append(os.environ['rapp'])
sys.path.append(os.environ['raco'])
from common import get_widest_range_file, strip_dirname,\
        get_iters_from_file, get_dict, rsun
from get_parameter import get_parameter
from compute_grid_info import compute_theta_grid
from time_scales import compute_Prot, compute_tdt
from translate_times import translate_times

# Get the run directory on which to perform the analysis
dirname = sys.argv[1]
dirname_stripped = strip_dirname(dirname)

# Directory with data and plots, make the plotting directory if it doesn't
# already exist    
datadir = dirname + '/data/'
plotdir = dirname + '/plots/eflux_radial_vs_time/'
if (not os.path.isdir(plotdir)):
    os.makedirs(plotdir)

# Get list of shell slices at all possible times to plot
file_list, int_file_list, nfiles = get_file_lists(radatadir)

# Get command-line arguments to adjust the interval of averaging files
minmax = None
rnorm = None
rvals = None
plot_enth_fluc = False
mark_bcz = False

args = sys.argv[2:]
nargs = len(args)

# get the time range to make plots
the_tuple = get_desired_range(int_file_list, args)
if the_tuple is None: # By default plot the last 10 Shell_Slices
    index_first, index_last = nfiles - 11, nfiles - 1  
else:
    index_first, index_last = the_tuple

for i in range(nargs):
    arg = args[i]
    if arg == '-minmax':
        minmax = float(args[i+1]), float(args[i+2])
    elif arg == '-rnorm':
        rnorm = float(args[i+1])
    elif arg == '-fluc':
        plot_enth_fluc = True
    elif arg == '-rvals':
        rvals_str = args[i+1].split()
        rvals = []
        for rval_str in rvals_str:
            rvals.append(float(rval_str))
    elif arg == '-bcz': # try to estimate the real BCZ (and mark it)
                        # from where the enthalpy flux first goes negative
        mark_bcz = True

#Create the plot
lw = 1. # regular lines
#lw = 1.5 # Bit thicker lines

# Read in the flux data
print ('Getting radial fluxes from ' + datadir + Shell_Avgs_file + ' ...')
di = get_dict(datadir + Shell_Avgs_file)
vals = di['vals']
lut = di['lut']
nq = di['nq']
iter1, iter2 = di['iter1'], di['iter2']
rr = di['rr']
nr = di['nr']

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

# Determine the simulation is magnetic
magnetism = get_parameter(dirname, 'magnetism')

# Make the plot name, labelling the first/last iterations we average over
savename = dirname_stripped + '_eflux_radial_' +\
    str(iter1).zfill(8) + '_' + str(iter2).zfill(8) + '.png'

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

if plot_enth_fluc:
    if lut[1458] < nq:
        eflux_fluc = vals[ :, lut[1458]]
        eflux_mean = eflux_enth - eflux_fluc
    else: # do the Reynolds decomposition "by hand"
        # Compute the enthalpy flux from mean flows (MER. CIRC.)
        ref = ReferenceState(dirname + '/reference', '')
        nt = get_parameter(dirname, 'n_theta')
        tt, tw = compute_theta_grid(nt)
        tw_2d = tw.reshape((nt, 1))

        AZ_Avgs_file = get_widest_range_file(datadir, 'AZ_Avgs')
        rho = (ref.density).reshape((1, nr))
        ref_prs = (ref.pressure).reshape((1, nr))
        ref_temp = (ref.temperature).reshape((1, nr))
        prs_spec_heat = get_parameter(dirname, 'pressure_specific_heat')
        gamma = 5./3.

        di_az = get_dict(datadir + AZ_Avgs_file)
        vals_az = di_az['vals']
        lut_az = di_az['lut']

        vr_av = vals_az[:, :, lut_az[1]]
        entropy_av = vals_az[:, :, lut_az[501]]
        prs_av = vals_az[:, :, lut_az[502]]

        # Calculate mean temp. from EOS
        temp_av = ref_temp*((1.-1./gamma)*(prs_av/ref_prs) +\
                entropy_av/prs_spec_heat)

        # And, finally, the enthalpy flux from mean/fluc flows
        eflux_mean_az = rho*prs_spec_heat*vr_av*temp_av
        
        eflux_mean = np.sum(eflux_mean_az*tw_2d, axis=0)
        eflux_fluc = eflux - eflux_mean

if magnetism:
    qindex_mflux = lut[2001] # this is actually (-4*pi) TIMES 
                        # the correct Poynting flux
    mflux = -vals[:, qindex_mflux]/(4*np.pi)
    tflux += mflux
    
# Compute the integrated fluxes
fpr = 4*np.pi*rr**2
hflux_int = hflux*fpr
eflux_int = eflux*fpr
if plot_enth_fluc:
    eflux_mean_int = eflux_mean*fpr
    eflux_fluc_int = eflux_fluc*fpr
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

# Make lstar = lsun unless otherwise specified in main_input
try:
    lstar = get_parameter(dirname, 'luminosity')
except:
    lstar = 3.845e33

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
if plot_enth_fluc:
    plt.plot(rr_n, eflux_fluc_int/lstar, 'm--',\
            label=r'$\rm{F}_{enth,\ pp}$', linewidth=lw)
    plt.plot(rr_n, eflux_mean_int/lstar, 'm:',\
            label=r'$\rm{F}_{enth,\ mm}$', linewidth=lw)

# Mark zero line
plt.plot(rr_n, np.zeros_like(rr_n), 'k--', linewidth=lw)

# Get the y-axis in scientific notation
plt.ticklabel_format(useMathText=True, axis='y', scilimits=(0,0))

# Get ticks everywhere
plt.minorticks_on()
plt.tick_params(top=True, right=True, direction='in', which='both')

# Set the x limits
xmin, xmax = np.min(rr_n), np.max(rr_n)
delta_x = xmax - xmin
plt.xlim(xmin, xmax)

# Set the y-limits (the following values seem to "work well" for my models
# so far...perhaps adjust this in the future. 

if minmax is None:
    minmax = -0.7, 1.3
plt.ylim(minmax[0], minmax[1])

# Label the axes
if rnorm is None:
    plt.xlabel(r'$r/R_\odot$',fontsize=12, **csfont)
else:
    plt.xlabel(r'r/(%.1e cm)' %rnorm, fontsize=12, **csfont)

# Try to find the BCZ from where enthalpy flux goes negative, if desired
# avoid the outer boundary
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

# Label trace interval
if rotation:
    time_string = ('t = %.1f to %.1f ' %(t1/time_unit, t2/time_unit))\
            + time_label + (r'$\ (\Delta t = %.1f\ $'\
            %((t2 - t1)/time_unit)) + time_label + ')'
else:
    time_string = ('t = %.3f to %.3f ' %(t1/time_unit, t2/time_unit))\
            + time_label + (r'$\ (\Delta t = %.3f\ $'\
            %((t2 - t1)/time_unit)) + time_label + ')'

# Make title
the_title = dirname_stripped + '\n' + 'radial energy flux, ' + time_string
if mark_bcz:
    the_title += ('\n' + r'$r_{BCZ}/\rm{rnorm} = %.3f$' %rbcz)

plt.title(the_title, **csfont)

# Create a see-through legend
plt.legend(loc='lower left', shadow=True, ncol=3, fontsize=10)

# Last command
plt.tight_layout()

# Save the plot
print ('Saving the eflux plot at ' + plotdir + savename + ' ...')
plt.savefig(plotdir + savename, dpi=300)

# Show the plot
plt.show()