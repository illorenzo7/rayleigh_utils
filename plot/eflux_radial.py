###############################################
# Author: Loren Matilsky
# Date created: 02/14/2018
#
# This script plots the radial energy fluxes as functions of
# radius using from the Shell_Avgs data

import matplotlib as mpl
mpl.use('TkAgg')
import matplotlib.pyplot as plt
plt.rcParams['mathtext.fontset'] = 'dejavuserif'
csfont = {'fontname':'DejaVu Serif'}
import numpy as np
import sys, os
sys.path.append(os.environ['rapp'])
sys.path.append(os.environ['raco'])
from common import *
from compute_grid_info import compute_theta_grid
from read_inner_vp import read_inner_vp
from read_eq_vp import read_eq_vp
from rayleigh_diagnostics import GridInfo
from reference_tools import equation_coefficients

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
plotdir = dirname + '/plots/'
if (not os.path.isdir(plotdir)):
    os.makedirs(plotdir)

# Find the Shell_Avgs file(s) in the data directory. If there are multiple, by
# default choose the one with widest range in the average
Shell_Avgs_file = get_widest_range_file(datadir, 'Shell_Avgs')
forced = False

# Get command-line arguments to adjust the interval of averaging files
minmax = None
xminmax = None
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
for i in range(nargs):
    arg = args[i]
    if arg == '-usefile':
        Shell_Avgs_file = args[i+1]
        Shell_Avgs_file = Shell_Avgs_file.split('/')[-1]
    elif arg == '-minmax':
        minmax = float(args[i+1]), float(args[i+2])
    elif arg == '-xminmax':
        xminmax = float(args[i+1]), float(args[i+2])
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
    elif arg == '-forced':
        forced = True
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

#Create the plot

# Read in the flux data
print ('Getting radial fluxes from ' + datadir + Shell_Avgs_file + ' ...')
di = get_dict(datadir + Shell_Avgs_file)
vals = di['vals']
lut = di['lut']
nq = di['nq']
iter1, iter2 = di['iter1'], di['iter2']
rr = di['rr']
nr = di['nr']

if forced:
    AZ_Avgs_file = get_widest_range_file(datadir, 'AZ_Avgs')
    print ('Getting forcing info from ' + datadir + AZ_Avgs_file + ' ...')
    di_az = get_dict(datadir + AZ_Avgs_file)
    vals_az = di_az['vals']
    lut_az = di_az['lut']

    gi = GridInfo(dirname + '/grid_info')
    rw = gi.rweights
    tw = gi.tweights
    nt = gi.ntheta
    ro = np.max(rr)

    mean_vp = vals_az[:, :, lut_az[3]]
    vp2 = vals_az[:, :, lut_az[416]]
    fluc_vp2 = vals_az[:, :, lut_az[424]]
    mean_vp2 = vp2 - fluc_vp2
    tacho_r = get_parameter(dirname, 'tacho_r')
    print ("read tacho_r = %1.2e" %tacho_r)
    tacho_dr = get_parameter(dirname, 'tacho_dr')
    tacho_tau = get_parameter(dirname, 'tacho_tau')
    work_forcing = np.zeros((nt, nr))
    eq = get_eq(dirname)
    rho = eq.rho

    if os.path.exists(dirname + '/eq_vp'): 
        print ("eq_vp file exists, so I assume you have a forcing function which\n quartically matches on to a CZ differential rotation\n with viscous-torque-free buffer zone")
        tacho_r2 = get_parameter(dirname, 'tacho_r2')
        i_tacho_r = np.argmin(np.abs(rr - tacho_r))
        print ("read tacho_r2 = %1.2e" %tacho_r2)
        eq_vp = read_eq_vp(dirname + '/eq_vp', nt, nr)
        for it in range(nt):
            for ir in range(nr):
                if rr[ir] <= tacho_r2:
                    if rr[ir] > tacho_r:
                        # Here is where the DR is forced differentially
                        # (a "buffer zone" to reduce shear)
                        desired_vp = eq_vp[it, ir]
                    elif rr[ir] > tacho_r - tacho_dr*rr[0]:
                        # Here is where the DR is forced to match 
                        # quartically from differential to solid-body
                        desired_vp = eq_vp[it, i_tacho_r]*(1.0 - ( (rr[ir] - tacho_r)/(tacho_dr*rr[0]) )**2)**2
                    else:
                        desired_vp = 0.0
                    work_forcing[it, ir] = -rho[ir]*(mean_vp2[it, ir] -\
                            desired_vp*mean_vp[it, ir])/tacho_tau
                else:
                    work_forcing[it, ir] = 0.
    else:
        forcing_coeff = -rho/tacho_tau*0.5*(1.0 - np.tanh((rr - tacho_r)/(tacho_dr*rr[0])))
        work_forcing = forcing_coeff.reshape((1, nr))*mean_vp2

    work_forcing_shav = np.sum(work_forcing*tw.reshape((nt, 1)), axis=0)
    forcing_flux_int = np.zeros(nr)
    for ir in range(nr):
        work_forcing_gav = np.sum((work_forcing_shav*rw)[:ir+1])/\
                np.sum(rw[:ir+1])
        volume = 4./3.*np.pi*(ro**3. - rr[ir]**3.)
        forcing_flux_int[ir] = work_forcing_gav*volume

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

hflux = vals[:, lut[1433]]
eflux = vals[:, lut[1455]]
cflux = vals[:, lut[1470]]
kflux = vals[:, lut[1923]]
vflux = -vals[:, lut[1935]]
tflux = hflux + eflux + cflux + kflux + vflux # compute the total flux

if plot_enth_fluc:
    if lut[1458] < nq:
        eflux_fluc = vals[ :, lut[1458]]
        eflux_mean = eflux - eflux_fluc
    else: # do the Reynolds decomposition "by hand"
        # Compute the enthalpy flux from mean flows (MER. CIRC.)
        eq = get_eq(dirname)
        nt = get_parameter(dirname, 'n_theta')
        tt, tw = compute_theta_grid(nt)
        tw_2d = tw.reshape((nt, 1))

        AZ_Avgs_file = get_widest_range_file(datadir, 'AZ_Avgs')
        rho = (eq.density).reshape((1, nr))
        ref_prs = (eq.pressure).reshape((1, nr))
        ref_temp = (eq.temperature).reshape((1, nr))
        prs_spec_heat = 3.5e8
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
    # A Space Oddysey is actually (-4*pi) TIMES the correct Poynting flux
    mflux = -vals[:, lut[2001]]/(4*np.pi)
    tflux += mflux

fpr = 4*np.pi*rr**2
if forced:
    tflux += forcing_flux_int/fpr

# Compute the integrated fluxes
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

lstar = get_lum(dirname)

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
if forced:
    plt.plot(rr_n, forcing_flux_int/lstar, 'g--',\
            label=r'$\rm{F}_{forcing}$', linewidth=lw)

# Mark zero line
plt.plot(rr_n, np.zeros_like(rr_n), 'k--', linewidth=lw)

# Get the y-axis in scientific notation
plt.ticklabel_format(useMathText=True, axis='y', scilimits=(0,0))

# Get ticks everywhere
plt.minorticks_on()
plt.tick_params(top=True, right=True, direction='in', which='both')

# Set the x limits
if xminmax is None:
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

if minmax is None:
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
plt.savefig(plotdir + savename, dpi=dpi)

# Show the plot
plt.show()
