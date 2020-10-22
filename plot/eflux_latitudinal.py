###############################################
# Author: Loren Matilsky
# Date created: 11/18
# Last modified: 11/18/2018
#
# This script computes the "cone"-averages latitudinal energy flux in the 
# spherical domain as a function of radius. 
# Plots various average energy fluxes, integrated over cones of opening angle theta
# Since Rayleigh has no "cone"-averages, we assume symmetry in phi and use the AZ_Avgs

import matplotlib as mpl
mpl.use('TkAgg')
import matplotlib.pyplot as plt
plt.rcParams['mathtext.fontset'] = 'dejavuserif'
csfont = {'fontname':'DejaVu Serif'}
import numpy as np
import sys, os
sys.path.append(os.environ['rapp'])
sys.path.append(os.environ['raco'])
from common import get_widest_range_file, strip_dirname, get_dict
from get_parameter import get_parameter
from rayleigh_diagnostics import GridInfo
from get_eq import get_eq
from time_scales import compute_Prot, compute_tdt
from translate_times import translate_times

# Get the run directory on which to perform the analysis
dirname = sys.argv[1]
dirname_stripped = strip_dirname(dirname)

# Directory with data and plots, make the plotting directory if it doesn't
# already exist    
datadir = dirname + '/data/'
plotdir = dirname + '/plots/'
if (not os.path.isdir(plotdir)):
    os.makedirs(plotdir)

# Find the Shell_Avgs file(s) in the data directory. If there are multiple, by
# default choose the one with widest range in the average
AZ_Avgs_file = get_widest_range_file(datadir, 'AZ_Avgs')
minmax = None

# Get command-line arguments to adjust the interval of averaging files
args = sys.argv[2:]
nargs = len(args)
for i in range(nargs):
    arg = args[i]
    if arg == '-usefile':
        AZ_Avgs_file = args[i+1]
        AZ_Avgs_file = AZ_Avgs_file.split('/')[-1]
    elif arg == '-minmax':
        minmax = float(args[i+1]), float(args[i+2])
       
print ('Getting data from ', datadir + AZ_Avgs_file, ' ...')
di = get_dict(datadir + AZ_Avgs_file)
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

# Make the plot name, labelling the first/last iterations we average over
savename = dirname_stripped + '_eflux_latitudinal_' + str(iter1).zfill(8) + '_' + str(iter2).zfill(8) + '.png'
print ("Will save plot at ", savename)

# Get grid info
gi = GridInfo(dirname + '/grid_info')
nr, nt = gi.nr, gi.ntheta
rr, tt, cost, sint = gi.radius, gi.theta, gi.costheta, gi.sintheta
ri, ro = np.min(rr), np.max(rr)
shell_depth = ro - ri
rweights = gi.rweights

#Create the plot
lw = 1.5 # Bit thicker lines

# Read in the flux data
#print ('Getting AZ_Avgs data from %s ...' %AZ_Avgs_file)
vals = di['vals']
lut = di['lut']
qindex_eflux = lut[1456]
qindex_eflux_pp = lut[1459]
qindex_cflux = lut[1471]
qindex_kflux = lut[1924]
qindex_vflux = lut[1936] # needs minus sign?

# Get the fluxes in the whole meridional plane
eflux = vals[:, :, qindex_eflux]
cflux = vals[:, :, qindex_cflux]
kflux = vals[:, :, qindex_kflux]
vflux = -vals[:, :, qindex_vflux]
tflux = eflux + cflux + kflux + vflux # compute the total flux

eq = get_eq(dirname)
try:
    eflux_pp = vals[:, :, qindex_eflux_pp]
except:
    # the fluctuating enthalpy flux wasn't computed directly
    # so compute it by subtracting out the flux due to meridional 
    # circulation
    rho_bar = eq.density
    T_bar = eq.temperature
    S_azav = vals[:, :, lut[501]]
    P_azav = vals[:, :, lut[502]]
    vt_azav = vals[:, :, lut[2]]
    eflux_mm = (rho_bar*T_bar*S_azav + P_azav)*vt_azav
    eflux_pp = eflux - eflux_mm

# Compute eflux_mm if the "try" succeeded
eflux_mm = eflux - eflux_pp

magnetism = get_parameter(dirname, 'magnetism')
if magnetism:
    qindex_mflux = lut[2002] # needs multiplication by -1/(4pi)?
    mflux = -1/(4*np.pi)*vals[:, :, qindex_mflux]
    tflux += mflux

# Compute the integrated fluxes
# At each point in the meridional plane we associate a "ring" of width dr and circumference 2 pi r sin(theta)
dr = rweights/rr**2/np.sum(rweights/rr**2)*shell_depth
areas = 2*np.pi*sint.reshape((nt, 1))*rr.reshape((1, nr))*\
        dr.reshape((1, nr))

eflux_int = np.sum(eflux*areas, axis=1)
eflux_pp_int = np.sum(eflux_pp*areas, axis=1)
eflux_mm_int = np.sum(eflux_mm*areas, axis=1)

cflux_int = np.sum(cflux*areas, axis=1)
kflux_int = np.sum(kflux*areas, axis=1)
vflux_int = np.sum(vflux*areas, axis=1)
if magnetism:
    mflux_int = np.sum(mflux*areas, axis=1)
tflux_int = np.sum(tflux*areas, axis=1)

# compute the "equilibrium flux" (latitudinal flux needed to balance out
# any differences between the inner and outer radial fluxes
cfluxr = vals[:, :, lut[1470]]
cflux_out = cfluxr[:, 0]
lum = eq.lum
rsq_Flux_in = lum/4/np.pi
integrand = -2*np.pi*(ro**2*cflux_out - rsq_Flux_in*np.ones(nt))
    
eqflux_int = np.zeros(nt)

# Get the latitudinal integration weights
gi = GridInfo(dirname + '/grid_info')
tw = gi.tweights
for it in range(nt):
    # Remember the variables are index "backwards" w.r.t. it (theta
    # runs from pi to 0)
    #eqflux_int[it] = 2*np.sum(tw[it:]*integrand[it:])
    if it <= nt//2:
        eqflux_int[it] = 2*np.sum(tw[it:nt//2]*integrand[it:nt//2])
    else:
        eqflux_int[it] = -2*np.sum(tw[nt//2:it]*integrand[nt//2:it])

# Create the plot; start with plotting all the energy fluxes
lats = 180*(np.pi/2 - tt)/np.pi
plt.plot(lats, eflux_int/lum, 'm', label=r'$\rm{F}_{enth}$',\
        linewidth=lw)
plt.plot(lats, eflux_pp_int/lum, 'm--', label=r'$\rm{F}_{enth,\ pp}$',\
        linewidth=lw)
plt.plot(lats, eflux_mm_int/lum, 'm:', label=r'$\rm{F}_{enth,\ mm}$',\
        linewidth=lw)
plt.plot(lats, cflux_int/lum, label=r'$\rm{F}_{cond}$', linewidth=lw)
plt.plot(lats, kflux_int/lum, label=r'$\rm{F}_{KE}$', linewidth=lw)
plt.plot(lats, vflux_int/lum, label=r'$\rm{F}_{visc}$', linewidth=lw)
if magnetism:
    plt.plot(lats, mflux_int/lum, label=r'$\rm{F}_{Poynting}$',\
            linewidth=lw)
plt.plot(lats, tflux_int/lum, label=r'$\rm{F}_{total}$',\
        linewidth=lw, color='black')
plt.plot(lats, eqflux_int/lum, 'k--', label=r'$\rm{F}_{eq}$',\
        linewidth=lw)

# Get the y-axis in scientific notation
#plt.ticklabel_format(useMathText=True, axis='y', scilimits=(0,0))

# Get ticks everywhere
plt.minorticks_on()
plt.tick_params(top=True, right=True, direction='in', which='both')

# Set the x limits
xmin, xmax = -90, 90
delta_x = xmax - xmin
plt.xlim(xmin, xmax)

# Set the y-limits 
maxabs = max((np.max(np.abs(eflux_int)), np.max(np.abs(cflux_int)),\
        np.max(np.abs(kflux_int)), np.max(np.abs(vflux_int)), np.max(np.abs(tflux_int))))

# Set y limits if user wanted you to
if not minmax is None:
    plt.ylim(minmax[0], minmax[1])
else:
    ymin, ymax = -1.2*maxabs/lum, 1.2*maxabs/lum
    delta_y = ymax - ymin
    plt.ylim(ymin, ymax)

# Label the axes
plt.xlabel(r'$\rm{Latitude\ (deg)}$', fontsize=12)
plt.ylabel('(Integrated Energy Flux)' + r'$/L_\odot$',\
        fontsize=12)

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
the_title = dirname_stripped + '\n' + 'latitudinal energy flux, ' +\
        time_string
plt.title(the_title, **csfont)

# Create a see-through legend
leg=plt.legend(loc='lower left',shadow=True, ncol=3,fontsize=10)
leg.get_frame().set_alpha(0.3)

# Last command
plt.tight_layout()

# Save the plot
print ('Saving the eflux plot at ' + plotdir + savename + ' ...')
plt.savefig(plotdir + savename, dpi=300)

# Show the plot
plt.show()
