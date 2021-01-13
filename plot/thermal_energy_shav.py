###############################################
# Author: Loren Matilsky
# Date created: 11/05/2019
#
# This script plots the volume heating terms as functions of
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
from rayleigh_diagnostics import GridInfo
from read_eq_vp import read_eq_vp

# Get the run directory on which to perform the analysis
dirname = sys.argv[1]
dirname_stripped = strip_dirname(dirname)

# Directory with data and plots, make the plotting directory if it doesn't
# already exist    
datadir = dirname + '/data/'
plotdir = dirname + '/plots/'
if not os.path.isdir(plotdir):
    os.makedirs(plotdir)

# Find the Shell_Avgs file(s) in the data directory. If there are multiple, by
# default choose the one with widest range in the average
Shell_Avgs_file = get_widest_range_file(datadir, 'Shell_Avgs')

# Get command-line arguments to adjust the interval of averaging files
xminmax = None
minmax = None
rnorm = None
rvals = []
entropy_equation = False
force_econs = False # By default, no tachocline econs term to worry about
sep_czrz = False # plots boundary line in between cz/rz and computes
        # heating separately in both domains

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
    elif arg == '-rvals':
        rvals_str = args[i+1].split()
        rvals = []
        for rval_str in rvals_str:
            rvals.append(float(rval_str))
        rvals = np.array(rvals)
    elif arg == '-s':
        entropy_equation = True
    elif arg == '-econs':
        force_econs = True
    elif arg == '-czrz':
        sep_czrz = True

lw = 1. # regular lines
#lw = 1.5 # Bit thicker lines

# Read in the flux data
print ('Getting heating terms from ' + datadir + Shell_Avgs_file + ' ...')
di = get_dict(datadir + Shell_Avgs_file)
vals = di['vals']
lut = di['lut']
nq = di['nq']
iter1, iter2 = di['iter1'], di['iter2']
rr = di['rr']
ro = di['ro']
ri = di['ri']
nr = di['nr']

# Get the rho*T
eq = get_eq(dirname)
rhot = eq.density*eq.temperature

# Determine the simulation is magnetic
magnetism = get_parameter(dirname, 'magnetism')

# Make the plot name, labelling the first/last iterations we average over
if entropy_equation:
    basename = '_entropy_equation_shav_'
else:
    basename = '_thermal_energy_shav_'

savename = dirname_stripped + basename +\
    str(iter1).zfill(8) + '_' + str(iter2).zfill(8) + '.png'

advec_tot = -vals[:, lut[1401]]
advec_fluc = -vals[:, lut[1402]]
advec_mean = advec_tot - advec_fluc
advec_vr = -vals[:, lut[1406]]
cond_heating = vals[:, lut[1421]]
int_heating = vals[:, lut[1434]]
visc_heating = vals[:, lut[1435]]*rhot
tot_heating = advec_tot + cond_heating + int_heating + visc_heating
if magnetism:
    joule_heating = vals[:, lut[1436]]*rhot
    have_joule_fluc = False
    try:
        joule_heating_fluc = vals[:, lut[1437]]*rhot
        joule_heating_mean = joule_heating - joule_heating_fluc
        have_joule_fluc = True
        print ("Joule heating fluc (1437) output, so plotting ")
        print ("Reynolds decomposition of Joule heating.")
    except:
        print ("Joule heating fluc (1437) not output, ")
        print ("only plotting total Joule heating.")
    tot_heating += joule_heating

if force_econs:
    AZ_Avgs_file = get_widest_range_file(datadir, 'AZ_Avgs')
    print ('Getting forcing info from ' + datadir + AZ_Avgs_file)
    di_az = get_dict(datadir + AZ_Avgs_file)
    nt = di_az['nt']
    vals_az = di_az['vals']
    lut_az = di_az['lut']
    mean_vp = vals_az[:, :, lut_az[3]]
    vp2 = vals_az[:, :, lut_az[416]]
    fluc_vp2 = vals_az[:, :, lut_az[424]]
    mean_vp2 = vp2 - fluc_vp2
    tacho_r = get_parameter(dirname, 'tacho_r')
    print ("read tacho_r = %1.2e" %tacho_r)
    tacho_dr = get_parameter(dirname, 'tacho_dr')
    tacho_tau = get_parameter(dirname, 'tacho_tau')
    force_heating_az = np.zeros((nt, nr))
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
                    force_heating_az[it, ir] = rho[ir]*(mean_vp2[it, ir] -\
                            desired_vp*mean_vp[it, ir])/tacho_tau
                else:
                    force_heating_az[it, ir] = 0.
    else:
        forcing_coeff = rho/tacho_tau*0.5*(1.0 - np.tanh((rr - tacho_r)/(tacho_dr*rr[0])))
        force_heating_az = forcing_coeff.reshape((1, nr))*mean_vp2

    # Calculate the latitudinally integrated forcing work
    gi = GridInfo(dirname + '/grid_info')
    tw = gi.tweights
    force_heating = np.sum(force_heating_az*tw.reshape((nt, 1)), axis=0)
    tot_heating += force_heating

# Compute the INTEGRATED total heating
gi = GridInfo(dirname + '/grid_info', '')
rw = gi.rweights
shell_volume = 4.0*np.pi/3.0*(ro**3.0 - ri**3.0)
tot_heating_integrated = shell_volume*np.sum(tot_heating*rw)

if sep_czrz:
    nr_cz = get_parameter(dirname, 'ncheby')[1]
    nr_rz = nr - nr_cz
    rbcz = rr[nr_cz-1]

    # get averaging weights for CZ and RZ separately
    rw_cz = np.copy(rw[:nr_cz])
    rw_rz = np.copy(rw[nr_cz:])
    rw_cz /= np.sum(rw_cz)
    rw_rz /= np.sum(rw_rz)

    shell_volume_cz = 4.0*np.pi/3.0*(ro**3.0 - rbcz**3.0)
    shell_volume_rz = 4.0*np.pi/3.0*(rbcz**3.0 - ri**3.0)
    tot_heating_integrated_cz =\
            shell_volume*np.sum(tot_heating[:nr_cz]*rw_cz)
    tot_heating_integrated_rz =\
            shell_volume*np.sum(tot_heating[nr_cz:]*rw_rz)


if entropy_equation:
    advec_tot /= rhot
    advec_fluc /= rhot
    advec_mean /= rhot
    advec_vr /= rhot
    cond_heating /= rhot
    int_heating /= rhot
    visc_heating /= rhot
    tot_heating /= rhot
    if magnetism:
        joule_heating /= rhot
        if have_joule_fluc:
            joule_heating_fluc /= rhot
            joule_heating_mean /= rhot

    if force_econs:
        force_heating /= rhot

# Create the plot

# User can specify what to normalize the radius by
# By default, normalize by the solar radius
if rnorm is None:
    rr_n = rr/rsun
else:
    rr_n = rr/rnorm                                           

plt.plot(rr_n, advec_tot, 'm', label='advection tot', linewidth=lw)
plt.plot(rr_n, advec_fluc, 'm--', label='advection fluc', linewidth=lw)
plt.plot(rr_n, advec_mean, 'm:', label='advection mean', linewidth=lw)
plt.plot(rr_n, cond_heating, 'r', label='conductive heating', linewidth=lw)
plt.plot(rr_n, int_heating, 'g', label='internal heating', linewidth=lw)
plt.plot(rr_n, visc_heating, 'c', label='viscous heating', linewidth=lw)
if magnetism:
    plt.plot(rr_n, joule_heating, 'b', label='Joule heating tot',\
            linewidth=lw)
    if have_joule_fluc:
        plt.plot(rr_n, joule_heating_fluc, 'b--',\
                label='Joule heating fluc', linewidth=lw)
        plt.plot(rr_n, joule_heating_mean, 'b:',\
                label='Joule heating mean', linewidth=lw)
if force_econs:
    plt.plot(rr_n, force_heating, 'y',\
        label='Heating from forcing (econs)', linewidth=lw)
plt.plot(rr_n, tot_heating, 'k', label='total heating')

# Get the y-axis in scientific notation
plt.ticklabel_format(useMathText=True, axis='y', scilimits=(0,0))

# Get ticks everywhere
plt.minorticks_on()
plt.tick_params(top=True, right=True, direction='in', which='both')

# Set the x limits
if xminmax is None:
    xminmax = np.min(rr_n), np.max(rr_n)
delta_x = xminmax[1] - xminmax[0]
plt.xlim(xminmax[0], xminmax[1])

# Set the y-limits

if minmax is None:
    ymin = min(np.min(advec_tot), np.min(advec_fluc),\
            np.min(advec_mean), np.min(cond_heating), np.min(int_heating),\
            np.min(visc_heating), np.min(tot_heating)),
    ymax = max(np.max(advec_tot), np.max(advec_fluc),\
            np.max(advec_mean), np.max(cond_heating), np.max(int_heating),\
            np.max(visc_heating), np.max(tot_heating))
    if magnetism:
        ymin = min(ymin, np.min(joule_heating))
        ymax = max(ymax, np.max(joule_heating))
        if have_joule_fluc:
            ymin = min(ymin, np.min(joule_heating_fluc),\
                    np.min(joule_heating_mean))
            ymax = max(ymax, np.max(joule_heating_fluc),\
                    np.max(joule_heating_mean))
    delta_y = ymax - ymin
    ybuffer = 0.1*delta_y
    minmax = ymin - 3*ybuffer, ymax + ybuffer
plt.ylim(minmax[0], minmax[1])

# Label the axes
if rnorm is None:
    plt.xlabel(r'$r/R_\odot$',fontsize=12, **csfont)
else:
    plt.xlabel(r'r/(%.1e cm)' %rnorm, fontsize=12, **csfont)

# Mark radii if desired
if sep_czrz:
    if rvals is None:
        rvals = np.array([rbcz])
    else:
        rvals = np.hstack((rvals, np.array([rbcz])))

if not rvals is None:
    yvals = np.linspace(minmax[0], minmax[1], 100)
    for rval in rvals:
        if rnorm is None:
            rval_n = rval/rsun
        else:
            rval_n = rval/rnorm
#        plt.ylim(ymin, ymax)
        plt.plot(rval_n + np.zeros(100), yvals, 'k--')

if entropy_equation:
    plt.ylabel(r'$\partial S/\partial t\ \rm{(erg\ g^{-1}\ K^{-1}\ s^{-1})}$', fontsize=12, **csfont)
else:
    plt.ylabel('heating (' + r'$\rm{erg\ cm^{-3}\ s^{-1}}$' + ')',\
            fontsize=12, **csfont)

# Make title
if entropy_equation:
    basetitle = 'entropy eqn., '
else:
    basetitle = 'thermal energy eqn., ' 

lum = 3.846e33

title = dirname_stripped + '\n' + basetitle +\
          str(iter1).zfill(8) + ' to ' + str(iter2).zfill(8) +\
          ('\nInteg. Tot. Heating = %9.3e lsun'\
          %(tot_heating_integrated/lum))
if sep_czrz:
    title += ('\nInteg. Tot. Heating (CZ) = %9.3e lsun'\
            %(tot_heating_integrated_cz/lum)) +\
        ('\nInteg. Tot. Heating (RZ) = %9.3e lsun'\
        %(tot_heating_integrated_rz/lum))
plt.title(title, **csfont)

# Create a see-through legend
plt.legend(loc='lower left', shadow=True, ncol=2, fontsize=8)

# Last command
plt.tight_layout()

# Save the plot
print ('Saving the energy eqn. plot at ' + plotdir + savename)
plt.savefig(plotdir + savename, dpi=300)

# Show the plot
plt.show()
