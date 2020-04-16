# Author: Loren Matilsky
# Date created: 03/31/2020
import matplotlib as mpl
mpl.use('TkAgg')
import matplotlib.pyplot as plt
import pickle
import numpy as np
import sys, os
sys.path.append(os.environ['raco'])
sys.path.append(os.environ['rapp'])
from subprocess import call
from common import get_file_lists, get_widest_range_file, strip_dirname, get_dict
from get_parameter import get_parameter
from time_scales import compute_Prot, compute_tdt

# Get the run directory on which to perform the analysis
dirname = sys.argv[1]

# Data and plot directories
datadir = dirname + '/data/'
plotdir = dirname + '/plots/'
if (not os.path.isdir(plotdir)):
    os.makedirs(plotdir)
dirname_stripped = strip_dirname(dirname)

# Find the etrace file(s) in the data directory. If there are multiple, by
# default choose the one with widest range in the trace.
the_file = get_widest_range_file(datadir, 'trace_2dom_G_Avgs')

# Set defaults
xiter = False
notfrom0 = False
magnetism = False
minmax = None
xminmax = None
plotall = False

# Get command-line arguments
args = sys.argv[2:]
nargs = len(args)
for i in range(nargs):
    arg = args[i]
    if arg == '-xiter': # plot w.r.t. iterations
        xiter = True
    elif arg == '-usefile':
        the_file = args[i+1]
        the_file = the_file.split('/')[-1]
    elif arg == '-notfrom0':
        notfrom0 = True
    elif arg == '-minmax':
        try: # first see if user wants to set all three ranges set
            # If "plotall", then use the first range for the total amom,
            # second for the convective, third for the mean
            minmax = float(args[i+1]), float(args[i+2]),\
                    float(args[i+3]), float(args[i+4]),\
                    float(args[i+5]), float(args[i+6])
        except: # if not, they want the same pair used for each subplot
            minmax = float(args[i+1]), float(args[i+2])
    elif arg == '-xminmax':
        xminmax = float(args[i+1]), float(args[i+2])
    elif arg == '-plotall':
        plotall = True

# Tag the plot by whether or not the x axis is in "time" or "iteration"
if xiter:
    tag = '_xiter'
else:
    tag = '_xtime'

# Read in the KE data (dictionary form)
di = get_dict(datadir + the_file)

vals_cz = di['vals_cz']
vals_rz = di['vals_rz']
lut = di['lut']
times = di['times']
iters = di['iters']
iter1 = di['iter1']
iter2 = di['iter2']

# Get the baseline time unit
rotation = get_parameter(dirname, 'rotation')
if rotation:
    time_unit = compute_Prot(dirname)
    time_label = r'$\rm{P_{rot}}$'
else:
    time_unit = compute_tdt(dirname)
    time_label = r'$\rm{TDT}$'

# Get the first and last indices in the time direction, depending on 
# xminmax
ntimes = len(times)
if xminmax is None:
    it1, it2 = 0, ntimes - 1 # Plot over whole range by default
else:
    if xiter:
        it1 = np.argmin(np.abs(iters - xminmax[0]))
        it2 = np.argmin(np.abs(iters - xminmax[1]))
    else:
        it1 = np.argmin(np.abs(times/time_unit - xminmax[0]))
        it2 = np.argmin(np.abs(times/time_unit - xminmax[1]))

# Make appropriate file name to save
savename = dirname_stripped + '_Ltrace_2dom_' + str(iter1).zfill(8) + '_' + str(iter2).zfill(8) + tag + '.png'

# Compute the volumes of the two shells
ri, rt, ro = get_parameter(dirname, 'domain_bounds')
Vcz = 4./3.*np.pi*(ro**3. - rt**3.)
Vrz = 4./3.*np.pi*(rt**3. - ri**3.)

# Get data in apropriate time range
times = times[it1:it2+1]
t1 = np.min(times)
t2 = np.max(times)
iters = iters[it1:it2+1]

# Compute the total angular momentum in the shell
Lz_cz = vals_cz[it1:it2+1, lut[1819]]*Vcz
Lz_rz = vals_rz[it1:it2+1, lut[1819]]*Vrz
if plotall:
    Lx_cz = vals_cz[it1:it2+1, lut[1820]]*Vcz
    Lx_rz = vals_rz[it1:it2+1, lut[1820]]*Vrz
    Ly_cz = vals_cz[it1:it2+1, lut[1821]]*Vcz
    Ly_rz = vals_rz[it1:it2+1, lut[1821]]*Vrz

Lpz_cz = vals_cz[it1:it2+1, lut[1822]]*Vcz
Lpz_rz = vals_rz[it1:it2+1, lut[1822]]*Vrz
if plotall:
    Lpx_cz = vals_cz[it1:it2+1, lut[1823]]*Vcz
    Lpx_rz = vals_rz[it1:it2+1, lut[1823]]*Vrz
    Lpy_cz = vals_cz[it1:it2+1, lut[1824]]*Vcz
    Lpy_rz = vals_rz[it1:it2+1, lut[1824]]*Vrz

Lmz_cz = vals_cz[it1:it2+1, lut[1825]]*Vcz
Lmz_rz = vals_rz[it1:it2+1, lut[1825]]*Vrz
if plotall:
    Lmx_cz = vals_cz[it1:it2+1, lut[1826]]*Vcz
    Lmx_rz = vals_rz[it1:it2+1, lut[1826]]*Vrz
    Lmy_cz = vals_cz[it1:it2+1, lut[1827]]*Vcz
    Lmy_rz = vals_rz[it1:it2+1, lut[1827]]*Vrz

# Get global min/max vals
mmin = min(np.min(Lz_cz), np.min(Lz_rz), np.min(Lz_cz + Lz_rz))
mmax = max(np.max(Lz_cz), np.max(Lz_rz), np.max(Lz_cz + Lz_rz))
mminm = min(np.min(Lmz_cz), np.min(Lmz_rz), np.min(Lmz_cz + Lmz_rz))
mmaxm = max(np.max(Lmz_cz), np.max(Lmz_rz), np.max(Lmz_cz + Lmz_rz))
mminp = min(np.min(Lpz_cz), np.min(Lpz_rz), np.min(Lpz_cz + Lpz_rz))
mmaxp = max(np.max(Lpz_cz), np.max(Lpz_rz), np.max(Lpz_cz + Lpz_rz))

print('mmin mmax: ', mmin, mmax)
print('mminm mmaxm: ', mminm, mmaxm)
print('mminp mmaxp: ', mminp, mmaxp)
if not xiter:
    xaxis = times/time_unit
else:
    xaxis = iters

# create figure with  3 panels in a row (total, mean and fluctuating amom)
fig, axs = plt.subplots(3, 1, figsize=(5, 10), sharex=True)
ax1 = axs[0]; ax2 = axs[1]; ax3 = axs[2]

# Make thin lines to see structure of variation
lw = 0.5

# Time string showing trace interval
if rotation:
    time_string = ('t = %.1f to %.1f ' %(t1/time_unit, t2/time_unit))\
            + time_label + (r'$\ (\Delta t = %.1f\ $'\
            %((t2 - t1)/time_unit)) + time_label + ')'
else:
    time_string = ('t = %.3f to %.3f ' %(t1/time_unit, t2/time_unit))\
            + time_label + (r'$\ (\Delta t = %.3f\ $'\
            %((t2 - t1)/time_unit)) + time_label + ')'

# Make first plot (total angular momentum trace)
ax1.plot(xaxis, Lz_cz, 'r', linewidth=lw, label=r'$\mathcal{L}_z\ \rm{(CZ)}$')
ax1.plot(xaxis, Lz_rz, 'b', linewidth=lw, label=r'$\mathcal{L}_z\ \rm{(RZ)}$')
ax1.plot(xaxis, Lz_cz + Lz_rz, 'k', linewidth=lw, label=r'$\mathcal{L}_z\ \rm{(tot)}$')
if plotall:
    ax1.plot(xaxis, Lx_cz, 'r--', linewidth=lw, label=r'$\mathcal{L}_x\ \rm{(CZ)}$')
    ax1.plot(xaxis, Lx_rz, 'b--', linewidth=lw, label=r'$\mathcal{L}_x\ \rm{(RZ)}$')
    ax1.plot(xaxis, Ly_cz, 'r:', linewidth=lw, label=r'$\mathcal{L}_y\ \rm{(CZ)}$')
    ax1.plot(xaxis, Ly_rz, 'b:', linewidth=lw, label=r'$\mathcal{L}_y\ \rm{(RZ)}$')

# Make title
ax1.set_title(dirname_stripped + '\n ' + time_string +\
          '\n\nangular momentum')

# Make the second plot (angular momentum of fluctuating motions)
ax2.plot(xaxis, Lpz_cz, 'r', linewidth=lw, label=r'$\mathcal{L}_z^\prime\ \rm{(CZ)}$')
ax2.plot(xaxis, Lpz_rz, 'b', linewidth=lw, label=r'$\mathcal{L}_z^\prime\ \rm{(RZ)}$')
ax2.plot(xaxis, Lpz_cz + Lpz_rz, 'k', linewidth=lw, label=r'$\mathcal{L}_z^\prime\ \rm{(tot)}$')
if plotall:
    ax2.plot(xaxis, Lpx_cz, 'r--', linewidth=lw, label=r'$\mathcal{L}_x^\prime\ \rm{(CZ)}$')
    ax2.plot(xaxis, Lpx_rz, 'b--', linewidth=lw, label=r'$\mathcal{L}_x^\prime\ \rm{(RZ)}$')
    ax2.plot(xaxis, Lpy_cz, 'r:', linewidth=lw, label=r'$\mathcal{L}_y^\prime\ \rm{(CZ)}$')
    ax2.plot(xaxis, Lpy_rz, 'b:', linewidth=lw, label=r'$\mathcal{L}_y^\prime\ \rm{(RZ)}$')

# Title and axis label
ax2.set_title('convection amom')
# Put the y-label on the middle plot
ax2.set_ylabel(r'$\rm{angular\ momentum\ density\ (g\ cm^{2}\ s^{-1})}$')

# Third plot: angular momentum of mean energies
ax3.plot(xaxis, Lmz_cz, 'r', linewidth=lw, label=r'$\langle\mathcal{L}\rangle_z\ \rm{(CZ)}$')
ax3.plot(xaxis, Lmz_rz, 'b', linewidth=lw, label=r'$\langle\mathcal{L}\rangle_z\ \rm{(RZ)}$')
ax3.plot(xaxis, Lmz_cz + Lmz_rz, 'k', linewidth=lw, label=r'$\langle\mathcal{L}_z\rangle\ \rm{(tot)}$')
if plotall:
    ax3.plot(xaxis, Lmx_cz, 'r--', linewidth=lw, label=r'$\langle\mathcal{L}\rangle_x\ \rm{(CZ)}$')
    ax3.plot(xaxis, Lmx_rz, 'b--', linewidth=lw, label=r'$\langle\mathcal{L}\rangle_x\ \rm{(RZ)}$')
    ax3.plot(xaxis, Lmy_cz, 'k', linewidth=lw, label=r'$\langle\mathcal{L}\rangle_y\ \rm{(CZ)}$')
    ax3.plot(xaxis, Lmy_rz, 'k--', linewidth=lw, label=r'$\langle\mathcal{L}\rangle_y\ \rm{(RZ)}$')

# title and x-axis label
ax3.set_title('mean-motion amom')
 
# Set x limits  
if xminmax is None:
    ax1.set_xlim((np.min(xaxis), np.max(xaxis)))
else:
    ax1.set_xlim(xminmax)

# Determine default minmax range
if minmax is None:
    diff = mmax - mmin
    buff = 0.05*diff
    diffm = mmaxm - mminm
    buffm = 0.05*diffm
    diffp = mmaxp - mminp
    buffp = 0.05*diffp
    ax1.set_ylim(mmin - buff, mmax + buff)
    ax2.set_ylim(mminp - buffp, mmaxp + buffp)
    ax3.set_ylim(mminm - buffm, mmaxm + buffm)
else:
    if len(minmax) == 2:
        ax1.set_ylim(minmax)
        ax2.set_ylim(minmax)
        ax3.set_ylim(minmax)
    elif len(minmax) == 6:
        ax1.set_ylim(minmax[0], minmax[1])
        ax1.set_ylim(minmax[2], minmax[3])
        ax1.set_ylim(minmax[4], minmax[5])

# legend
ax1.legend(ncol=2, fontsize=8)

# Make axes use scientific notation
ax1.ticklabel_format(scilimits = (-3,4), useMathText=True)
ax2.ticklabel_format(scilimits = (-3,4), useMathText=True)
ax3.ticklabel_format(scilimits = (-3,4), useMathText=True)


# Put the x-axis label on the bottom
if (xiter):
    ax3.set_xlabel('iteration #')
else:
    if rotation:
        ax3.set_xlabel(r'$\rm{t\ (P_{rot})}$')
    else:
        ax3.set_xlabel(r'$\rm{t\ (T_{diff})}$')

# Get ticks everywhere
plt.sca(ax1)
plt.minorticks_on()
plt.tick_params(top=True, right=True, direction='in', which='both')

plt.sca(ax2)
plt.minorticks_on()
plt.tick_params(top=True, right=True, direction='in', which='both')

plt.sca(ax3)
plt.minorticks_on()
plt.tick_params(top=True, right=True, direction='in', which='both')

# Space the subplots to make them look pretty
plt.tight_layout
plt.subplots_adjust(left=0.15, bottom=0.08, top=0.85, wspace=0.4)

# Save the plot
print ('Saving the etrace plot at ' + plotdir + savename + ' ...')
plt.savefig(plotdir + savename, dpi=300)

# Show the plot
plt.show()
