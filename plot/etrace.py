# Author: Loren Matilsky
# Date created: 06/08/2017
import matplotlib as mpl
mpl.use('TkAgg')
import matplotlib.pyplot as plt
import pickle
import numpy as np
import sys, os
sys.path.append(os.environ['raco'])
sys.path.append(os.environ['rapp'])
from subprocess import call
from common import get_file_lists, get_widest_range_file, strip_dirname,\
        get_dict
from get_parameter import get_parameter
from rayleigh_diagnostics import GridInfo
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
the_file = get_widest_range_file(datadir, 'trace_G_Avgs')

# Set defaults
xiter = False
from0 = False
magnetism = False
ylog = False
minmax = None
xminmax = None
xmin = None
xmax = None
plot_inte = False
plot_tote = False
savename = None

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
    elif arg == '-from0':
        from0 = True
    elif arg == '-mag':
        magnetism = True
    elif arg == '-log':
        ylog = True
    elif arg == '-minmax':
        minmax = float(args[i+1]), float(args[i+2])
    elif arg == '-xminmax':
        xminmax = float(args[i+1]), float(args[i+2])
    elif arg == '-xmin':
        xmin = float(args[i+1])
    elif arg == '-xmax':
        xmax = float(args[i+1])
    elif arg == '-inte':
        plot_inte = True
    elif arg == '-tote':
        plot_tote = True
    elif arg == '-name':
        savename = args[i+1] + '.png'

# Tag the plot by whether or not the x axis is in "time" or "iteration"
if (xiter):
    tag = '_xiter'
else:
    tag = '_xtime'

# Read in the KE data (dictionary form)
print ('Getting energy trace from ' + datadir + the_file + ' ...')
di = get_dict(datadir + the_file)

vals = di['vals']
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

# Make appropriate file name to save
basename = '_etrace_'
if plot_inte:
    basename += 'inte_'
if plot_tote:
    basename += 'tote_'

if savename is None:
    savename = dirname_stripped + basename + str(iter1).zfill(8) + '_' + str(iter2).zfill(8) + tag + '.png'

# Take slices based on what xminmax is
if not xiter:
    xaxis = times/time_unit
else:
    xaxis = iters

if from0:
    xmin = 0.
else:
    xmin = np.min(xaxis)

if xminmax is None:
    xminmax = xmin, np.max(xaxis)

# Change JUST xmin or xmax, if desired
if not xmin is None:
    xminmax = xmin, xminmax[1]
if not xmax is None:
    xminmax = xminmax[0], xmax

ixmin = np.argmin(np.abs(xaxis - xminmax[0]))
ix_max = np.argmin(np.abs(xaxis - xminmax[1]))
t1 = times[ixmin]
t2 = times[ix_max]

xaxis = xaxis[ixmin:ix_max + 1]
#ke = vals[lut[401]][ixmin:ix_max + 1]
rke = vals[lut[402]][ixmin:ix_max + 1]
tke = vals[lut[403]][ixmin:ix_max + 1]
pke = vals[lut[404]][ixmin:ix_max + 1]
ke = rke + tke + pke

#fke = vals[lut[409]][ixmin:ix_max + 1]
frke = vals[lut[410]][ixmin:ix_max + 1]
ftke = vals[lut[411]][ixmin:ix_max + 1]
fpke = vals[lut[412]][ixmin:ix_max + 1]
fke = frke + ftke + fpke

#mke = vals[lut[405]][ixmin:ix_max + 1]
#mrke = vals[lut[406]][ixmin:ix_max + 1]
#mtke = vals[lut[407]][ixmin:ix_max + 1]
#mpke = vals[lut[408]][ixmin:ix_max + 1]
mrke = rke - frke
mtke = tke - ftke
mpke = pke - fpke
mke = mrke + mtke + mpke

# Get the magnetic energies if they are available
if magnetism:
    #me = vals[lut[1101]][ixmin:ix_max + 1]
    rme = vals[lut[1102]][ixmin:ix_max + 1]
    tme = vals[lut[1103]][ixmin:ix_max + 1]
    pme = vals[lut[1104]][ixmin:ix_max + 1]
    me = rme + tme + pme

    #mme = vals[lut[1105]][ixmin:ix_max + 1]
    #mrme = vals[lut[1106]][ixmin:ix_max + 1]
    #mtme = vals[lut[1107]][ixmin:ix_max + 1]
    #mpme = vals[lut[1108]][ixmin:ix_max + 1]
    #mme = mrme + mtme + mpme

    #fme = vals[lut[1109]][ixmin:ix_max + 1]
    frme = vals[lut[1110]][ixmin:ix_max + 1]
    ftme = vals[lut[1111]][ixmin:ix_max + 1]
    fpme = vals[lut[1112]][ixmin:ix_max + 1]
    fme = frme + ftme + fpme

    mrme = rme - frme
    mtme = tme - ftme
    mpme = pme - fpme
    mme = mrme + mtme + mpme

if plot_inte or plot_tote:
    try: # get the internal energy if it is available
        try: # First try to get it from Shell_Avgs data
            the_file = get_widest_range_file(datadir,\
                    'inte_from_Shell_Avgs')
            di_inte = get_dict(datadir + the_file)
            inte = di_inte['inte'][ixmin:ix_max + 1]
            print("Got internal energy from Shell_Avgs")
        except:
            inte = vals[lut[701]][ixmin:ix_max + 1]
            print("Got internal energy from G_Avgs")
    except:
        print ("Internal energy not available; setting to 0")
        int_e = np.zeros_like(times[ixmin:ix_max + 1])

if plot_tote:
    tote = ke + inte
    if magnetism:
        tote += me

# Get global min/max vals
if magnetism:
    mmax = np.max((np.max(ke), np.max(me)))
    mmin = np.min((np.min(mrke), np.min(mtke), np.min(mpke), np.min(frke), np.min(ftke), np.min(fpke),\
     np.min(mrme), np.min(mtme), np.min(mpme), np.min(frke), np.min(ftme), np.min(fpme)))
else:
    mmax = np.max(ke)
    mmin = np.min((np.min(mrke), np.min(mtke), np.min(mpke), np.min(frke), np.min(ftke), np.min(fpke)))   

if plot_inte:
    mmax = max(np.max(inte), mmax)
    mmin = min(np.min(inte), mmin)

if plot_tote:
    mmax = max(np.max(tote), mmax)
    mmin = min(np.min(tote), mmin)

# create figure with  3 panels in a row (total, mean and fluctuating energy)
fig, axs = plt.subplots(3, 1, figsize=(5, 10), sharex=True, sharey=True)
ax1 = axs[0]; ax2 = axs[1]; ax3 = axs[2]

# Make thin lines to see structure of variation
lw = 0.5

# first plot: total kinetic energy trace      
ax1.plot(xaxis, ke, 'k', linewidth=lw, label=r'$\rm{KE_{tot}}$')
ax1.plot(xaxis, rke, 'r', linewidth=lw, label=r'$\rm{KE}_r$')
ax1.plot(xaxis, tke, 'g', linewidth=lw, label=r'$\rm{KE}_\theta$')
ax1.plot(xaxis, pke, 'b', linewidth=lw, label=r'$\rm{KE}_\phi$')

# If magnetic, plot magnetic energies!
if magnetism:
    ax1.plot(xaxis, me, 'k--', linewidth=lw, label=r'$\rm{ME_{tot}}$')
    ax1.plot(xaxis, rme, 'r--', linewidth=lw, label=r'$\rm{ME}_r$')
    ax1.plot(xaxis, tme, 'g--', linewidth=lw, label=r'$\rm{ME}_\theta$')
    ax1.plot(xaxis, pme, 'b--', linewidth=lw, label=r'$\rm{ME}_\phi$')


if plot_inte:
    ax1.plot(xaxis, inte, 'm', linewidth=lw, label='INT E')

if plot_tote:
    ax1.plot(xaxis, tote, 'c', linewidth=lw, label='TOT E')

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
title = dirname_stripped + '\n ' + time_string +\
          '\ntotal energy'
if plot_inte:
    # Compute change in energy over time, to add to label
    gi = GridInfo(dirname + '/grid_info')
    ri, ro = np.min(gi.radius), np.max(gi.radius)
    shell_volume = 4/3*np.pi*(ro**3 - ri**3)
    dE = (inte[-1] - inte[0])*shell_volume
    dt = times[-1] - times[0]
    title += (('\nINT E: ' + r'$\rm{\Delta E/\Delta t = %1.3e\ cgs}$')\
            %(dE/dt))
if plot_tote:
    gi = GridInfo(dirname + '/grid_info')
    ri, ro = np.min(gi.radius), np.max(gi.radius)
    shell_volume = 4/3*np.pi*(ro**3 - ri**3)
    dE = (tote[-1] - tote[0])*shell_volume
    dt = times[-1] - times[0]
    title += (('\nTOT E: ' + r'$\rm{\Delta E/\Delta t = %1.3e\ cgs}$')\
            %(dE/dt))

ax1.set_title(title)

if minmax is None:
    if ylog:
        minmax = mmin/3.0, mmax*3.0
    else:
        ydiff = mmax - mmin
        ybuffer = 0.05*ydiff 
        minmax = mmin - ybuffer, mmax + ybuffer

#  set axis limits
ax1.set_xlim((xminmax[0], xminmax[1]))
ax1.set_ylim((minmax[0], minmax[1]))

# legend
ax1.legend(ncol=2, fontsize=8)

# Make axes use scientific notation
# Only x-axis if on log scale
if ylog:
    ax1.set_yscale('log')
    ax1.ticklabel_format(axis='x', scilimits = (-3,4), useMathText=True)
else:
    ax1.ticklabel_format(scilimits = (-3,4), useMathText=True)

# Make the second plot (kinetic energy of the mean motions)
ax2.plot(xaxis, mke, 'k', linewidth=lw)
ax2.plot(xaxis, mrke, 'r', linewidth=lw)
ax2.plot(xaxis, mtke, 'g', linewidth=lw)
ax2.plot(xaxis, mpke, 'b', linewidth=lw)

# If magnetic, plot magnetic energies!
if magnetism:
    ax2.plot(xaxis, mme, 'k--', linewidth=lw)
    ax2.plot(xaxis, mrme, 'r--', linewidth=lw)
    ax2.plot(xaxis, mtme, 'g--', linewidth=lw)
    ax2.plot(xaxis, mpme, 'b--', linewidth=lw)

# Title and axis label
ax2.set_title('mean energy')

# Put the y-label on the middle plot
ax2.set_ylabel(r'$\rm{energy\ density\ (erg}\ cm^{-3})$')

# Third plot: fluctuating energy
ax3.plot(xaxis, fke, 'k', linewidth=lw)
ax3.plot(xaxis, frke, 'r', linewidth=lw)
ax3.plot(xaxis, ftke, 'g', linewidth=lw)
ax3.plot(xaxis, fpke, 'b', linewidth=lw)

# If magnetic, plot magnetic energies!
if magnetism:
    ax3.plot(xaxis, fme, 'k--', linewidth=lw)
    ax3.plot(xaxis, frme, 'r--', linewidth=lw)
    ax3.plot(xaxis, ftme, 'g--', linewidth=lw)
    ax3.plot(xaxis, fpme, 'b--', linewidth=lw)

# title and x-axis label
ax3.set_title('fluctuating energy')

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
