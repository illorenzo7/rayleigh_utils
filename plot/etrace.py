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
from reference_tools import equation_coefficients
from rayleigh_diagnostics import ReferenceState, GridInfo

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
trace_G_Avgs_file = get_widest_range_file(datadir, 'trace_G_Avgs')

# Set defaults
xiter = False
notfrom0 = False
magnetism = False
ylog = False
minmax = None
xminmax = None
plot_inte = False
plot_tote = False

# Get command-line arguments
args = sys.argv[2:]
nargs = len(args)
for i in range(nargs):
    arg = args[i]
    if arg == '-xiter': # plot w.r.t. iterations
        xiter = True
    elif arg == '-usefile':
        trace_G_Avgs_file = args[i+1]
        trace_G_Avgs_file = trace_G_Avgs_file.split('/')[-1]
    elif arg == '-notfrom0':
        notfrom0 = True
    elif arg == '-mag':
        magnetism = True
    elif arg == '-log':
        ylog = True
    elif arg == '-minmax':
        minmax = float(args[i+1]), float(args[i+2])
    elif arg == '-xminmax':
        xminmax = float(args[i+1]), float(args[i+2])
    elif arg == '-inte':
        plot_inte = True
    elif arg == '-tote':
        plot_tote = True

# Tag the plot by whether or not the x axis is in "time" or "iteration"
if (xiter):
    tag = '_xiter'
else:
    tag = '_xtime'

# Read in the KE data (dictionary form)
print ('Getting energy trace from ' + datadir + trace_G_Avgs_file + ' ...')
di = get_dict(datadir + trace_G_Avgs_file)

vals = di['vals']
lut = di['lut']
times = di['times']
iters = di['iters']
iter1 = di['iter1']
iter2 = di['iter2']

# Get global rotation rate, if present
rotation = get_parameter(dirname, 'rotation')
if rotation:
    angular_velocity = get_parameter(dirname, 'angular_velocity')
    Prot = 2*np.pi/angular_velocity
    tnorm = Prot # normalize the time axis by rotation period if applicable
else:
    try:
        trans = TransportCoeffs(dirname + '/reference', '')
        ktop = trans.kappa[0]
    except:
        eq = equation_coefficients()
        eq.read(dirname + '/equation_coefficients')
        ktop = eq.functions[4][0]
#    ktop = get_parameter(dirname, 'kappa_top')
    try:
        rmin = get_parameter(dirname, 'rmin')
        rmax = get_parameter(dirname, 'rmax')
    except: # two domains stitched together
        domain_bounds = get_parameter(dirname, 'domain_bounds')
        rmin = np.min(domain_bounds)
        rmax = np.max(domain_bounds)
    depth = rmax - rmin
    tdt = depth**2/ktop
    tnorm = tdt

# Make appropriate file name to save
savename = dirname_stripped + '_etrace_' + str(iter1).zfill(8) + '_' + str(iter2).zfill(8) + tag + '.png'

ke = vals[lut[401]]
rke = vals[lut[402]]
tke = vals[lut[403]]
pke = vals[lut[404]]

mke = vals[lut[405]]
mrke = vals[lut[406]]
mtke = vals[lut[407]]
mpke = vals[lut[408]]

fke = vals[lut[409]]
frke = vals[lut[410]]
ftke = vals[lut[411]]
fpke = vals[lut[412]]

# Get the magnetic energies if they are available
if magnetism:
    me = vals[lut[1101]]
    rme = vals[lut[1102]]
    tme = vals[lut[1103]]
    pme = vals[lut[1104]]

    mme = vals[lut[1105]]
    mrme = vals[lut[1106]]
    mtme = vals[lut[1107]]
    mpme = vals[lut[1108]]

    fme = vals[lut[1109]]
    frme = vals[lut[1110]]
    ftme = vals[lut[1111]]
    fpme = vals[lut[1112]]

if plot_inte or plot_tote:
    try: # get the internal energy if it is available
        try: # First try to get it from Shell_Avgs data
            the_file = get_widest_range_file(datadir,\
                    'inte_from_Shell_Avgs')
            di_inte = get_dict(datadir + the_file)
            inte = di_inte['inte']
            print("Got internal energy from Shell_Avgs")
        except:
            inte = vals[lut[701]]
            print("Got internal energy from G_Avgs")
    except:
        print ("Internal energy not available; setting to 0")
        int_e = np.zeros_like(times)

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
    mmax = max(np.max(np.abs(inte)), mmax)
    mmin = min(np.min(np.abs(inte)), mmin)

if plot_tote:
    mmax = max(np.max(np.abs(tote)), mmax)
    mmin = min(np.min(np.abs(tote)), mmin)

if not xiter:
    xaxis = times/tnorm
else:
    xaxis = iters

if notfrom0:
    x_min = np.min(xaxis)
else:
    x_min = 0

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
    sign = np.sign(np.mean(inte))
    if sign == -1:
        sign_str1 = 'negative'
    else:
        sign_str1 = 'positive'
    ax1.plot(xaxis, np.abs(inte), 'm', linewidth=lw, label='|INT E|')

if plot_tote:
    sign = np.sign(np.mean(tote))
    if sign == -1:
        sign_str2 = 'negative'
    else:
        sign_str2 = 'positive'
    ax1.plot(xaxis, np.abs(tote), 'c', linewidth=lw, label='|TOT E|')

title = dirname_stripped + '\n ' +\
          str(iter1).zfill(8) + ' to ' + str(iter2).zfill(8) +\
          '\ntotal energy'
if plot_inte:
    # Compute change in energy over time, to add to label
    gi = GridInfo(dirname + '/grid_info')
    ri, ro = np.min(gi.radius), np.max(gi.radius)
    shell_volume = 4/3*np.pi*(ro**3 - ri**3)
    dE = (inte[-1] - inte[0])*shell_volume
    dt = times[-1] - times[0]
    title += (('\n' + r'$\rm{\ sign(\overline{INT\ E})\ = %s}$' + '\n' +\
            r'$\rm{\Delta E/\Delta t = %1.3e\ cgs}$') %(sign_str1, dE/dt))
if plot_tote:
    gi = GridInfo(dirname + '/grid_info')
    ri, ro = np.min(gi.radius), np.max(gi.radius)
    shell_volume = 4/3*np.pi*(ro**3 - ri**3)
    dE = (tote[-1] - tote[0])*shell_volume
    dt = times[-1] - times[0]
    title += (('\n' + r'$\rm{\ sign(\overline{TOT\ E})\ = %s}$' + '\n' +\
            r'$\rm{\Delta E/\Delta t = %1.3e\ cgs}$') %(sign_str2, dE/dt))

ax1.set_title(title)

if minmax is None:
    if ylog:
        minmax = mmin/3.0, mmax*3.0
    else:
        ydiff = mmax - mmin
        ybuffer = 0.05*ydiff 
        minmax = mmin - ybuffer, mmax + ybuffer

ax1.set_ylim((minmax[0], minmax[1]))
    
# Set x limits  
if xminmax is None:
    xminmax = x_min, np.max(xaxis)

ax1.set_xlim((xminmax[0], xminmax[1]))

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
