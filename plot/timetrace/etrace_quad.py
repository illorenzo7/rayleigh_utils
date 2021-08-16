# Author: Loren Matilsky
# plots energies in different quadrants (hence, "quad") of the meridional plane
import matplotlib.pyplot as plt
import numpy as np
import sys, os
sys.path.append(os.environ['raco'])
sys.path.append(os.environ['rapp'])
sys.path.append(os.environ['rapl'])
from common import *
from plotcommon import *
from cla_util import *
from rayleigh_diagnostics import GridInfo

# Get the run directory on which to perform the analysis
args = sys.argv
clas0, clas = read_clas(args)
dirname = clas0['dirname']
dirname_stripped = strip_dirname(dirname)

# See if magnetism is "on"
magnetism = get_parameter(dirname, 'magnetism')

# SPECIFIC ARGS for etrace_quad:
kwargs_default = dict({'the_file': None, 'xminmax': None, 'xmin': None, 'xmax': None, 'minmax': None, 'min': None, 'max': None, 'coords': None, 'ntot': 500, 'xiter': False, 'log': False, 'nodyn': False, 'dynfrac': 0.5, 'xvals': np.array([]), 'inte': False})

# update these defaults from command-line
kwargs = update_dict(kwargs_default, clas)

fontsize = default_titlesize
the_file = kwargs.the_file
xminmax = kwargs.xminmax
xmin = kwargs.xmin
xmax = kwargs.xmax
minmax = kwargs.minmax
ymin = kwargs.min
ymax = kwargs.max
coords = kwargs.coords
ntot = kwargs.ntot
xiter = kwargs.xiter
logscale = kwargs.log
nodyn = kwargs.nodyn
dynfrac = kwargs.dynfrac
xvals = make_array(kwargs.xvals)
plot_inte = kwargs.inte

# deal with coords (if user wants minmax to only apply to certain subplots)
if not coords is None:
    numpanels = len(coords)//2
    acopy = np.copy(coords)
    coords = []
    for i in range(numpanels):
        coords.append((acopy[2*i], acopy[2*i + 1]))

# Find the etrace file(s) in the data directory. If there are multiple, by
# default choose the one with widest range in the trace.
the_file = None

# Set defaults
ntot = 500 # default number of x-axis points to use in plt.plot
xiter = False
from0 = False
magnetism = None
ylog = False
nodyn = False # by default don't adjust the min val to ignore super small 
    # magnetic energies during dynamo growth when ylog=True (i.e., plot
    # the dynamo growth phase by default)
    # to change, use -nodyn/ -dynfrac [val] to ignore the ME values over
    # the last dynfrac of the simulation
dyn_frac = 1./2.
mes = None
subinte = True # by default shift the internal energies by a constant
    # so they aren't so huge
leak_frac = 1./4. # compute leak averaged over last quarter of simulation
tag = ''

xminmax = None
xmin = None
xmax = None

plot_inte = False
inte_subt = False # subtracts top value of S for inte
inte_subb = False # subtracts bot value of S for inte
plot_tote = False

# Get command-line arguments
plotdir = None

args = sys.argv[2:]
nargs = len(args)
for i in range(nargs):
    arg = args[i]
    if arg == '-plotdir':
        plotdir = args[i+1]
    if arg == '-ntot': # plot w.r.t. iterations
        ntot = int(args[i+1])
    elif arg == '-xiter': # plot w.r.t. iterations
        xiter = True
    elif arg == '-usefile':
        the_file = args[i+1]
        the_file = the_file.split('/')[-1]
    elif arg == '-from0':
        from0 = True
    elif arg == '-mag':
        magnetism = bool(args[i+1])
    elif arg == '-log':
        ylog = True
    elif arg == '-nodyn':
        nodyn = True
    elif arg == '-dynfrac':
        dyn_frac = float(args[i+1])
    elif arg == '-fullinte':
        subinte = False
    elif arg == '-frac':
        leak_frac = float(args[i+1])

    elif arg == '-xminmax':
        xminmax = float(args[i+1]), float(args[i+2])
    elif arg == '-xmin':
        xmin = float(args[i+1])
    elif arg == '-xmax':
        xmax = float(args[i+1])
    elif arg == '-inte':
        plot_inte = True
    elif arg == '-subt':
        plot_inte = True
        inte_subt = True
    elif arg == '-subb':
        plot_inte = True
        inte_subb = True
    elif arg == '-tote':
        plot_tote = True
    elif arg == '-tag':
        tag = '_' + args[i+1]

# by default determine magnetism from main_input
if magnetism is None:
    magnetism = get_parameter(dirname, 'magnetism')

# Tag the plot by whether or not the x axis is in "time" or "iteration"
if xiter and tag == '':
    tag = '_xiter'

# read the data
if the_file is None:
    the_file = get_widest_range_file(datadir, 'trace_quad_G_Avgs')
di = get_dict(datadir + the_file)
vals = di['vals']
vals_full = di['vals_full']
volumes = di['volumes']
latbounds = di['latbounds']
rbounds = di['rbounds']
nsep_t = di['nsep_t']
nsep_r = di['nsep_r']
nquad = di['nquad']
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

if plotdir is None:
    plotdir = dirname + '/plots/'
    if not os.path.isdir(plotdir):
        os.makedirs(plotdir)

# Take slices based on what xminmax is
if not xiter:
    xaxis = times/time_unit
else:
    xaxis = iters

# only need to do stuff here of xminmax was not set by user
if xminmax is None:
    # set xmin possibly
    if xmin is None:
        if from0:
            xmin = 0.
        else:
            xmin = np.min(xaxis)
    # set xmax possibly
    if xmax is None:
        xmax = np.max(xaxis)
    xminmax = xmin, xmax

ixmin = np.argmin(np.abs(xaxis - xminmax[0]))
ixmax = np.argmin(np.abs(xaxis - xminmax[1]))

# Now shorten all the "x" arrays
xaxis = xaxis[ixmin:ixmax+1]
times = times[ixmin:ixmax+1]
iters = iters[ixmin:ixmax+1]
t1 = times[0]
t2 = times[-1]

vals = vals[ixmin:ixmax+1, :, :, :]

# Thin out the arrays to not deal obscene quantities of data 
# (and unreadable "curves")
def thin_data(vals, ntot):
    nx = np.shape(vals)[0]
    nskip = nx//ntot
    if not nskip in [0, 1]: #for ntot < 2*nx, do nothing
        vals_new = vals[::nskip]
    else:
        vals_new = vals
    return vals_new
print ("ntot = %i" %ntot)
print ("before thin_data: len(xaxis) = %i" %len(xaxis))
xaxis = thin_data(xaxis, ntot)
times = thin_data(times, ntot)
iters = thin_data(iters, ntot)
vals = thin_data(vals, ntot)
print ("after thin_data: len(xaxis) = %i" %len(xaxis))

# Get reference state in case it's needed
eq = get_eq(dirname)
rhot = eq.density*eq.temperature
# radial integration weights for integrals of rho * T
gi = GridInfo(dirname + '/grid_info', '')
rr = gi.radius
nr = gi.nr
# luminosity
lstar = get_lum(dirname)

# set up figure axes
nrow = nsep_r + 1
ncol = nsep_t + 1

# figure for total energies
fig, axs = plt.subplots(nrow, ncol, figsize=(5.*ncol, 10),\
    sharex=True)
# if it's not an array, make it one
axs = np.array(axs)
if nrow == 1:
    axs = np.expand_dims(axs, 0)
if ncol == 1: # need the axis array to consistently be doubly indexed
    axs = np.expand_dims(axs, 1)

# figure for mean energies
fig_mean, axs_mean = plt.subplots(nrow, ncol, figsize=(5.*ncol, 10),\
    sharex=True)
# if it's not an array, make it one
axs_mean = np.array(axs_mean)
if nrow == 1:
    axs_mean = np.expand_dims(axs_mean, 0)
if ncol == 1: # need the axis array to consistently be doubly indexed
    axs_mean = np.expand_dims(axs_mean, 1)

# figure for fluc energies
fig_fluc, axs_fluc = plt.subplots(nrow, ncol, figsize=(5.*ncol, 10),\
    sharex=True)
# if it's not an array, make it one
axs_fluc = np.array(axs_fluc)
if nrow == 1:
    axs_fluc = np.expand_dims(axs_fluc, 0)
if ncol == 1: # need the axis array to consistently be doubly indexed
    axs_fluc = np.expand_dims(axs_fluc, 1)

# Make thin lines to see structure of variation for ME
lw = 0.5
lw_ke = 1. # bit thicker for KE to differentiate between ME

# loop over quadrants and make plots
for ir in range(nrow):
    for it in range(ncol):
        # Get energy densities (averaged over whole shell)
        rke = vals[:, lut[402], it, ir]
        tke = vals[:, lut[403], it, ir]
        pke = vals[:, lut[404], it, ir]
        ke = rke + tke + pke

        frke = vals[:, lut[410], it, ir]
        ftke = vals[:, lut[411], it, ir]
        fpke = vals[:, lut[412], it, ir]
        fke = frke + ftke + fpke

        mrke = rke - frke
        mtke = tke - ftke
        mpke = pke - fpke
        mke = mrke + mtke + mpke

        # Get the magnetic energies if they are available
        if magnetism:
            rme = vals[:, lut[1102], it, ir]
            tme = vals[:, lut[1103], it, ir]
            pme = vals[:, lut[1104], it, ir]
            me = rme + tme + pme

            frme = vals[:, lut[1110], it, ir]
            ftme = vals[:, lut[1111], it, ir]
            fpme = vals[:, lut[1112], it, ir]
            fme = frme + ftme + fpme

            mrme = rme - frme
            mtme = tme - ftme
            mpme = pme - fpme
            mme = mrme + mtme + mpme

        # Check what internal energies we might need 
        if inte_subt:
            inte = vals[:, lut[4001], it, ir]
            print("Got SUBT internal energy trace")
            inte_label = "IE SUBT"
        elif inte_subb:
            inte = vals[:, lut[4002], it, ir]
            print("Got SUBB internal energy trace")
            inte_label = "IE SUBB"
        elif plot_inte or plot_tote:
            # inte not from trace_G_Avgs
            inte = vals[:, lut[4000], it, ir]
            print("Got FULL inte (including drift)")
            inte_label = "IE W/ DRIFT"

        # possibly shift inte values to lie just above max ke
        if plot_inte and subinte:
            min_inte = np.min(inte)
            max_ke = np.max(ke)
            diff = min_inte - max_ke
            buff = max_ke*0.05
            sub = diff - buff
            inte -= sub
         
        # Get total energy if desired
        if plot_tote:
            tote = ke + inte
            ftote = np.copy(fke)
            mtote = mke + inte
            if magnetism:
                tote += me
                mtote += mme
                ftote += fme

            # Will need to compute "energy leaks":
            it_leak = np.argmin(np.abs((times - t1) - (1. - leak_frac)*(t2 - t1)))

            # Compute leaks (full shell)
            # best fit line to last part of trace:
            times_leak = np.copy(times[it_leak:])
            tote_leak = np.copy(tote[it_leak:])
            mtote_leak = np.copy(mtote[it_leak:])
            ftote_leak = np.copy(ftote[it_leak:])
           
            m_leak, b_leak = np.polyfit(times_leak, tote_leak, 1)
            mm_leak, mb_leak = np.polyfit(times_leak, mtote_leak, 1)
            fm_leak, fb_leak = np.polyfit(times_leak, ftote_leak, 1)

            volume = volumes[it, ir]
            dEdt = m_leak*volume/lstar 
            mdEdt = mm_leak*volume/lstar
            fdEdt = fm_leak*volume/lstar

       # See if y-axis should be on log scale (must do this before setting limits)
        # Make all axes use scientific notation (except for y if ylog=True)
        for ax in np.hstack((axs.flatten(), axs_mean.flatten(),\
                axs_fluc.flatten())):
            if ylog:
                ax.set_yscale('log')
                ax.ticklabel_format(axis='x', scilimits=(-3,4), useMathText=True)
            else:
                ax.ticklabel_format(scilimits = (-3,4), useMathText=True)

        # plot total energies
        ax = axs[ir, it]
        ax.plot(xaxis, ke, 'm', linewidth=lw_ke, label=r'$\rm{KE_{tot}}$')
        ax.plot(xaxis, rke, 'r', linewidth=lw_ke, label=r'$\rm{KE_r}$')
        ax.plot(xaxis, tke, 'g', linewidth=lw_ke, label=r'$\rm{KE_\theta}$')
        ax.plot(xaxis, pke, 'b', linewidth=lw_ke, label=r'$\rm{KE_\phi}$')

        # plot mean KE
        ax_mean = axs_mean[ir, it]
        ax_mean.plot(xaxis, mke, 'm', linewidth=lw_ke)
        ax_mean.plot(xaxis, mrke, 'r', linewidth=lw_ke)
        ax_mean.plot(xaxis, mtke, 'g', linewidth=lw_ke)
        ax_mean.plot(xaxis, mpke, 'b', linewidth=lw_ke)

        # plot fluc KE
        ax_fluc = axs_fluc[ir, it]
        ax_fluc.plot(xaxis, fke, 'm', linewidth=lw_ke)
        ax_fluc.plot(xaxis, frke, 'r', linewidth=lw_ke)
        ax_fluc.plot(xaxis, ftke, 'g', linewidth=lw_ke)
        ax_fluc.plot(xaxis, fpke, 'b', linewidth=lw_ke)

        # If magnetic, plot magnetic energies!
        if magnetism:
            # plot total ME
            ax.plot(xaxis, me, 'm--', linewidth=lw, label=r'$\rm{ME_{tot}}$')
            ax.plot(xaxis, rme, 'r--', linewidth=lw, label=r'$\rm{ME_r}$')
            ax.plot(xaxis, tme, 'g--', linewidth=lw,\
                    label=r'$\rm{ME_\theta}$')
            ax.plot(xaxis, pme, 'b--', linewidth=lw, label=r'$\rm{ME_\phi}$')

            # plot mean ME
            ax_mean.plot(xaxis, mme, 'm--', linewidth=lw)
            ax_mean.plot(xaxis, mrme, 'r--', linewidth=lw)
            ax_mean.plot(xaxis, mtme, 'g--', linewidth=lw)
            ax_mean.plot(xaxis, mpme, 'b--', linewidth=lw)

            # plot fluc ME
            ax_fluc.plot(xaxis, fme, 'm--', linewidth=lw)
            ax_fluc.plot(xaxis, frme, 'r--', linewidth=lw)
            ax_fluc.plot(xaxis, ftme, 'g--', linewidth=lw)
            ax_fluc.plot(xaxis, fpme, 'b--', linewidth=lw)

        # See if various internal/total energies should be plotted
        if plot_inte:
            ax.plot(xaxis, inte, 'c', linewidth=lw_ke, label=inte_label)
            ax_mean.plot(xaxis, inte, 'c', linewidth=lw_ke)
        if plot_tote:
            ax.plot(xaxis, tote, 'k', linewidth=lw_ke, label=r'$\rm{TE}$')
            ax_mean.plot(xaxis, mtote, 'k', linewidth=lw_ke)
            ax_fluc.plot(xaxis, ftote, 'k', linewidth=lw_ke)

            # Label energy leaks with rate of leak and dashed red line
            # full
            leak_label = r'$\rm{(\overline{dE/dt})/L_* = %09.3e}$'
            leak_label_x = 0.99
            leak_label_y = 0.01
            ax.text(leak_label_x, leak_label_y, leak_label %dEdt,\
                    va='bottom', ha='right', transform=ax.transAxes)
            ax.plot(xaxis[it_leak:], m_leak*times_leak + b_leak, 'r--') 
            # mean
            ax_mean.text(leak_label_x, leak_label_y, leak_label %mdEdt,\
                    va='bottom', ha='right', transform=ax.transAxes)
            ax_mean.plot(xaxis[it_leak:], mm_leak*times_leak + mb_leak, 'r--') 
            # fluc
            ax_fluc.text(leak_label_x, leak_label_y, leak_label %fdEdt,\
                    va='bottom', ha='right', transform=ax.transAxes)
            ax_fluc.plot(xaxis[it_leak:], fm_leak*times_leak + fb_leak, 'r--') 

# put a legend on the upper left axis
axs[0,0].legend(loc='lower left', ncol=3, fontsize=8, columnspacing=1)

# Set some parameters defining all subplots
# x limits and label
axs[0,0].set_xlim((xminmax[0], xminmax[1]))
if xiter:
    axs[nrow-1,0].set_xlabel('iteration #')
else:
    if rotation:
        axs[nrow-1,0].set_xlabel(r'$\rm{t\ (P_{rot})}$')
    else:
        axs[nrow-1,0].set_xlabel(r'$\rm{t\ (T_{diff})}$')

# label size
fs = 12

# y labels
for ir in range(nrow):
    r1 = rbounds[ir]/rsun
    r2 = rbounds[ir+1]/rsun
    for ax in axs[ir,0], axs_mean[ir,0], axs_fluc[ir,0]:
        ax.set_ylabel('rbounds = [%.3f, %.3f]' %(r1, r2), fontsize=fs)

# titles
for it in range(ncol):
    lat1 = latbounds[it]
    lat2 = latbounds[it+1]
    for ax in axs[0, it], axs_mean[0, it], axs_fluc[0, it]:
        ax.set_title('latbounds = [%.1f, %.1f]' %(lat1, lat2), fontsize=fs)

# Get ticks everywhere
for ax in np.hstack((axs.flatten(), axs_mean.flatten(), axs_fluc.flatten())):
    plt.sca(ax)
    plt.minorticks_on()
    plt.tick_params(top=True, right=True, direction='in', which='both')

rootnames = ['etrace_quad', 'etrace_quad_mean', 'etrace_quad_fluc']
count = 0
for fi in fig, fig_mean, fig_fluc:
    # Space the subplots to make them look pretty
    fi.tight_layout()
    #plt.subplots_adjust(left=0.15, bottom=0.08, top=0.85, wspace=0.4)
    savename = dirname_stripped + '_' + rootnames[count] + '_' + str(iter1).zfill(8) + '_' + str(iter2).zfill(8) + tag + '.png'

    # Save the plot
    print ('Saving ' + plotdir + savename)
    fi.savefig(plotdir + savename, dpi=300)
    plt.close(fi)
    count += 1
