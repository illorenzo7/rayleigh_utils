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
the_file = None

# Set defaults
xiter = False
from0 = False
magnetism = False
ylog = False

minmax = None
fminmax = None
mminmax = None
ymin = None
ymax = None

minmax_cz = None
fminmax_cz = None
mminmax_cz = None
ymin_cz = None
ymax_cz = None

minmax_rz = None
fminmax_rz = None
mminmax_rz = None
ymin_rz = None
ymax_rz = None

xminmax = None
xmin = None
xmax = None

plot_inte = False
inte_subt = False # subtracts top value of S for inte
inte_subb = False # subtracts bot value of S for inte
sep_czrz = False # plots two more columns with energies in CZ and RZ 
    # separately 
inte_gtr2 = False # plots just the regular inte but using the 
    # trace_2dom_G_Avgs file 
    # (just to compare and make sure nothing's wonky)
plot_tote = False
savename = None
tol = 0.05
mtol = 0.9 # for exponential growth need a HIGH tolerance
plot_equil_time = False
plot_mag_equil_time = False
chosen_eqtime = None # choose eq time a priori
leg_loc = 'best'
tag = ''

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
    elif arg == '-fminmax':
        fminmax = float(args[i+1]), float(args[i+2])
    elif arg == '-mminmax':
        mminmax = float(args[i+1]), float(args[i+2])
    elif arg == '-min':
        ymin = float(args[i+1])
    elif arg == '-max':
        ymax = float(args[i+1])

    elif arg == '-minmaxcz':
        minmax_cz = float(args[i+1]), float(args[i+2])
    elif arg == '-fminmaxcz':
        fminmax_cz = float(args[i+1]), float(args[i+2])
    elif arg == '-mminmaxcz':
        mminmax_cz = float(args[i+1]), float(args[i+2])
    elif arg == '-mincz':
        ymin_cz = float(args[i+1])
    elif arg == '-maxcz':
        ymax_cz = float(args[i+1])

    elif arg == '-minmaxrz':
        minmax_rz = float(args[i+1]), float(args[i+2])
    elif arg == '-fminmaxrz':
        fminmax_rz = float(args[i+1]), float(args[i+2])
    elif arg == '-mminmaxrz':
        mminmax_rz = float(args[i+1]), float(args[i+2])
    elif arg == '-minrz':
        ymin_rz = float(args[i+1])
    elif arg == '-maxrz':
        ymax_rz = float(args[i+1])

    elif arg == '-xminmax':
        xminmax = float(args[i+1]), float(args[i+2])
    elif arg == '-xmin':
        xmin = float(args[i+1])
    elif arg == '-xmax':
        xmax = float(args[i+1])
    elif arg == '-inte':
        plot_inte = True
    elif arg == '-gtr2':
        plot_inte = True
        inte_gtr2 = True
    elif arg == '-subt':
        plot_inte = True
        inte_subt = True
    elif arg == '-subb':
        plot_inte = True
        inte_subb = True
    elif arg == '-tote':
        plot_tote = True
    elif arg == '-czrz':
        sep_czrz = True
    elif arg == '-name':
        savename = args[i+1] + '.png'
    elif arg == '-tol':
        tol = float(args[i+1])
    elif arg == '-mtol':
        mtol = float(args[i+1])
    elif arg == '-equil':
        plot_equil_time = True
    elif arg == '-mequil':
        plot_mag_equil_time = True
    elif arg == '-eqtime':
        chosen_eqtime = float(args[i+1])
    elif arg == '-legloc':
        leg_loc = args[i+1]
    elif arg == '-tag':
        tag = '_' + args[i+1]

# Tag the plot by whether or not the x axis is in "time" or "iteration"
if xiter and tag == '':
    tag = '_xiter'

# Might need to use trace_2dom_G_Avgs file instead of trace_G_Avgs
if inte_gtr2 or inte_subt or inte_subb or sep_czrz:
    if the_file is None:
        the_file = get_widest_range_file(datadir, 'trace_2dom_G_Avgs')
    print ('Using 2dom trace from ' + datadir + the_file)
    di = get_dict(datadir + the_file) 
    vals = di['vals']
    vals_cz = di['vals_cz']
    vals_rz = di['vals_rz']
    lut = di['lut']
    times = di['times']
    iters = di['iters']
    iter1 = di['iter1']
    iter2 = di['iter2']

else: # otherwise get it from trace_G_Avgs
    if the_file is None:
        the_file = get_widest_range_file(datadir, 'trace_G_Avgs')
    print ('Using G_Avgs trace from ' + datadir + the_file)
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
ixmax = np.argmin(np.abs(xaxis - xminmax[1]))

# Now shorten all the "x" arrays
xaxis = xaxis[ixmin:ixmax+1]
times = times[ixmin:ixmax+1]
iters = iters[ixmin:ixmax+1]
t1 = times[0]
t2 = times[-1]

vals = vals[:, ixmin:ixmax+1]
if sep_czrz:
    vals_cz = vals_cz[:, ixmin:ixmax+1]
    vals_cz = vals_cz[:, ixmin:ixmax+1]

# Get energy densities (averaged over whole shell)
rke = vals[:, lut[402]]
tke = vals[:, lut[403]]
pke = vals[:, lut[404]]
ke = rke + tke + pke

frke = vals[:, lut[410]]
ftke = vals[:, lut[411]]
fpke = vals[:, lut[412]]
fke = frke + ftke + fpke

mrke = rke - frke
mtke = tke - ftke
mpke = pke - fpke
mke = mrke + mtke + mpke

# Get global min/max vals
#mmax = np.max(ke)
#mmin = np.min((np.min(mrke), np.min(mtke), np.min(mpke), np.min(frke),\
#        np.min(ftke), np.min(fpke)))   

# Get the magnetic energies if they are available
if magnetism:
    rme = vals[:, lut[1102]]
    tme = vals[:, lut[1103]]
    pme = vals[:, lut[1104]]
    me = rme + tme + pme

    frme = vals[:, lut[1110]]
    ftme = vals[:, lut[1111]]
    fpme = vals[:, lut[1112]]
    fme = frme + ftme + fpme

    mrme = rme - frme
    mtme = tme - ftme
    mpme = pme - fpme
    mme = mrme + mtme + mpme

    # Update global min/max vals
#    mmax = max(mmax, np.max(me))
#    mmin = min(mmin, np.min(mrme), np.min(mtme), np.min(mpme),\
#            np.min(frke), np.min(ftme), np.min(fpme))

# Check what internal energies we might need (sometimes CZ/RZ separated ones
# too just in case) these will always be available for some options)
if inte_subt:
    inte = vals[:, lut[4001]]
    inte_cz = vals_cz[:, lut[4001]]
    inte_rz = vals_rz[lut[4001]]
    print("Got SUBT internal energy trace from trace_2dom_G_Avgs")
    inte_label = "INTE SUBT"
elif inte_subb:
    inte = vals[:, lut[4002]]
    inte_cz = vals_cz[:, lut[4002]]
    inte_rz = vals_rz[lut[4002]]
    print("Got SUBB internal energy trace from trace_2dom_G_Avgs")
    inte_label = "INTE SUBB"
elif inte_gtr2: # inte not from trace_G_Avgs
    inte = vals[:, lut[4000]]
    inte_cz = vals_cz[:, lut[4000]]
    inte_rz = vals_rz[:, lut[4000]]
    print("Got internal energy trace from trace_2dom_G_Avgs")
    inte_label = "INTE W/ DRIFT"
elif plot_inte or plot_tote: 
    # try to get inte from trace_G_Avgs
    try:
        inte = vals[:, lut[701]][ixmin:ixmax + 1]
        print("Got internal energy trace from trace_G_Avgs")
        inte_label = "INTE W/ DRIFT"
    except:
        # if not available in trace_G_Avgs, set inte to zero 
        print ("Internal energy trace not available; setting to 0")
        inte = np.zeros_like(times)
        inte_label = "INTE NOT FOUND"

#if plot_inte:
    # Update global min/max vals
#    mmax = max(mmax, np.max(inte))
#    mmin = min(mmin, np.min(inte))

if plot_tote:
    tote = ke + inte
    ftote = np.copy(fke)
    mtote = mke + inte
    if magnetism:
        tote += me
        mtote += mme
        ftote += fme
    # Update global min/max vals
#    mmax = max(mmax, np.max(tote))

if sep_czrz:
    # Get energy densities (CZ and RZ separately)
    rke_cz = vals_cz[lut_gtr2[402]][ixmin:ixmax + 1]
    tke_cz = vals_cz[lut_gtr2[403]][ixmin:ixmax + 1]
    pke_cz = vals_cz[lut_gtr2[404]][ixmin:ixmax + 1]
    ke_cz = rke_cz + tke_cz + pke_cz

    frke_cz = vals_cz[lut_gtr2[410]][ixmin:ixmax + 1]
    ftke_cz = vals_cz[lut_gtr2[411]][ixmin:ixmax + 1]
    fpke_cz = vals_cz[lut_gtr2[412]][ixmin:ixmax + 1]
    fke_cz = frke_cz + ftke_cz + fpke_cz

    mrke_cz = rke_cz - frke_cz
    mtke_cz = tke_cz - ftke_cz
    mpke_cz = pke_cz - fpke_cz
    mke_cz = mrke_cz + mtke_cz + mpke_cz

    # Get global min/max vals
#    mmax_cz = np.max(ke_cz)
#    mmin_cz = np.min((np.min(mrke_cz), np.min(mtke_cz), np.min(mpke_cz),\
#            np.min(frke_cz), np.min(ftke_cz), np.min(fpke_cz)))

    rke_rz = vals_rz[lut_gtr2[402]][ixmin:ixmax + 1]
    tke_rz = vals_rz[lut_gtr2[403]][ixmin:ixmax + 1]
    pke_rz = vals_rz[lut_gtr2[404]][ixmin:ixmax + 1]
    ke_rz = rke_rz + tke_rz + pke_rz

    frke_rz = vals_rz[lut_gtr2[410]][ixmin:ixmax + 1]
    ftke_rz = vals_rz[lut_gtr2[411]][ixmin:ixmax + 1]
    fpke_rz = vals_rz[lut_gtr2[412]][ixmin:ixmax + 1]
    fke_rz = frke_rz + ftke_rz + fpke_rz

    mrke_rz = rke_rz - frke_rz
    mtke_rz = tke_rz - ftke_rz
    mpke_rz = pke_rz - fpke_rz
    mke_rz = mrke_rz + mtke_rz + mpke_rz

    # Get global min/max vals
#    mmax_rz = np.max(ke_rz)
#    mmin_rz = np.min((np.min(mrke_rz), np.min(mtke_rz), np.min(mpke_rz),\
#            np.min(frke_rz), np.min(ftke_rz), np.min(fpke_rz)))

    # Get the magnetic energies if they are available
    if magnetism:
        rme_cz = vals_cz[lut_gtr2[1102]][ixmin:ixmax + 1]
        tme_cz = vals_cz[lut_gtr2[1103]][ixmin:ixmax + 1]
        pme_cz = vals_cz[lut_gtr2[1104]][ixmin:ixmax + 1]
        me_cz = rme_cz + tme_cz + pme_cz

        frme_cz = vals_cz[lut_gtr2[1110]][ixmin:ixmax + 1]
        ftme_cz = vals_cz[lut_gtr2[1111]][ixmin:ixmax + 1]
        fpme_cz = vals_cz[lut_gtr2[1112]][ixmin:ixmax + 1]
        fme_cz = frme_cz + ftme_cz + fpme_cz

        mrme_cz = rme_cz - frme_cz
        mtme_cz = tme_cz - ftme_cz
        mpme_cz = pme_cz - fpme_cz
        mme_cz = mrme_cz + mtme_cz + mpme_cz

        # Update global min/max vals
#        mmax_cz = max(mmax_cz, np.max(me_cz))
#        mmin_cz = min(mmin_cz, np.min(mrme_cz), np.min(mtme_cz),\
#            np.min(mpme_cz), np.min(frke_cz), np.min(ftme), np.min(fpme))

        rme_rz = vals_rz[lut_gtr2[1102]][ixmin:ixmax + 1]
        tme_rz = vals_rz[lut_gtr2[1103]][ixmin:ixmax + 1]
        pme_rz = vals_rz[lut_gtr2[1104]][ixmin:ixmax + 1]
        me_rz = rme_rz + tme_rz + pme_rz

        frme_rz = vals_rz[lut_gtr2[1110]][ixmin:ixmax + 1]
        ftme_rz = vals_rz[lut_gtr2[1111]][ixmin:ixmax + 1]
        fpme_rz = vals_rz[lut_gtr2[1112]][ixmin:ixmax + 1]
        fme_rz = frme_rz + ftme_rz + fpme_rz

        mrme_rz = rme_rz - frme_rz
        mtme_rz = tme_rz - ftme_rz
        mpme_rz = pme_rz - fpme_rz
        mme_rz = mrme_rz + mtme_rz + mpme_rz

        # Update global min/max vals
#        mmax_rz = max(mmax_rz, np.max(me_rz))
#        mmin_rz = min(mmin_rz, np.min(mrme_rz), np.min(mtme_rz),\
#            np.min(mpme_rz), np.min(frke_rz), np.min(ftme), np.min(fpme))

    # Separated internal energies already taken care of
    # But should update min/max vals
    if plot_inte:
        mmax_rz = max(mmax_rz, np.max(inte_rz))
        mmax_cz = max(mmax_cz, np.max(inte_cz))

    # Get total energies if desired
    if plot_tote:
        tote_cz = ke_cz + inte_cz
        tote_rz = ke_rz + inte_rz
        ftote_cz = np.copy(fke_cz)
        ftote_rz = np.copy(fke_rz)
        mtote_cz = mke_cz + inte_cz
        mtote_rz = mke_rz + inte_rz

        if magnetism:
            tote_cz += me_cz
        # Update global min/max vals
        mmax_cz = max(mmax_cz, np.max(tote_cz))
        mmax_rz = max(mmax_rz, np.max(tote_rz))

# create figure with 3-panel columns (total, mean and fluctuating energy)
# 1 column if only full energies desired
# 3 columns if CZ/RZ separation desired
if sep_czrz:
    ncol = 3
else:
    ncol = 1
fig, axs = plt.subplots(3, ncol, figsize=(5.*ncol, 10),\
        sharex=True)
if ncol == 1: # need the axis array to consistently be doubly indexed
    axs = np.expand_dims(axs, 1)

# Make thin lines to see structure of variation
lw = 0.5

# first plot: total kinetic energy trace      
axs[0,0].plot(xaxis, ke, 'm', linewidth=lw, label=r'$\rm{KE_{tot}}$')
axs[0,0].plot(xaxis, rke, 'r', linewidth=lw, label=r'$\rm{KE}_r$')
axs[0,0].plot(xaxis, tke, 'g', linewidth=lw, label=r'$\rm{KE}_\theta$')
axs[0,0].plot(xaxis, pke, 'b', linewidth=lw, label=r'$\rm{KE}_\phi$')

# If magnetic, plot magnetic energies!
if magnetism:
    axs[0,0].plot(xaxis, me, 'm--', linewidth=lw, label=r'$\rm{ME_{tot}}$')
    axs[0,0].plot(xaxis, rme, 'r--', linewidth=lw, label=r'$\rm{ME}_r$')
    axs[0,0].plot(xaxis, tme, 'g--', linewidth=lw,\
            label=r'$\rm{ME}_\theta$')
    axs[0,0].plot(xaxis, pme, 'b--', linewidth=lw, label=r'$\rm{ME}_\phi$')

# See if various internal/total energies should be plotted
if plot_inte:
    axs[0,0].plot(xaxis, inte, 'c', linewidth=lw, label=inte_label)
if plot_tote:
    axs[0,0].plot(xaxis, tote, 'k', linewidth=lw, label='TOT E')

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
if plot_tote:
    # Compute change in energy over time, to add to label
    # Average over the last hundred rotations 
    # (1 diffusion time, if nonrotating)
    # or all time, whichever is shorter
    if rotation:
        num = 100.
    else:
        num = 1.
    it_first = np.argmin(np.abs(times/time_unit - (t2/time_unit - num)))

    gi = GridInfo(dirname + '/grid_info')
    ri, ro = np.min(gi.radius), np.max(gi.radius)
    shell_volume = 4/3*np.pi*(ro**3 - ri**3)
    dE = (tote[-1] - tote[it_first])*shell_volume
    dt = times[-1] - times[it_first]
    title += (('\nTOT E: ' + r'$\rm{\Delta E/\Delta t = %1.3e\ cgs}$')\
            %(dE/dt))

axs[0,0].set_title(title)

#if minmax is None:
#    if ylog:
#        minmax = mmin/3.0, mmax*3.0
#    else:
#        ydiff = mmax - mmin
#        ybuffer = 0.05*ydiff 
#        minmax = mmin - ybuffer, mmax + ybuffer

# See if y-axis should be on log scale
# Make axes use scientific notation
if ylog:
    for ax in axs[:, 0]:
        ax.set_yscale('log')
        ax.ticklabel_format(axis='x', scilimits = (-3,4), useMathText=True)
else:
    for ax in axs[:, 0]:
        ax.ticklabel_format(scilimits = (-3,4), useMathText=True)

#  set axis limits
axs[0,0].set_xlim((xminmax[0], xminmax[1]))
if minmax is None:
    minmax = axs[0, 0].get_ylim()
# Change JUST ymin or ymax, if desired
if not ymin is None:
    minmax = ymin, minmax[1]
if not ymax is None:
    minmax = minmax[0], ymax
axs[0,0].set_ylim((minmax[0], minmax[1]))

# Calculate equilibration time (for KE) if desired
if plot_equil_time:
    print ("=============================")
    print ("calculating equil time for KE")
    print ("tol = %0.03f" %tol)
    print ("to use a different value, type -tol [val]")
    max_ke = np.max(ke)
    foundit = False
    i_equil = 0
    if chosen_eqtime is None:
        for i in range(len(ke)):
            if ke[i] > (1.0 - tol)*max_ke and not foundit:
                foundit = True
                i_equil = np.copy(i)
    else:
        i_equil = np.argmin(np.abs(xaxis - chosen_eqtime))
    x_equil = xaxis[i_equil]
    ke_equil = np.mean(ke[i_equil:])
    xvals = np.linspace(xminmax[0], xminmax[1], 100)
    yvals = np.linspace(minmax[0], minmax[1], 100)
    axs[0,0].plot(x_equil + np.zeros(100), yvals, 'k--')
    axs[0,0].plot(xvals, ke_equil + np.zeros(100), 'k--')
    axs[0,0].text(x_equil + 0.05*(x_equil - xminmax[0]),\
            0.25*ke_equil, ('t = %1.2e\ntol = %.03f' %(x_equil, tol)),\
            ha='left', va='center')
    axs[0,0].text(0.05*(x_equil - xminmax[0]),\
            0.95*ke_equil, ('KE = %1.2e' %ke_equil), ha='left', va='top')

if plot_mag_equil_time:
    print ("=============================")
    print ("calculating equil time for mag. E")
    print ("mtol = %0.03f" %mtol)
    print ("to use a different value, type -tol [val]")
    max_me = np.max(me)
    foundit = False
    i_equil = 0
    for i in range(len(me)):
        if me[i] > (1.0 - mtol)*max_me and not foundit:
            foundit = True
            i_equil = np.copy(i)
    x_equil = xaxis[i_equil]
    Dx_equil = x_equil - xminmax[0]
    me_equil = np.mean(me[i_equil:])
    xvals = np.linspace(xminmax[0], xminmax[1], 100)
    yvals = np.linspace(minmax[0], minmax[1], 100)
    axs[0,0].plot(x_equil + np.zeros(100), yvals, 'k--')
    axs[0,0].plot(xvals, me_equil + np.zeros(100), 'k--')
    if ylog:
        dynrange = me_equil/minmax[0]
        ycoord = minmax[0]*dynrange**0.25
    else:
        ycoord = 0.25*me_equil
    axs[0,0].text(x_equil + 0.05*(x_equil - xminmax[0]),\
            ycoord, ('t = %1.2e\n' + r'$\Delta t$' +\
            ' = %1.2e\nmtol = %.03f')\
            %(x_equil, Dx_equil, mtol), ha='left', va='center')
    axs[0,0].text(0.05*(x_equil - xminmax[0]),\
            0.95*me_equil, ('KE = %1.2e' %me_equil), ha='left', va='top')

# legend
axs[0,0].legend(ncol=2, fontsize=8, loc=leg_loc)

# Make the second plot (energy of the mean motions)
axs[1,0].plot(xaxis, mke, 'm', linewidth=lw)
axs[1,0].plot(xaxis, mrke, 'r', linewidth=lw)
axs[1,0].plot(xaxis, mtke, 'g', linewidth=lw)
axs[1,0].plot(xaxis, mpke, 'b', linewidth=lw)

# If magnetic, plot magnetic energies!
if magnetism:
    axs[1,0].plot(xaxis, mme, 'm--', linewidth=lw)
    axs[1,0].plot(xaxis, mrme, 'r--', linewidth=lw)
    axs[1,0].plot(xaxis, mtme, 'g--', linewidth=lw)
    axs[1,0].plot(xaxis, mpme, 'b--', linewidth=lw)

# See if various internal/total energies should be plotted
if plot_inte:
    axs[1,0].plot(xaxis, inte, 'c', linewidth=lw, label=inte_label)
if plot_tote:
    axs[1,0].plot(xaxis, mtote, 'k', linewidth=lw, label='TOT E')

#  set axis limits
if mminmax is None:
    mminmax = axs[1,0].get_ylim()
# Change JUST ymin or ymax, if desired
if not ymin is None:
    mminmax = ymin, mminmax[1]
if not ymax is None:
    mminmax = mminmax[0], ymax
axs[1,0].set_ylim((mminmax[0], mminmax[1]))

# Title and axis label
axs[1,0].set_title('mean energy')

# Put the y-label on the middle plot
axs[1,0].set_ylabel(r'$\rm{energy\ density\ (erg}\ cm^{-3})$')

# Third plot: fluctuating energy
axs[2,0].plot(xaxis, fke, 'k', linewidth=lw)
axs[2,0].plot(xaxis, frke, 'r', linewidth=lw)
axs[2,0].plot(xaxis, ftke, 'g', linewidth=lw)
axs[2,0].plot(xaxis, fpke, 'b', linewidth=lw)

# If magnetic, plot magnetic energies!
if magnetism:
    axs[2,0].plot(xaxis, fme, 'k--', linewidth=lw)
    axs[2,0].plot(xaxis, frme, 'r--', linewidth=lw)
    axs[2,0].plot(xaxis, ftme, 'g--', linewidth=lw)
    axs[2,0].plot(xaxis, fpme, 'b--', linewidth=lw)

# See if various internal/total energies should be plotted
if plot_tote:
    axs[2,0].plot(xaxis, ftote, 'k', linewidth=lw, label='TOT E')

# title and x-axis label
axs[2,0].set_title('fluctuating energy')

#  set axis limits
if fminmax is None:
    fminmax = axs[2,0].get_ylim()
# Change JUST ymin or ymax, if desired
if not ymin is None:
    fminmax = ymin, fminmax[1]
if not ymax is None:
    fminmax = fminmax[0], ymax
axs[2,0].set_ylim((fminmax[0], fminmax[1]))

# Put the x-axis label on the bottom
if (xiter):
    axs[2,0].set_xlabel('iteration #')
else:
    if rotation:
        axs[2,0].set_xlabel(r'$\rm{t\ (P_{rot})}$')
    else:
        axs[2,0].set_xlabel(r'$\rm{t\ (T_{diff})}$')

# Get ticks everywhere
for ax in axs.flatten():
    plt.sca(ax)
    plt.minorticks_on()
    plt.tick_params(top=True, right=True, direction='in', which='both')

# Space the subplots to make them look pretty
plt.tight_layout
plt.subplots_adjust(left=0.15, bottom=0.08, top=0.85, wspace=0.4)

if savename is None:
    savename = dirname_stripped + '_etrace_' + str(iter1).zfill(8) + '_' + str(iter2).zfill(8) + tag + '.png'

# Save the plot
print ('Saving the etrace plot at ' + plotdir + savename)
plt.savefig(plotdir + savename, dpi=300)

# Show the plot
plt.show()
