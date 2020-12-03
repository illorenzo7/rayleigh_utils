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
        get_dict, get_lum
from get_parameter import get_parameter
from rayleigh_diagnostics import GridInfo
from time_scales import compute_Prot, compute_tdt
from get_eq import get_eq

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
magnetism = None
ylog = False
dyn = False # by default adjust the min val to ignore super small 
    # magnetic energies during dynamo growth when ylog=True
mes = None
subinte = True # by default shift the internal energies by a constant
    # so they aren't so huge
leak_frac = 1./4. # compute leak averaged over last quarter of simulation
minmax = None
ymin = None
ymax = None
mminmax = None
mymin = None
mymax = None
fminmax = None
fymin = None
fymax = None

minmax_cz = None
ymin_cz = None
ymax_cz = None
mminmax_cz = None
mymin_cz = None
mymax_cz = None
fminmax_cz = None
fymin_cz = None
fymax_cz = None

minmax_rz = None
ymin_rz = None
ymax_rz = None
mminmax_rz = None
mymin_rz = None
mymax_rz = None
fminmax_rz = None
fymin_rz = None
fymax_rz = None

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
        magnetism = bool(args[i+1])
    elif arg == '-log':
        ylog = True
    elif arg == '-dyn':
        dyn = True
    elif arg == '-fullinte':
        subinte = False
    elif arg == '-frac':
        leak_frac = float(args[i+1])
    elif arg == '-minmax':
        minmax = float(args[i+1]), float(args[i+2])
    elif arg == '-min':
        ymin = float(args[i+1])
    elif arg == '-max':
        ymax = float(args[i+1])
    elif arg == '-mminmax':
        mminmax = float(args[i+1]), float(args[i+2])
    elif arg == '-mmin':
        mymin = float(args[i+1])
    elif arg == '-mmax':
        mymax = float(args[i+1])
    elif arg == '-fminmax':
        fminmax = float(args[i+1]), float(args[i+2])
    elif arg == '-fmin':
        fymin = float(args[i+1])
    elif arg == '-fmax':
        fymax = float(args[i+1])

    elif arg == '-minmaxcz':
        minmax_cz = float(args[i+1]), float(args[i+2])
    elif arg == '-mincz':
        ymin_cz = float(args[i+1])
    elif arg == '-maxcz':
        ymax_cz = float(args[i+1])
    elif arg == '-mminmaxcz':
        mminmax_cz = float(args[i+1]), float(args[i+2])
    elif arg == '-mmincz':
        mymin_cz = float(args[i+1])
    elif arg == '-mmaxcz':
        mymax_cz = float(args[i+1])
    elif arg == '-fminmaxcz':
        fminmax_cz = float(args[i+1]), float(args[i+2])
    elif arg == '-fmincz':
        fymin_cz = float(args[i+1])
    elif arg == '-fmaxcz':
        fymax_cz = float(args[i+1])

    elif arg == '-minmaxrz':
        minmax_rz = float(args[i+1]), float(args[i+2])
    elif arg == '-minrz':
        ymin_rz = float(args[i+1])
    elif arg == '-maxrz':
        ymax_rz = float(args[i+1])
    elif arg == '-mminmaxrz':
        mminmax_rz = float(args[i+1]), float(args[i+2])
    elif arg == '-mminrz':
        mymin_rz = float(args[i+1])
    elif arg == '-mmaxrz':
        mymax_rz = float(args[i+1])
    elif arg == '-fminmaxrz':
        fminmax_rz = float(args[i+1]), float(args[i+2])
    elif arg == '-fminrz':
        fymin_rz = float(args[i+1])
    elif arg == '-fmaxrz':
        fymax_rz = float(args[i+1])

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
    elif arg == '-tag':
        tag = '_' + args[i+1]

# by default determine magnetism from main_input
if magnetism is None:
    magnetism = get_parameter(dirname, 'magnetism')

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

vals = vals[ixmin:ixmax+1, :]
if sep_czrz:
    vals_cz = vals_cz[ixmin:ixmax+1, :]
    vals_rz = vals_rz[ixmin:ixmax+1, :]

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

# Check what internal energies we might need (sometimes CZ/RZ separated ones
# too just in case) these will always be available for some options)
if inte_subt:
    inte = vals[:, lut[4001]]
    inte_cz = vals_cz[:, lut[4001]]
    inte_rz = vals_rz[:, lut[4001]]
    print("Got SUBT internal energy trace from trace_2dom_G_Avgs")
    inte_label = "IE SUBT"
elif inte_subb:
    inte = vals[:, lut[4002]]
    inte_cz = vals_cz[:, lut[4002]]
    inte_rz = vals_rz[:, lut[4002]]
    print("Got SUBB internal energy trace from trace_2dom_G_Avgs")
    inte_label = "IE SUBB"
elif inte_gtr2 or (plot_inte and sep_czrz): # inte not from trace_G_Avgs
    inte = vals[:, lut[4000]]
    inte_cz = vals_cz[:, lut[4000]]
    inte_rz = vals_rz[:, lut[4000]]
    print("Got internal energy trace from trace_2dom_G_Avgs")
    inte_label = "IE W/ DRIFT"
elif plot_inte or plot_tote: 
    # try to get inte from trace_G_Avgs
    try:
        inte = vals[:, lut[701]]
        print("Got internal energy trace from trace_G_Avgs")
        inte_label = "IE W/ DRIFT"
    except:
        # if not available in trace_G_Avgs, set inte to zero 
        print ("Internal energy trace not available; setting to 0")
        inte = np.zeros_like(times)
        inte_label = "IE NOT FOUND"

# Get reference state in case it's needed
eq = get_eq(dirname)
rhot = eq.density*eq.temperature
# radial integration weights for integrals of rho * T
gi = GridInfo(dirname + '/grid_info', '')
nr = gi.nr
if sep_czrz:
    nr_cz = get_parameter(dirname, 'ncheby')[1]
    nr_rz = nr - nr_cz
    rw = gi.rweights
    rw_cz = np.copy(rw[:nr_cz])
    rw_rz = np.copy(rw[nr_cz:])
    rw_cz /= np.sum(rw_cz)
    rw_rz /= np.sum(rw_rz)

# possibly shift inte values to lie just above max ke
if plot_inte and subinte:
    min_inte = np.min(inte)
    max_ke = np.max(ke)
    diff = min_inte - max_ke
    buff = max_ke*0.05
    sub = diff - buff
    inte -= sub
    if sep_czrz:
        # integrate rho * T
        integ = np.sum(rw*rhot)
        S0 = sub/integ # entropy associated with subtraction
        integ_cz = np.sum(rw_cz*rhot[:nr_cz])
        integ_rz = np.sum(rw_rz*rhot[nr_cz:])
        inte_cz -= (S0*integ_cz)
        inte_rz -= (S0*integ_rz)
 
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

    # use units of the stellar luminosity
    lstar = get_lum(dirname)
    # m_leak represents leak in energy DENSITY (multiply by volume)
    ri, ro = np.min(gi.radius), np.max(gi.radius)
    shell_volume = 4./3.*np.pi*(ro**3. - ri**3.)
    dEdt = m_leak*shell_volume/lstar 
    mdEdt = mm_leak*shell_volume/lstar
    fdEdt = fm_leak*shell_volume/lstar

if sep_czrz:
    # Get energy densities (CZ and RZ separately)
    rke_cz = vals_cz[:, lut[402]]
    tke_cz = vals_cz[:, lut[403]]
    pke_cz = vals_cz[:, lut[404]]
    ke_cz = rke_cz + tke_cz + pke_cz

    frke_cz = vals_cz[:, lut[410]]
    ftke_cz = vals_cz[:, lut[411]]
    fpke_cz = vals_cz[:, lut[412]]
    fke_cz = frke_cz + ftke_cz + fpke_cz

    mrke_cz = rke_cz - frke_cz
    mtke_cz = tke_cz - ftke_cz
    mpke_cz = pke_cz - fpke_cz
    mke_cz = mrke_cz + mtke_cz + mpke_cz

    rke_rz = vals_rz[:, lut[402]]
    tke_rz = vals_rz[:, lut[403]]
    pke_rz = vals_rz[:, lut[404]]
    ke_rz = rke_rz + tke_rz + pke_rz

    frke_rz = vals_rz[:, lut[410]]
    ftke_rz = vals_rz[:, lut[411]]
    fpke_rz = vals_rz[:, lut[412]]
    fke_rz = frke_rz + ftke_rz + fpke_rz

    mrke_rz = rke_rz - frke_rz
    mtke_rz = tke_rz - ftke_rz
    mpke_rz = pke_rz - fpke_rz
    mke_rz = mrke_rz + mtke_rz + mpke_rz

    # Get the magnetic energies if they are available
    if magnetism:
        rme_cz = vals_cz[:, lut[1102]]
        tme_cz = vals_cz[:, lut[1103]]
        pme_cz = vals_cz[:, lut[1104]]
        me_cz = rme_cz + tme_cz + pme_cz

        frme_cz = vals_cz[:, lut[1110]]
        ftme_cz = vals_cz[:, lut[1111]]
        fpme_cz = vals_cz[:, lut[1112]]
        fme_cz = frme_cz + ftme_cz + fpme_cz

        mrme_cz = rme_cz - frme_cz
        mtme_cz = tme_cz - ftme_cz
        mpme_cz = pme_cz - fpme_cz
        mme_cz = mrme_cz + mtme_cz + mpme_cz

        rme_rz = vals_rz[:, lut[1102]]
        tme_rz = vals_rz[:, lut[1103]]
        pme_rz = vals_rz[:, lut[1104]]
        me_rz = rme_rz + tme_rz + pme_rz

        frme_rz = vals_rz[:, lut[1110]]
        ftme_rz = vals_rz[:, lut[1111]]
        fpme_rz = vals_rz[:, lut[1112]]
        fme_rz = frme_rz + ftme_rz + fpme_rz

        mrme_rz = rme_rz - frme_rz
        mtme_rz = tme_rz - ftme_rz
        mpme_rz = pme_rz - fpme_rz
        mme_rz = mrme_rz + mtme_rz + mpme_rz

    # Separated internal energies already taken care of

    # see if we desire total energies
    if plot_tote:
        tote_cz = ke_cz + inte_cz
        tote_rz = ke_rz + inte_rz
        ftote_cz = np.copy(fke_cz)
        ftote_rz = np.copy(fke_rz)
        mtote_cz = mke_cz + inte_cz
        mtote_rz = mke_rz + inte_rz

        if magnetism:
            tote_cz += me_cz

        # zone volumes
        rbcz = gi.radius[nr_cz - 1]
        shell_volume_cz = 4./3.*np.pi*(ro**3. - rbcz**3.)
        shell_volume_rz = 4./3.*np.pi*(rbcz**3. - ri**3.)

        # Compute leaks (CZ)
        tote_leak_cz = np.copy(tote_cz[it_leak:])
        mtote_leak_cz = np.copy(mtote_cz[it_leak:])
        ftote_leak_cz = np.copy(ftote_cz[it_leak:])
        
        m_leak_cz, b_leak_cz = np.polyfit(times_leak, tote_leak_cz, 1)
        mm_leak_cz, mb_leak_cz = np.polyfit(times_leak, mtote_leak_cz, 1)
        fm_leak_cz, fb_leak_cz = np.polyfit(times_leak, ftote_leak_cz, 1)

        dEdt_cz = m_leak_cz*shell_volume_cz/lstar 
        mdEdt_cz = mm_leak_cz*shell_volume_cz/lstar
        fdEdt_cz = fm_leak_cz*shell_volume_cz/lstar

        # Compute leaks (RZ)
        tote_leak_rz = np.copy(tote_rz[it_leak:])
        mtote_leak_rz = np.copy(mtote_rz[it_leak:])
        ftote_leak_rz = np.copy(ftote_rz[it_leak:])
        
        m_leak_rz, b_leak_rz = np.polyfit(times_leak, tote_leak_rz, 1)
        mm_leak_rz, mb_leak_rz = np.polyfit(times_leak, mtote_leak_rz, 1)
        fm_leak_rz, fb_leak_rz = np.polyfit(times_leak, ftote_leak_rz, 1)

        dEdt_rz = m_leak_rz*shell_volume_rz/lstar 
        mdEdt_rz = mm_leak_rz*shell_volume_rz/lstar
        fdEdt_rz = fm_leak_rz*shell_volume_rz/lstar

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

# See if y-axis should be on log scale (must do this before setting limits)
# Make all axes use scientific notation (except for y if ylog=True)
if ylog:
    for ax in axs.flatten():
        ax.set_yscale('log')
        ax.ticklabel_format(axis='x', scilimits=(-3,4), useMathText=True)
else:
    for ax in axs.flatten():
        ax.ticklabel_format(scilimits = (-3,4), useMathText=True)


# ===============TOTAL ENERGIES (IN BOTH ZONES)========================
# plot total KE
axs[0,0].plot(xaxis, ke, 'm', linewidth=lw, label=r'$\rm{E_{tot}}$')
axs[0,0].plot(xaxis, rke, 'r', linewidth=lw, label=r'$\rm{E_r}$')
axs[0,0].plot(xaxis, tke, 'g', linewidth=lw, label=r'$\rm{E_\theta}$')
axs[0,0].plot(xaxis, pke, 'b', linewidth=lw, label=r'$\rm{E_\phi}$')

# plot mean KE
axs[1,0].plot(xaxis, mke, 'm', linewidth=lw)
axs[1,0].plot(xaxis, mrke, 'r', linewidth=lw)
axs[1,0].plot(xaxis, mtke, 'g', linewidth=lw)
axs[1,0].plot(xaxis, mpke, 'b', linewidth=lw)

# plot fluc KE
axs[2,0].plot(xaxis, fke, 'm', linewidth=lw)
axs[2,0].plot(xaxis, frke, 'r', linewidth=lw)
axs[2,0].plot(xaxis, ftke, 'g', linewidth=lw)
axs[2,0].plot(xaxis, fpke, 'b', linewidth=lw)

# If magnetic, plot magnetic energies!
if magnetism:
    # plot total ME
    axs[0,0].plot(xaxis, me, 'm--', linewidth=lw)
    axs[0,0].plot(xaxis, rme, 'r--', linewidth=lw)
    axs[0,0].plot(xaxis, tme, 'g--', linewidth=lw)
    axs[0,0].plot(xaxis, pme, 'b--', linewidth=lw)

    # plot mean ME
    axs[1,0].plot(xaxis, mme, 'm--', linewidth=lw)
    axs[1,0].plot(xaxis, mrme, 'r--', linewidth=lw)
    axs[1,0].plot(xaxis, mtme, 'g--', linewidth=lw)
    axs[1,0].plot(xaxis, mpme, 'b--', linewidth=lw)

    # plot fluc ME
    axs[2,0].plot(xaxis, fme, 'm--', linewidth=lw)
    axs[2,0].plot(xaxis, frme, 'r--', linewidth=lw)
    axs[2,0].plot(xaxis, ftme, 'g--', linewidth=lw)
    axs[2,0].plot(xaxis, fpme, 'b--', linewidth=lw)

# See if various internal/total energies should be plotted
if plot_inte:
    axs[0,0].plot(xaxis, inte, 'c', linewidth=lw, label=inte_label)
    axs[1,0].plot(xaxis, inte, 'c', linewidth=lw)
if plot_tote:
    axs[0,0].plot(xaxis, tote, 'k', linewidth=lw, label=r'$\rm{TE}$')
    axs[1,0].plot(xaxis, mtote, 'k', linewidth=lw)
    axs[2,0].plot(xaxis, ftote, 'k', linewidth=lw)

    # Label energy leaks with rate of leak and dashed red line
    # full
    leak_label = r'$\rm{(\overline{dE/dt})/L_* = %09.3e}$'
    leak_label_x = 0.99
    leak_label_y = 0.01
    ax = axs[0,0]
    ax.text(leak_label_x, leak_label_y, leak_label %dEdt,\
            va='bottom', ha='right', transform=ax.transAxes)
    ax.plot(xaxis[it_leak:], m_leak*times_leak + b_leak, 'r--') 
    # mean
    ax = axs[1,0]
    ax.text(leak_label_x, leak_label_y, leak_label %mdEdt,\
            va='bottom', ha='right', transform=ax.transAxes)
    ax.plot(xaxis[it_leak:], mm_leak*times_leak + mb_leak, 'r--') 
    # fluc
    ax = axs[2,0]
    ax.text(leak_label_x, leak_label_y, leak_label %fdEdt,\
            va='bottom', ha='right', transform=ax.transAxes)
    ax.plot(xaxis[it_leak:], fm_leak*times_leak + fb_leak, 'r--') 

# put a legend on the upper left axis
axs[0,0].legend(loc='lower left', ncol=3, fontsize=8, columnspacing=1)

# set y-axis limits
# take into account buffer for legend and leak labels
def get_minmax(ax, ylog=False, withleg=False,\
        dyn=True, mes=None):
    if withleg: # extra room for legend
        buff_frac = 0.3
    else:
        buff_frac = 0.1
    ymin_current, ymax_current = ax.get_ylim()
    if ylog:
        if dyn:
            yratio = ymax_current/ymin_current
            ymin_new = ymin_current/(yratio**buff_frac)
        else:
            # Get minimum mag energy over last half of dynamo
            rme_loc, tme_loc, pme_loc = mes
            ntimes_loc = len(rme_loc)
            it_half = ntimes_loc//2
            minme = min(np.min(rme_loc[it_half]),\
                    np.min(rme_loc[it_half]), np.min(rme_loc[it_half]))
            yratio = ymax_current/minme
            ymin_new = minme/(yratio**buff_frac)
    else:
        ydiff = ymax_current - ymin_current
        ymin_new = ymin_current - buff_frac*ydiff
    return ymin_new, ymax_current

if minmax is None:
    if magnetism and not dyn: mes = rme, tme, pme
    minmax = get_minmax(axs[0,0], ylog=ylog, withleg=True, dyn=dyn, mes=mes)
if not ymin is None:
    minmax = ymin, minmax[1]
if not ymax is None:
    minmax = minmax[0], ymax
axs[0,0].set_ylim((minmax[0], minmax[1]))

# mean flows
if mminmax is None:
    if magnetism and not dyn: mes = mrme, mtme, mpme
    mminmax = get_minmax(axs[1,0], ylog=ylog, withleg=True,\
            dyn=dyn, mes=mes)
if not mymin is None:
    mminmax = mymin, mminmax[1]
if not mymax is None:
    mminmax = mminmax[0], mymax
axs[1,0].set_ylim((mminmax[0], mminmax[1]))

# fluc flows
if fminmax is None:
    if magnetism and not dyn: mes = frme, ftme, fpme
    fminmax = get_minmax(axs[2,0], ylog=ylog, withleg=True,\
            dyn=dyn, mes=mes)
if not fymin is None:
    fminmax = fymin, fminmax[1]
if not fymax is None:
    fminmax = fminmax[0], fymax
axs[2,0].set_ylim((fminmax[0], fminmax[1]))

# NOW POSSIBLY, SEPARATE ENERGIES IN CZ AND RZ
if sep_czrz:
    # =============CONVECTION ZONE=================
    # plot total KE
    axs[0,1].plot(xaxis, ke_cz, 'm', linewidth=lw)
    axs[0,1].plot(xaxis, rke_cz, 'r', linewidth=lw)
    axs[0,1].plot(xaxis, tke_cz, 'g', linewidth=lw)
    axs[0,1].plot(xaxis, pke_cz, 'b', linewidth=lw)

    # plot mean KE
    axs[1,1].plot(xaxis, mke_cz, 'm', linewidth=lw)
    axs[1,1].plot(xaxis, mrke_cz, 'r', linewidth=lw)
    axs[1,1].plot(xaxis, mtke_cz, 'g', linewidth=lw)
    axs[1,1].plot(xaxis, mpke_cz, 'b', linewidth=lw)

    # plot fluc KE
    axs[2,1].plot(xaxis, fke_cz, 'm', linewidth=lw)
    axs[2,1].plot(xaxis, frke_cz, 'r', linewidth=lw)
    axs[2,1].plot(xaxis, ftke_cz, 'g', linewidth=lw)
    axs[2,1].plot(xaxis, fpke_cz, 'b', linewidth=lw)

    # If magnetic, plot magnetic energies!
    if magnetism:
        # plot total ME
        axs[0,1].plot(xaxis, me_cz, 'm--', linewidth=lw)
        axs[0,1].plot(xaxis, rme_cz, 'r--', linewidth=lw)
        axs[0,1].plot(xaxis, tme_cz, 'g--', linewidth=lw)
        axs[0,1].plot(xaxis, pme_cz, 'b--', linewidth=lw)

        # plot mean ME
        axs[1,1].plot(xaxis, mme_cz, 'm--', linewidth=lw)
        axs[1,1].plot(xaxis, mrme_cz, 'r--', linewidth=lw)
        axs[1,1].plot(xaxis, mtme_cz, 'g--', linewidth=lw)
        axs[1,1].plot(xaxis, mpme_cz, 'b--', linewidth=lw)

        # plot fluc ME
        axs[2,1].plot(xaxis, fme_cz, 'k--', linewidth=lw)
        axs[2,1].plot(xaxis, frme_cz, 'r--', linewidth=lw)
        axs[2,1].plot(xaxis, ftme_cz, 'g--', linewidth=lw)
        axs[2,1].plot(xaxis, fpme_cz, 'b--', linewidth=lw)

    # See if various internal/total energies should be plotted
    if plot_inte:
        axs[0,1].plot(xaxis, inte_cz, 'c', linewidth=lw)
        axs[1,1].plot(xaxis, inte_cz, 'c', linewidth=lw)
    if plot_tote:
        axs[0,1].plot(xaxis, tote_cz, 'k', linewidth=lw)
        axs[1,1].plot(xaxis, mtote_cz, 'k', linewidth=lw)
        axs[2,1].plot(xaxis, ftote_cz, 'k', linewidth=lw)

        # Label energy leaks with rate of leak and dashed line
        # full
        ax = axs[0,1]
        ax.text(leak_label_x, leak_label_y, leak_label %dEdt_cz,\
                va='bottom', ha='right', transform=ax.transAxes)
        ax.plot(xaxis[it_leak:], m_leak_cz*times_leak + b_leak_cz, 'r--') 
        # mean
        ax = axs[1,1]
        ax.text(leak_label_x, leak_label_y, leak_label %mdEdt_cz,\
                va='bottom', ha='right', transform=ax.transAxes)
        ax.plot(xaxis[it_leak:], mm_leak_cz*times_leak + mb_leak_cz, 'r--') 
        # fluc
        ax = axs[2,1]
        ax.text(leak_label_x, leak_label_y, leak_label %fdEdt_cz,\
                va='bottom', ha='right', transform=ax.transAxes)
        ax.plot(xaxis[it_leak:], fm_leak_cz*times_leak + fb_leak_cz, 'r--') 

    # set y-axis limits
    # take into account buffer for legend and leak labels
    # full flows
    if minmax_cz is None:
        if magnetism and not dyn: mes = rme_cz, tme_cz, pme_cz
        minmax_cz = get_minmax(axs[0,1], ylog=ylog, withleg=True,\
                dyn=dyn, mes=mes)
    if not ymin_cz is None:
        minmax_cz = ymin_cz, minmax_cz[1]
    if not ymax_cz is None:
        minmax_cz = minmax_cz[0], ymax_cz
    axs[0,1].set_ylim((minmax_cz[0], minmax_cz[1]))
    # mean flows
    if mminmax_cz is None:
        if magnetism and not dyn: mes = mrme_cz, mtme_cz, mpme_cz
        mminmax_cz = get_minmax(axs[1,1], ylog=ylog, withleg=True,\
                dyn=dyn, mes=mes)
    if not mymin_cz is None:
        mminmax_cz = mymin_cz, mminmax_cz[1]
    if not mymax_cz is None:
        mminmax_cz = mminmax_cz[0], mymax_cz
    axs[1,1].set_ylim((mminmax_cz[0], mminmax_cz[1]))
    # fluc flows
    if fminmax_cz is None:
        if magnetism and not dyn: mes = frme_cz, ftme_cz, fpme_cz
        fminmax_cz = get_minmax(axs[2,1], ylog=ylog, withleg=True,\
                dyn=dyn, mes=mes)
    if not fymin_cz is None:
        fminmax_cz = fymin_cz, fminmax_cz[1]
    if not fymax_cz is None:
        fminmax_cz = fminmax_cz[0], fymax_cz
    axs[2,1].set_ylim((fminmax_cz[0], fminmax_cz[1]))

    # ==================RADIATIVE ZONE=============
    # plot total KE
    axs[0,2].plot(xaxis, ke_rz, 'm', linewidth=lw)
    axs[0,2].plot(xaxis, rke_rz, 'r', linewidth=lw)
    axs[0,2].plot(xaxis, tke_rz, 'g', linewidth=lw)
    axs[0,2].plot(xaxis, pke_rz, 'b', linewidth=lw)

    # plot mean KE
    axs[1,2].plot(xaxis, mke_rz, 'm', linewidth=lw)
    axs[1,2].plot(xaxis, mrke_rz, 'r', linewidth=lw)
    axs[1,2].plot(xaxis, mtke_rz, 'g', linewidth=lw)
    axs[1,2].plot(xaxis, mpke_rz, 'b', linewidth=lw)

    # plot fluc KE
    axs[2,2].plot(xaxis, fke_rz, 'm', linewidth=lw)
    axs[2,2].plot(xaxis, frke_rz, 'r', linewidth=lw)
    axs[2,2].plot(xaxis, ftke_rz, 'g', linewidth=lw)
    axs[2,2].plot(xaxis, fpke_rz, 'b', linewidth=lw)

    # If magnetic, plot magnetic energies!
    if magnetism:
        # plot total ME
        axs[0,2].plot(xaxis, me_rz, 'm--', linewidth=lw)
        axs[0,2].plot(xaxis, rme_rz, 'r--', linewidth=lw)
        axs[0,2].plot(xaxis, tme_rz, 'g--', linewidth=lw)
        axs[0,2].plot(xaxis, pme_rz, 'b--', linewidth=lw)

        # plot mean ME
        axs[1,2].plot(xaxis, mme_rz, 'm--', linewidth=lw)
        axs[1,2].plot(xaxis, mrme_rz, 'r--', linewidth=lw)
        axs[1,2].plot(xaxis, mtme_rz, 'g--', linewidth=lw)
        axs[1,2].plot(xaxis, mpme_rz, 'b--', linewidth=lw)

        # plot fluc ME
        axs[2,2].plot(xaxis, fme_rz, 'k--', linewidth=lw)
        axs[2,2].plot(xaxis, frme_rz, 'r--', linewidth=lw)
        axs[2,2].plot(xaxis, ftme_rz, 'g--', linewidth=lw)
        axs[2,2].plot(xaxis, fpme_rz, 'b--', linewidth=lw)

    # See if various internal/total energies should be plotted
    if plot_inte:
        axs[0,2].plot(xaxis, inte_rz, 'c', linewidth=lw)
        axs[1,2].plot(xaxis, inte_rz, 'c', linewidth=lw)
    if plot_tote:
        axs[0,2].plot(xaxis, tote_rz, 'k', linewidth=lw)
        axs[1,2].plot(xaxis, mtote_rz, 'k', linewidth=lw)
        axs[2,2].plot(xaxis, ftote_rz, 'k', linewidth=lw)

        # Label energy leaks with rate of leak and dashed line
        # full
        ax = axs[0,2]
        ax.text(leak_label_x, leak_label_y, leak_label %dEdt_rz,\
                va='bottom', ha='right', transform=ax.transAxes)
        ax.plot(xaxis[it_leak:], m_leak_rz*times_leak + b_leak_rz, 'r--') 
        # mean
        ax = axs[1,2]
        ax.text(leak_label_x, leak_label_y, leak_label %mdEdt_rz,\
                va='bottom', ha='right', transform=ax.transAxes)
        ax.plot(xaxis[it_leak:], mm_leak_rz*times_leak + mb_leak_rz, 'r--') 
        # fluc
        ax = axs[2,2]
        ax.text(leak_label_x, leak_label_y, leak_label %fdEdt_rz,\
                va='bottom', ha='right', transform=ax.transAxes)
        ax.plot(xaxis[it_leak:], fm_leak_rz*times_leak + fb_leak_rz, 'r--') 

    # set y-axis limits
    # take into account buffer for legend and leak labels
    # full flows
    if minmax_rz is None:
        if magnetism and not dyn: mes = rme_rz, tme_rz, pme_rz
        minmax_rz = get_minmax(axs[0,2], ylog=ylog, withleg=True,\
                dyn=dyn, mes=mes)
    if not ymin_rz is None:
        minmax_rz = ymin_rz, minmax_rz[1]
    if not ymax_rz is None:
        minmax_rz = minmax_rz[0], ymax_rz
    axs[0,2].set_ylim((minmax_rz[0], minmax_rz[1]))
    # mean flows
    if mminmax_rz is None:
        if magnetism and not dyn: mes = mrme_rz, mtme_rz, mpme_rz
        mminmax_rz = get_minmax(axs[1,2], ylog=ylog, withleg=True,\
                dyn=dyn, mes=mes)
    if not mymin_rz is None:
        mminmax_rz = mymin_rz, mminmax_rz[1]
    if not mymax_rz is None:
        mminmax_rz = mminmax_rz[0], mymax_rz
    axs[1,2].set_ylim((mminmax_rz[0], mminmax_rz[1]))
    # fluc flows
    if fminmax_rz is None:
        if magnetism and not dyn: mes = frme_rz, ftme_rz, fpme_rz
        fminmax_rz = get_minmax(axs[2,2], ylog=ylog, withleg=True,\
                dyn=dyn, mes=mes)
    if not fymin_rz is None:
        fminmax_rz = fymin_rz, fminmax_rz[1]
    if not fymax_rz is None:
        fminmax_rz = fminmax_rz[0], fymax_rz
    axs[2,2].set_ylim((fminmax_rz[0], fminmax_rz[1]))

# Set some parameters defining all subplots
# x limits and label
axs[0,0].set_xlim((xminmax[0], xminmax[1]))
if xiter:
    axs[2,0].set_xlabel('iteration #')
else:
    if rotation:
        axs[2,0].set_xlabel(r'$\rm{t\ (P_{rot})}$')
    else:
        axs[2,0].set_xlabel(r'$\rm{t\ (T_{diff})}$')

# y labels
fs = 14
axs[0,0].set_ylabel('full energy', fontsize=fs)
axs[1,0].set_ylabel('mean energy', fontsize=fs)
axs[2,0].set_ylabel('fluc energy', fontsize=fs)

# Make titles
if ncol == 1:
    icol_mid = 0
elif ncol == 3:
    icol_mid = 1
titles = ["ALL ZONES", "CZ", "RZ"]
for icol in range(ncol):
    if icol == icol_mid:
        title = dirname_stripped + '\n\n ' + titles[icol]
    else:
        title = titles[icol]
    axs[0, icol].set_title(title)

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


#  set axis limits
if fminmax is None:
    fminmax = axs[2,0].get_ylim()
# Change JUST ymin or ymax, if desired
if not ymin is None:
    fminmax = ymin, fminmax[1]
if not ymax is None:
    fminmax = fminmax[0], ymax
axs[2,0].set_ylim((fminmax[0], fminmax[1]))

# Get ticks everywhere
for ax in axs.flatten():
    plt.sca(ax)
    plt.minorticks_on()
    plt.tick_params(top=True, right=True, direction='in', which='both')

# Space the subplots to make them look pretty
plt.tight_layout()
#plt.subplots_adjust(left=0.15, bottom=0.08, top=0.85, wspace=0.4)

if savename is None:
    savename = dirname_stripped + '_etrace_' + str(iter1).zfill(8) + '_' + str(iter2).zfill(8) + tag + '.png'

# Save the plot
print ('Saving the etrace plot at ' + plotdir + savename)
plt.savefig(plotdir + savename, dpi=300)

# Show the plot
plt.show()
