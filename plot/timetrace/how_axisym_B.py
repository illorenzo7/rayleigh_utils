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
from common import *
from rayleigh_diagnostics import GridInfo

# Get the run directory on which to perform the analysis
dirname = sys.argv[1]

# Data and plot directories
datadir = dirname + '/data/'
dirname_stripped = strip_dirname(dirname)

# mark times
plottimes = None

# domain bounds
ncheby, domain_bounds = get_domain_bounds(dirname)
ri = np.min(domain_bounds)
ro = np.max(domain_bounds)
d = ro - ri

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

# bunch of minmax options
# use -minmax [val1] [val2] to set both minmax values (combined zones)
# use -min [val] or -max [val] to set just the min or max (combined zones)
# append cz or rz for values in cz or rz
# prepend a, m, f, to set all min/max vals (full, mean, fluc), m for mean,
# f for fluc, 
# or aa for full, mean, and fluc for each zone (combined, RZ, CZ)
# example:
# aminmaxcz 0 5e6 sets the full, mean, fluc bounds in the CZ to 0, 5e6
# aaminmax 0 5e6 does the same for combined zones, CZ, and RZ

aminmax = None
amin = None
amax = None
aminmax_cz = None
amin_cz = None
amax_cz = None
aminmax_rz = None
amin_rz = None
amax_rz = None
aaminmax = None
aamin = None
aamax = None

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
    elif arg == '-times':
        strings = args[i+1].split()
        plottimes = []
        for string in strings:
            plottimes.append(float(string))
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

    elif arg == '-aminmax':
        aminmax = float(args[i+1]), float(args[i+2])
    elif arg == '-amin':
        amin = float(args[i+1])
    elif arg == '-amax':
        amax = float(args[i+1])
    elif arg == '-aminmaxcz':
        aminmax_cz = float(args[i+1]), float(args[i+2])
    elif arg == '-amincz':
        amin_cz = float(args[i+1])
    elif arg == '-amaxcz':
        amax_cz = float(args[i+1])
    elif arg == '-aminmaxrz':
        aminmax_rz = float(args[i+1]), float(args[i+2])
    elif arg == '-aminrz':
        amin_rz = float(args[i+1])
    elif arg == '-amaxrz':
        amax_rz = float(args[i+1])
    elif arg == '-aaminmax':
        aaminmax = float(args[i+1]), float(args[i+2])
    elif arg == '-aamin':
        aamin = float(args[i+1])
    elif arg == '-aamax':
        aamax = float(args[i+1])

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

vals = vals[ixmin:ixmax+1, :]
if sep_czrz:
    vals_cz = vals_cz[ixmin:ixmax+1, :]
    vals_rz = vals_rz[ixmin:ixmax+1, :]

print ("ntot = %i" %ntot)
print ("before thin_data: len(xaxis) = %i" %len(xaxis))
xaxis = thin_data(xaxis, ntot)
times = thin_data(times, ntot)
iters = thin_data(iters, ntot)
vals = thin_data(vals, ntot)
if sep_czrz:
    vals_cz = thin_data(vals_cz, ntot)
    vals_rz = thin_data(vals_rz, ntot)
print ("after thin_data: len(xaxis) = %i" %len(xaxis))

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
elif inte_gtr2 or ((plot_inte or plot_tote) and sep_czrz): 
    # inte not from trace_G_Avgs
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
rr = gi.radius
nr = gi.nr
if sep_czrz:
    nr_cz = get_parameter(dirname, 'ncheby')[1]
    rbcz = rr[nr_cz - 1]
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
    ri, ro = np.min(rr), np.max(rr)
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
fig, axs = plt.subplots(1, ncol, figsize=(5.*ncol, 5.),\
        sharex=True)
# first, expand it row wise
axs = np.expand_dims(axs, 0)
if ncol == 1: # need the axis array to consistently be doubly indexed
    axs = np.expand_dims(axs, 1)

# Make thin lines to see structure of variation for ME
lw = 0.5
lw_ke = 1. # bit thicker for KE to differentiate between ME

# See if y-axis should be on log scale (must do this before setting limits)
# Make all axes use scientific notation (except for y if ylog=True)
if ylog:
    for ax in axs.flatten():
        ax.set_yscale('log')
        ax.ticklabel_format(axis='x', scilimits=(-3,4), useMathText=True)
else:
    for ax in axs.flatten():
        ax.ticklabel_format(scilimits = (-3,4), useMathText=True)

# ===============axisym ratio (IN BOTH ZONES)========================
axs[0,0].plot(xaxis, np.sqrt(mme/fme), 'm', linewidth=lw, label=r'$\rm{ME_{tot}}$')
axs[0,0].plot(xaxis, np.sqrt(mrme/frme), 'r', linewidth=lw, label=r'$\rm{ME_r}$')
axs[0,0].plot(xaxis, np.sqrt(mtme/ftme), 'g', linewidth=lw, label=r'$\rm{ME_\theta}$')
axs[0,0].plot(xaxis, np.sqrt(mpme/fpme), 'b', linewidth=lw, label=r'$\rm{ME_\phi}$')

# NOW POSSIBLY, SEPARATE ENERGIES IN CZ AND RZ
if sep_czrz:
    # ===============axisym ratio (CONVECTION ZONE)========================
    axs[0,1].plot(xaxis, np.sqrt(mme_cz/fme_cz), 'm', linewidth=lw, label=r'$\rm{ME_{tot}}$')
    axs[0,1].plot(xaxis, np.sqrt(mrme_cz/frme_cz), 'r', linewidth=lw, label=r'$\rm{ME_r}$')
    axs[0,1].plot(xaxis, np.sqrt(mtme_cz/ftme_cz), 'g', linewidth=lw, label=r'$\rm{ME_\theta}$')
    axs[0,1].plot(xaxis, np.sqrt(mpme_cz/fpme_cz), 'b', linewidth=lw, label=r'$\rm{ME_\phi}$')

    # ===============axisym ratio (RADIATION ZONE)========================
    axs[0,2].plot(xaxis, np.sqrt(mme_rz/fme_rz), 'm', linewidth=lw, label=r'$\rm{ME_{tot}}$')
    axs[0,2].plot(xaxis, np.sqrt(mrme_rz/frme_rz), 'r', linewidth=lw, label=r'$\rm{ME_r}$')
    axs[0,2].plot(xaxis, np.sqrt(mtme_rz/ftme_rz), 'g', linewidth=lw, label=r'$\rm{ME_\theta}$')
    axs[0,2].plot(xaxis, np.sqrt(mpme_rz/fpme_rz), 'b', linewidth=lw, label=r'$\rm{ME_\phi}$')

# Set some parameters defining all subplots
# x limits and label
axs[0,0].set_xlim((xminmax[0], xminmax[1]))
if xiter:
    axs[0,0].set_xlabel('iteration #')
else:
    if rotation:
        axs[0,0].set_xlabel(r'$\rm{t\ (P_{rot})}$')
    else:
        axs[0,0].set_xlabel(r'$\rm{t\ (T_{diff})}$')

# y labels
fs = 14
axs[0,0].set_ylabel('sqrt(ratio): m = 0 to m != 0', fontsize=fs)

# Make titles
if ncol == 1:
    icol_mid = 0
elif ncol == 3:
    icol_mid = 1
titles = ["ALL ZONES", "CZ", "RZ"]
for icol in range(ncol):
    if icol == icol_mid:
        if sep_czrz:
            title = dirname_stripped + '\n' r'$r_{bcz} =\ $' +\
                sci_format(rbcz, 2) + ' cm\n' + titles[icol]
        else:
            title = dirname_stripped + '\n\n' + titles[icol]
    else:
        title = titles[icol]
    axs[0, icol].set_title(title)

if not plottimes is None:
    for ax in axs.flatten():
        y1, y2 = ax.get_ylim()
        yvals = np.linspace(y1, y2, 100)
        for time in plottimes:
            ax.plot(time + np.zeros(100), yvals,'k--')

# Get ticks everywhere
for ax in axs.flatten():
    plt.sca(ax)
    plt.minorticks_on()
    plt.tick_params(top=True, right=True, direction='in', which='both')

# Space the subplots to make them look pretty
plt.tight_layout()
#plt.subplots_adjust(left=0.15, bottom=0.08, top=0.85, wspace=0.4)

if savename is None:
    savename = dirname_stripped + '_how_axisym_B_' + str(iter1).zfill(8) + '_' + str(iter2).zfill(8) + tag + '.png'

# Save the plot
print ('Saving the etrace plot at ' + plotdir + savename)
plt.savefig(plotdir + savename, dpi=300)

# Show the plot
plt.show()
