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
from common import *
from cla_util import *
from rayleigh_diagnostics import GridInfo

# Get the run directory on which to perform the analysis
dirname = sys.argv[1]

# Data and plot directories
datadir = dirname + '/data/'
dirname_stripped = strip_dirname(dirname)

# Read command-line arguments (CLAs)
args = sys.argv[2:]
clas = read_clas(dirname, args)
minmax_user = clas['minmax']
ymin = clas['ymin']
ymax = clas['ymax']

# See if magnetism is "on"
try:
    magnetism = get_parameter(dirname, 'magnetism')
except:
    magnetism = False # if magnetism wasn't specified, it must be "off"

# SPECIFIC ARGS for etrace:
coords = read_cla_arbitrary(args, 'coords', None, 'int') 
if not coords is None:
    numpanels = len(coords)//2
    acopy = np.copy(coords)
    coords = []
    for i in range(numpanels):
        coords.append((acopy[2*i], acopy[2*i + 1]))

# choose one panel (or panels) to apply minmax to
ntot = read_cla_arbitrary(args, 'ntot', 500, 'int') 
# total no. points to use (thin out really long arrays)
xiter = read_cla_arbitrary(args, 'xiter', False)
plottimes = read_cla_arbitrary(args, 'times')
from0 = read_cla_arbitrary(args, 'from0', False)
ylog = read_cla_arbitrary(args, 'log', False)
nodyn = read_cla_arbitrary(args, 'nodyn', False)
# by default don't adjust the min val to ignore super small 
# magnetic energies during dynamo growth when ylog=True (i.e., plot
# the dynamo growth phase by default)
# to change, use --nodyn / --dynfrac [val=0.5] to ignore the ME values over
# the last dynfrac of the simulation
dyn_frac = read_cla_arbitrary(args, 'dynfrac', 0.5)
leak_frac = read_cla_arbitrary(args, 'leakfrac', 0.25)

plot_inte = read_cla_arbitrary(args, 'inte', False)
plot_tote = read_cla_arbitrary(args, 'tote', False)
subinte = True # by default shift the internal energies by a constant
    # so they aren't so huge
if read_cla_arbitrary(args, 'nosub', False):
    subinte = False
inte_subt = read_cla_arbitrary(args, 'subt', False)
# subtracts top value of S for inte
inte_subb = read_cla_arbitrary(args, 'subb', False)
# subtracts bot value of S for inte
sep_czrz = read_cla_arbitrary(args, 'czrz', False)
# plots two more columns with energies in CZ and RZ separately 
inte_gtr2 = read_cla_arbitrary(args, 'gtr2', False)
# plots just the regular inte but using the 
# G_Avgs_trace_2dom file 
if inte_subb or inte_subt or inte_gtr2: 
    plot_inte = True

# Might need to use 2dom trace instead of regular trace
the_file = clas['the_file']
if inte_gtr2 or inte_subt or inte_subb or sep_czrz:
    if the_file is None:
        the_file = get_widest_range_file(datadir, 'G_Avgs_trace_2dom')
    print ('Using 2dom trace from ' + datadir + the_file)
    di = get_dict(datadir + the_file) 
    vals_gav = di['vals']
    vals_cz = di['vals_cz']
    vals_rz = di['vals_rz']
else: 
    if the_file is None:
        the_file = get_widest_range_file(datadir, 'G_Avgs_trace')
    print ('Using G_Avgs trace from ' + datadir + the_file)
    di = get_dict(datadir + the_file)
    vals_gav = di['vals']
# get lut, times, iters no matter the file type
lut = di['lut']
times = di['times']
iters = di['iters']

# get the x axis
time_unit, time_label, rotation = get_time_unit(dirname)
if not xiter:
    xaxis = times/time_unit
else:
    xaxis = iters

# set xminmax if not set by user
xminmax = clas['xminmax']
xmin = clas['xmin']
xmax = clas['xmax']
if xminmax is None:
    # set xmin possibly
    if xmin is None:
        if from0:
            xmin = 0.0
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
tmin, tmax = times[0], times[-1]

vals_gav = vals_gav[ixmin:ixmax+1, :]
if sep_czrz:
    vals_cz = vals_cz[ixmin:ixmax+1, :]
    vals_rz = vals_rz[ixmin:ixmax+1, :]

# Check what internal energies we might need (sometimes CZ/RZ separated ones
# too just in case) these will always be available for some options)
if inte_subt:
    inte_gav = vals_gav[:, lut[4001]]
    inte_cz = vals_cz[:, lut[4001]]
    inte_rz = vals_rz[:, lut[4001]]
    print("Got SUBT internal energy trace from trace_2dom_G_Avgs")
    inte_label = "IE SUBT"
elif inte_subb:
    inte_gav = vals_gav[:, lut[4002]]
    inte_cz = vals_cz[:, lut[4002]]
    inte_rz = vals_rz[:, lut[4002]]
    print("Got SUBB internal energy trace from trace_2dom_G_Avgs")
    inte_label = "IE SUBB"
elif inte_gtr2 or ((plot_inte or plot_tote) and sep_czrz): 
    # inte not from trace_G_Avgs
    inte_gav = vals_gav[:, lut[4000]]
    inte_cz = vals_cz[:, lut[4000]]
    inte_rz = vals_rz[:, lut[4000]]
    print("Got internal energy trace from trace_2dom_G_Avgs")
    inte_label = "IE W/ DRIFT"
elif plot_inte or plot_tote: 
    # try to get inte from trace_G_Avgs
    try:
        inte_gav = vals_gav[:, lut[701]]
        print("Got internal energy trace from trace_G_Avgs")
        inte_label = "IE W/ DRIFT"
    except:
        # if not available in trace_G_Avgs, set inte to zero 
        print ("Internal energy trace not available; setting to 0")
        inte_gav = np.zeros_like(times)
        inte_label = "IE NOT FOUND"

# Get rho*T and integrate it, + get volume
eq = get_eq(dirname)
rhot = eq.density*eq.temperature
integs = integrate_in_r(rhot, dirname)
volumes = get_volumes(dirname)
# get luminosity
lstar = get_lum(dirname)

# deal with x axis
print ("ntot = %i" %ntot)
print ("before thin_data: len(xaxis) = %i" %len(xaxis))
xaxis = thin_data(xaxis, ntot)
times = thin_data(times, ntot)
iters = thin_data(iters, ntot)
vals_gav = thin_data(vals_gav, ntot)
if plot_inte or plot_tote:
    inte_gav = thin_data(inte_gav, ntot)
if sep_czrz:
    vals_cz = thin_data(vals_cz, ntot)
    vals_rz = thin_data(vals_rz, ntot)
    if plot_inte or plot_tote:
        inte_cz = thin_data(inte_cz, ntot)
        inte_rz = thin_data(inte_rz, ntot)

print ("after thin_data: len(xaxis) = %i" %len(xaxis))

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

# start making plots
# loop over tot, fluc, mean and the different domains
vals_list = [vals_gav]
if sep_czrz:
    vals_list.append(vals_rz)
    vals_list.append(vals_cz)

if plot_inte or plot_tote:
    inte_list = [inte_gav]
    if sep_czrz:
        inte_list.append(inte_rz)
        inte_list.append(inte_cz)


# specialized function to get axis limits
# take into account buffer for legend and leak labels
def get_minmax(ax, ylog=False, withleg=False, nodyn=False, mes=None):
    ymin, ymax = ax.get_ylim()
    # possibly adjust ymin
    if nodyn:
        # Get minimum mag energy over last dyn_frac of dynamo
        rme_loc, tme_loc, pme_loc = mes
        ntimes_loc = len(rme_loc)
        it_frac = int(ntimes_loc*(1.0 - dyn_frac))
        minme = min(np.min(rme_loc[it_frac:]),\
                np.min(rme_loc[it_frac:]), np.min(rme_loc[it_frac:]))
        ymin = minme
        buff = 0.05
    if withleg:
        buff = 0.3

    # only need to adjust things for nodyn or withleg
    if nodyn or withleg:
        if ylog:
            yratio = ymax/ymin
            ymin = ymin/(yratio**buff)
        else:
            ydiff = ymax - ymin
            ymin = ymin - buff*ydiff
    return ymin, ymax

for irow in range(3): # tot, fluc, mean
    for icol in range(ncol):
        vals = vals_list[icol]
        if plot_inte or plot_tote:
            inte = inte_list[icol]
            if irow == 1: # fluc, no inte for fluctuating S'
                inte = np.zeros(len(xaxis))
        if irow == 0: # tot
            rke = vals[:, lut[402]]
            tke = vals[:, lut[403]]
            pke = vals[:, lut[404]]
        if irow == 1: # fluc
            rke = vals[:, lut[410]]
            tke = vals[:, lut[411]]
            pke = vals[:, lut[412]]
        if irow == 2: # mean
            rke = vals[:, lut[402]] - vals[:, lut[410]]
            tke = vals[:, lut[403]] - vals[:, lut[411]]
            pke = vals[:, lut[404]] - vals[:, lut[412]]
        ke = rke + tke + pke
        if magnetism:
            if irow == 0: # tot
                rme = vals[:, lut[1102]]
                tme = vals[:, lut[1103]]
                pme = vals[:, lut[1104]]
            if irow == 1: # fluc
                rme = vals[:, lut[1110]]
                tme = vals[:, lut[1111]]
                pme = vals[:, lut[1112]]
            if irow == 2: # mean
                rme = vals[:, lut[1102]] - vals[:, lut[1110]]
                tme = vals[:, lut[1103]] - vals[:, lut[1111]]
                pme = vals[:, lut[1104]] - vals[:, lut[1112]]
            me = rme + tme + pme

        # possibly shift inte values to lie just above max ke
        if plot_inte and subinte:
            min_inte = np.min(inte)
            max_ke = np.max(ke)
            diff = min_inte - max_ke
            buff = max_ke*0.05
            sub = diff - buff
            if icol == 0: # full domain
                inte -= sub*integs[icol]/integs[0]
         
        # Get total energy if desired
        if plot_tote:
            tote = ke + inte
            if magnetism:
                tote += me

            # Will need to compute "energy leaks":
            it_leak = np.argmin(np.abs((times - tmin) - (1.0 - leak_frac)*(tmax - tmin)))
            times_leak = np.copy(times[it_leak:])
            tote_leak = np.copy(tote[it_leak:])
           
            m_leak, b_leak = np.polyfit(times_leak, tote_leak, 1)

            # m_leak represents leak in energy DENSITY (multiply by volume)
            dEdt = m_leak*volumes[icol]/lstar 

        # make line plots
        ax = axs[irow, icol]
        ax.plot(xaxis, ke, color_order[0],\
                linewidth=lw_ke, label=r'$\rm{KE_{tot}}$')
        ax.plot(xaxis, rke, color_order[1],\
                linewidth=lw_ke, label=r'$\rm{KE_r}$')
        ax.plot(xaxis, tke, color_order[2],\
                linewidth=lw_ke, label=r'$\rm{KE_\theta}$')
        ax.plot(xaxis, pke, color_order[3],\
                linewidth=lw_ke, label=r'$\rm{KE_\phi}$')
        # If magnetic, plot magnetic energies!
        if magnetism:
            ax.plot(xaxis, me, color_order[0] + '--',\
                    linewidth=lw, label=r'$\rm{ME_{tot}}$')
            ax.plot(xaxis, rme, color_order[1] + '--',\
                    linewidth=lw, label=r'$\rm{ME_r}$')
            ax.plot(xaxis, tme, color_order[2] + '--',\
                    linewidth=lw, label=r'$\rm{ME_\theta}$')
            ax.plot(xaxis, pme, color_order[3] + '--',\
                    linewidth=lw, label=r'$\rm{ME_\phi}$')

        # See if various internal/total energies should be plotted
        if plot_inte:
            ax.plot(xaxis, inte, color_order[4],\
                    linewidth=lw_ke, label=inte_label)
        if plot_tote:
            ax.plot(xaxis, tote, color_order[5],\
                    linewidth=lw_ke, label=r'$\rm{TE}$')

            # Label energy leak with rate of leak and dashed red line
            leak_label = r'$\rm{(\overline{dE/dt})/L_* = %09.3e}$'
            leak_label_x = 0.99
            leak_label_y = 0.01
            ax.text(leak_label_x, leak_label_y,\
                    leak_label %dEdt,\
                    va='bottom', ha='right', transform=ax.transAxes)
            ax.plot(xaxis[it_leak:], m_leak*times_leak + b_leak, 'r--') 

        if irow == 0 and icol == 0: # put a legend on the upper left axis
            ax.legend(loc='lower left', ncol=3, fontsize=8, columnspacing=1)

        # set the y limits
        if coords is None:
            minmax = minmax_user
        else:
            if (irow, icol) in coords:
                minmax = minmax_user
            else: 
                minmax = None
        if minmax is None:
            if nodyn: 
                mes = rme, tme, pme
            else:
                mes = None
            minmax = get_minmax(ax, ylog=ylog, withleg=True, nodyn=nodyn, mes=mes)
        if not ymin is None:
            minmax = ymin, minmax[1]
        if not ymax is None:
            minmax = minmax[0], ymax
        ax.set_ylim((minmax[0], minmax[1]))

# Set some parameters defining all subplots
# x limits and label
if ncol == 1:
    icol_mid = 0
elif ncol == 3:
    icol_mid = 1
axs[0,icol_mid].set_xlim((xminmax[0], xminmax[1]))
if xiter:
    axs[2,icol_mid].set_xlabel('iteration #')
else:
    axs[2,icol_mid].set_xlabel(r'$t\ ($' + time_label + r'$)$')

# y labels
fs = 14
axs[0,0].set_ylabel('full energy', fontsize=fs)
axs[1,0].set_ylabel('fluc energy', fontsize=fs)
axs[2,0].set_ylabel('mean energy', fontsize=fs)

# Make titles
titles = ["ALL ZONES", "RZ", "CZ"]
for icol in range(ncol):
    if icol == icol_mid:
        title = dirname_stripped + '\n' + titles[icol]
    else:
        title = titles[icol]
    axs[0, icol].set_title(title)

# mark times if desired
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

# Save the plot
iter1, iter2 = get_iters_from_file(the_file)
# Tag the plot by whether or not the x axis is in "time" or "iteration"
tag = clas['tag']
if xiter and tag == '':
    tag = '_xiter'

savename = 'etrace' + tag + '-' + str(iter1).zfill(8) + '_' + str(iter2).zfill(8) + '.png'

# get plotdir and possibly create it
plotdir = clas['plotdir']
if plotdir is None:
    plotdir = dirname + '/plots/'
make_plotdir(plotdir)

print ('Saving the etrace plot at ' + plotdir + savename)
plt.savefig(plotdir + savename, dpi=300)

# Show the plot
plt.show()
