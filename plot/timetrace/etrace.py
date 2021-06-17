# Author: Loren Matilsky
# Date created: 06/08/2017
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
try:
    magnetism = get_parameter(dirname, 'magnetism')
except:
    magnetism = False # if magnetism wasn't specified, it must be "off"

# SPECIFIC ARGS for etrace:
kwargs_default = dict({'the_file': None, 'xminmax': None, 'xmin': None, 'xmax': None, 'minmax': None, 'min': None, 'max': None, 'coords': None, 'ntot': 500, 'xiter': False, 'log': False, 'nodyn': False, 'dynfrac': 0.5, 'xvals': np.array([]), 'czrz': False, 'inte': False})
# plots two more columns with energies in CZ and RZ separately 
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
sep_czrz = kwargs.czrz
plot_inte = kwargs.inte

# deal with coords (if user wants minmax to only apply to certain subplots)
if not coords is None:
    numpanels = len(coords)//2
    acopy = np.copy(coords)
    coords = []
    for i in range(numpanels):
        coords.append((acopy[2*i], acopy[2*i + 1]))

# Might need to use 2dom trace instead of regular trace
if the_file is None:
    if sep_czrz:
        the_file = get_widest_range_file(clas0['datadir'], 'G_Avgs_trace_2dom')
    else: 
        the_file = get_widest_range_file(clas0['datadir'], 'G_Avgs_trace')

print ('Getting data from ' + the_file)
di = get_dict(the_file)
vals_gav = di['vals']
if sep_czrz:
    vals_cz = di['vals_cz']
    vals_rz = di['vals_rz']
lut = di['lut']
times = di['times']
iters = di['iters']

# get the x axis
time_unit, time_label, rotation, simple_label = get_time_unit(dirname)
if not xiter:
    xaxis = times/time_unit
else:
    xaxis = iters

# set xminmax if not set by user
if xminmax is None:
    # set xmin possibly
    if xmin is None:
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

# deal with x axis, maybe thinning data
print ("ntot = %i" %ntot)
print ("before thin_data: len(xaxis) = %i" %len(xaxis))
xaxis = thin_data(xaxis, ntot)
times = thin_data(times, ntot)
iters = thin_data(iters, ntot)
vals_gav = thin_data(vals_gav, ntot)
if sep_czrz:
    vals_cz = thin_data(vals_cz, ntot)
    vals_rz = thin_data(vals_rz, ntot)
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
# Make all axes use scientific notation (except for y if logscale=True)
if logscale:
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

for irow in range(3): # tot, fluc, mean
    for icol in range(ncol):
        vals = vals_list[icol]
        # KINETIC ENERGY
        if irow == 0: # tot
            rke = vals[:, lut[402]]
            tke = vals[:, lut[403]]
            pke = vals[:, lut[404]]
            if plot_inte:
                inte = vals[:, lut[701]]
        if irow == 1: # fluc
            rke = vals[:, lut[410]]
            tke = vals[:, lut[411]]
            pke = vals[:, lut[412]]
        if irow == 2: # mean
            rke = vals[:, lut[402]] - vals[:, lut[410]]
            tke = vals[:, lut[403]] - vals[:, lut[411]]
            pke = vals[:, lut[404]] - vals[:, lut[412]]
        ke = rke + tke + pke

        # INTERNAL ENERGY
        if plot_inte:
            if irow == 1: # fluc, no inte for fluctuating S'
                inte = np.zeros(len(xaxis))
            else:
                inte = vals[:, lut[701]]
        
        # MAGNETIC ENERGY
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

        # make line plots
        # KINETIC
        # collect all the total energies together for min/max vals
        all_e = [rke, tke, pke, ke]

        ax = axs[irow, icol]
        ax.plot(xaxis, ke, color_order[0],\
                linewidth=lw_ke, label=r'$\rm{KE_{tot}}$')
        ax.plot(xaxis, rke, color_order[1],\
                linewidth=lw_ke, label=r'$\rm{KE_r}$')
        ax.plot(xaxis, tke, color_order[2],\
                linewidth=lw_ke, label=r'$\rm{KE_\theta}$')
        ax.plot(xaxis, pke, color_order[3],\
                linewidth=lw_ke, label=r'$\rm{KE_\phi}$')
        # INTERNAL
        if plot_inte:
            all_e += [inte]
            ax.plot(xaxis, inte, color_order[4], linewidth=lw_ke,\
                    label='INTE')

        # MAGNETIC
        if magnetism:
            if nodyn:
                tcut = tmin + dynfrac*(tmax - tmin)
                itcut = np.argmin(np.abs(times - tcut))
            else:
                itcut = 0
            all_e += [rme[itcut:], tme[itcut:], pme[itcut:], me[itcut:]]

            ax.plot(xaxis, me, color_order[0] + '--',\
                    linewidth=lw, label=r'$\rm{ME_{tot}}$')
            ax.plot(xaxis, rme, color_order[1] + '--',\
                    linewidth=lw, label=r'$\rm{ME_r}$')
            ax.plot(xaxis, tme, color_order[2] + '--',\
                    linewidth=lw, label=r'$\rm{ME_\theta}$')
            ax.plot(xaxis, pme, color_order[3] + '--',\
                    linewidth=lw, label=r'$\rm{ME_\phi}$')

        if irow == 0 and icol == 0: # put a legend on the upper left axis
            legfrac = 1/4
            ax.legend(loc='lower left', ncol=4, fontsize=0.8*fontsize, columnspacing=1)
        else:
            legfrac = None

        # set the y limits
        minmax_loc = minmax
        if not coords is None:
            if not (irow, icol) in coords: # reset minmax_loc to None
                # (will become default) if not in desired coordinates
                minmax_loc = None
        if minmax_loc is None:
            minmax_loc = lineplot_minmax(all_e, logscale=logscale, legfrac=legfrac)
        if not ymin is None:
            minmax_loc = ymin, minmax_loc[1]
        if not ymax is None:
            minmax_loc = minmax_loc[0], ymax
        ax.set_ylim((minmax_loc[0], minmax_loc[1]))

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
    axs[2,icol_mid].set_xlabel('time [' + time_label + ']')

# y labels
axs[0,0].set_ylabel('full energy', fontsize=fontsize)
axs[1,0].set_ylabel('fluc energy', fontsize=fontsize)
axs[2,0].set_ylabel('mean energy', fontsize=fontsize)

# Make titles
titles = ["ALL ZONES", "RZ", "CZ"]
for icol in range(ncol):
    if icol == icol_mid:
        title = dirname_stripped + '\n' + titles[icol]
    else:
        title = titles[icol]
    axs[0, icol].set_title(title)

# mark times if desired
for ax in axs.flatten():
    y1, y2 = ax.get_ylim()
    yvals = np.linspace(y1, y2, 100)
    for time in xvals:
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
tag = clas0['tag']
if xiter and tag == '':
    tag = '_xiter'
plotdir = my_mkdir(clas0['plotdir']) 
savename = 'etrace' + tag + '-' + str(iter1).zfill(8) + '_' + str(iter2).zfill(8) + '.png'

if clas0['saveplot']:
    print ('Saving the etrace plot at ' + plotdir + savename)
    plt.savefig(plotdir + savename, dpi=300)

# Show the plot
if clas0['showplot']:
    plt.show()
plt.close()
