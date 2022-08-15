# Author: Loren Matilsky
# plot the energy in different quadrants
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

# SPECIFIC ARGS for etrace:
kwargs_default = dict({'the_file': None, 'xminmax': None, 'xmin': None, 'xmax': None, 'minmax': None, 'min': None, 'max': None, 'coords': None, 'ntot': 500, 'xiter': False, 'xvals': np.array([]), 'nquadr': 1, 'legfrac': None})

# make figure kwargs
lineplot_fig_dimensions['margin_top_inches'] = 3/4
make_figure_kwargs_default.update(lineplot_fig_dimensions)
kwargs_default.update(make_figure_kwargs_default)

# plots two more columns with energies in CZ and RZ separately 
# update these defaults from command-line
kwargs = update_dict(kwargs_default, clas)
kw_make_figure = update_dict(make_figure_kwargs_default, clas)

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
xvals = make_array(kwargs.xvals)
nquadr = kwargs.nquadr
legfrac = kwargs.legfrac

# deal with coords (if user wants minmax to only apply to certain subplots)
if not coords is None:
    numpanels = len(coords)//2
    acopy = np.copy(coords)
    coords = []
    for i in range(numpanels):
        coords.append((acopy[2*i], acopy[2*i + 1]))

# get desired data file
dataname = 'poynt_trace'
if the_file is None:
    dataname += '_nquadr%i' %nquadr
    the_file = get_widest_range_file(clas0['datadir'], dataname)

print ('Getting data from ' + the_file)
di = get_dict(the_file)
vals = di['vals']
rbounds = di['rbounds']
volumes = di['volumes']
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
vals = vals[ixmin:ixmax+1, :]

# deal with x axis, maybe thinning data
if np.all(ntot == 'full'):
    print ('ntot = full')
    ntot = len(times)
print ("ntot = %i" %ntot)
print ("before thin_data: len(xaxis) = %i" %len(xaxis))
xaxis = thin_data(xaxis, ntot)
times = thin_data(times, ntot)
iters = thin_data(iters, ntot)
vals = thin_data(vals, ntot)
print ("after thin_data: len(xaxis) = %i" %len(xaxis))

# now finally get the shape of the "vals" array
ntimes, nq, nquadr = np.shape(vals)
# nq should = 7
nplots = nquadr

# create the figure dimensions
kw_make_figure.nplots = nplots
kw_make_figure.ncol = nquadr
fig, axs, fpar = make_figure(**kw_make_figure)

# Make thin lines to see structure of variation
lw = 0.5

# Make all axes use scientific notation 
for ax in axs.flatten():
    ax.ticklabel_format(scilimits = (-3,4), useMathText=True)

# loop over different domains
for ir in range(nquadr):
    vals_loc = vals[:, :, ir]
    ax = axs[0, ir]

    # Get 5 quantities in each quadrant:
    # total ME, generation by flows, dissipation by eta,
    # and two Poynting flluxes
    me = vals_loc[:, 0]
    dmedt = np.gradient(me, times)
    v_work = vals_loc[:, 4]
    diss = vals_loc[:, 5]
    poynt_bot = -vals[:, 6, ir] # remember bottom one needs negative
    poynt_top = vals[:, 6, ir+1]
    # also need to multiply poynting flux by surface area, then 
    # divide by the volume of the shell to get energy change per unit vol
    poynt_bot *= (4*np.pi*rbounds[ir]**2)
    poynt_top *= (4*np.pi*rbounds[ir+1]**2)
    poynt_bot /= volumes[ir]
    poynt_top /= volumes[ir]
    the_sum = v_work + diss + poynt_bot + poynt_top

    # make line plots
    # need to collect everything for profiles
    all_terms = [dmedt, diss, v_work, poynt_bot, poynt_top, the_sum]

    ax.plot(xaxis, dmedt, color_order[0],\
            linewidth=lw, label='dME/dt')
    ax.plot(xaxis, v_work, color_order[1],\
            linewidth=lw, label='v work')
    ax.plot(xaxis, diss, color_order[2],\
            linewidth=lw, label='ohm loss')
    ax.plot(xaxis, poynt_bot, color_order[3],\
            linewidth=lw, label='flux bot')
    ax.plot(xaxis, poynt_top, color_order[4],\
            linewidth=lw, label='flux top')
    ax.plot(xaxis, the_sum, color_order[5],\
            linewidth=lw, label='sum RHS')

    if ir == 0: # put a legend on the upper left axis
        plotleg = True
        ax.legend(loc='lower left', ncol=3, fontsize=0.7*fontsize, columnspacing=1)
    else:
        plotleg = False

    # set the y limits
    minmax_loc = minmax
    if not coords is None:
        if not (0, ir) in coords: # reset minmax_loc to None
            # (will become default) if not in desired coordinates
            minmax_loc = None
    if minmax_loc is None:
        minmax_loc = lineplot_minmax(xaxis, all_terms, legfrac=legfrac, plotleg=plotleg)
    if not ymin is None:
        minmax_loc = ymin, minmax_loc[1]
    if not ymax is None:
        minmax_loc = minmax_loc[0], ymax
    ax.set_ylim((minmax_loc[0], minmax_loc[1]))

    ax.set_xlim((xminmax[0], xminmax[1]))
if xiter:
    axs[-1, 0].set_xlabel('iteration #')
else:
    axs[-1, 0].set_xlabel('time [' + time_label + ']')

# x titles
for ir in range(nquadr):
    r1 = rbounds[ir]
    r2 = rbounds[ir+1]
    title = 'rad. range = [%.3f, %.3f]' %(r1, r2)
    if ir == 0:
        title = dirname_stripped + '\n' + title
    axs[0, ir].set_title(title, fontsize=fontsize)

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

# Save the plot
iter1, iter2 = get_iters_from_file(the_file)
# Tag the plot by whether or not the x axis is in "time" or "iteration"
tag = clas0['tag']
if xiter and tag == '':
    tag = '_xiter'
plotdir = my_mkdir(clas0['plotdir']) 
basename = dataname
savename = basename + tag + '-' + str(iter1).zfill(8) + '_' + str(iter2).zfill(8) + '.png'

if clas0['saveplot']:
    print ('Saving the etrace plot at ' + plotdir + savename)
    plt.savefig(plotdir + savename, dpi=300)

# Show the plot
if clas0['showplot']:
    plt.show()
plt.close()
