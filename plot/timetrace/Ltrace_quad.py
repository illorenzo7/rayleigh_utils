# Author: Loren Matilsky
# plot the amom. in different quadrants
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

# SPECIFIC ARGS for etrace:
kwargs_default = dict({'the_file': None, 'xminmax': None, 'xmin': None, 'xmax': None, 'minmax': None, 'min': None, 'max': None, 'coords': None, 'ntot': 500, 'xiter': False, 'log': False, 'xvals': np.array([]), 'plotall': False, 'vol': False})
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
xvals = make_array(kwargs.xvals)
plotall = kwargs.plotall
vol = kwargs.vol

# deal with coords (if user wants minmax to only apply to certain subplots)
if not coords is None:
    numpanels = len(coords)//2
    acopy = np.copy(coords)
    coords = []
    for i in range(numpanels):
        coords.append((acopy[2*i], acopy[2*i + 1]))

# Might need to use 2dom trace instead of regular trace
if the_file is None:
    the_file = get_widest_range_file(clas0['datadir'], 'G_Avgs_trace_quad')

print ('Getting data from ' + the_file)
di = get_dict(the_file)
vals = di['vals']
rbounds = di['rbounds']
latbounds = di['latbounds']
lut = di['lut']
times = di['times']
iters = di['iters']
volumes = di['volumes']

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
print ("ntot = %i" %ntot)
print ("before thin_data: len(xaxis) = %i" %len(xaxis))
xaxis = thin_data(xaxis, ntot)
times = thin_data(times, ntot)
iters = thin_data(iters, ntot)
vals = thin_data(vals, ntot)
print ("after thin_data: len(xaxis) = %i" %len(xaxis))

# now finally get the shape of the "vals" array
ntimes, nq, nquadlat, nquadr = np.shape(vals)
nplots = nquadlat*nquadr

# create figure with nquadr columns and nquadlat rows
fig, axs = plt.subplots(nquadlat, nquadr, figsize=(3.5*nquadr, 10), sharex=True)
if nquadr == 1: # need the axis array to consistently be doubly indexed
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

# loop over different domains
for ilat in range(nquadlat):
    for ir in range(nquadr):
        vals_loc = vals[:, :, ilat, ir]
        if vol: # multiply by the volume of the quadrant
            print (" vol = ", volumes[ilat,ir])
            vals_loc *= volumes[ilat, ir]
        ax = axs[ilat, ir]

        all_L = []
        for i in range(3): # tot, fluc, mean of energies
            # z-mom:
            if i == 0: # tot
                Lz = vals_loc[:, lut[1819]]
                if plotall:
                    Lx = vals_loc[:, lut[1820]]
                    Ly = vals_loc[:, lut[1821]]
                linestyle = '-'
                label_pre = ''
                label_app = ''

            if i == 1: # fluc
                Lx = vals_loc[:, lut[1822]]
                if plotall:
                    Lx = vals_loc[:, lut[1823]]
                    Ly = vals_loc[:, lut[1824]]
                linestyle = ':'
                label_pre = ''
                label_app = '\''
            if i == 2: # mean
                Lz = vals_loc[:, lut[1819]] - vals_loc[:, lut[1822]]
                if plotall:
                    Lx = vals_loc[:, lut[1820]] - vals_loc[:, lut[1823]]
                    Ly = vals_loc[:, lut[1821]] - vals_loc[:, lut[1824]]
                linestyle = '--'
                label_pre = '<'
                label_app = '>'

            # make line plots

            ax.plot(xaxis, Lz, color_order[0], linewidth=lw, label=label_pre+r'$\rm{\mathcal{L}_{z}}$'+label_app)
            # collect all the momenta for min/max values
            all_L = [Lz]
            if plotall:
                ax.plot(xaxis, Lx, color_order[0], linewidth=lw, label=label_pre+r'$\rm{\mathcal{L}_{x}}$'+label_app)
                ax.plot(xaxis, Ly, color_order[0], linewidth=lw, label=label_pre+r'$\rm{\mathcal{L}_{y}}$'+label_app)
                all_L += [Lx, Ly]

        if ilat == 0 and ir == 0: # put a legend on the upper left axis
            #legfrac = 1/4
            legfrac = 1/2
            ax.legend(loc='lower left', ncol=3, fontsize=0.7*fontsize, columnspacing=1)
        else:
            legfrac = None

        # set the y limits
        minmax_loc = minmax
        if not coords is None:
            if not (it, ir) in coords: # reset minmax_loc to None
                # (will become default) if not in desired coordinates
                minmax_loc = None
        if minmax_loc is None:
            minmax_loc = lineplot_minmax(xaxis, all_L, logscale=logscale, legfrac=legfrac)
        if not ymin is None:
            minmax_loc = ymin, minmax_loc[1]
        if not ymax is None:
            minmax_loc = minmax_loc[0], ymax
        ax.set_ylim((minmax_loc[0], minmax_loc[1]))

# Set some parameters defining all subplots
# x limits and label
axs[0, 0].set_xlim((xminmax[0], xminmax[1]))
if xiter:
    axs[2, 0].set_xlabel('iteration #')
else:
    axs[2, 0].set_xlabel('time [' + time_label + ']')

# x titles
for ir in range(nquadr):
    r1 = rbounds[ir]
    r2 = rbounds[ir+1]
    title = 'rad. range = [%.3f, %.3f]' %(r1, r2) + '\n'
    if ir == 0:
        title = dirname_stripped + '\n' + title
    axs[0, ir].set_title(title, fontsize=fontsize)

# y labels
for it in range(nquadlat):
    lat1 = latbounds[it]
    lat2 = latbounds[it+1]
    axs[it, 0].set_ylabel('lat. range = [%.1f, %.1f]' %(lat1, lat2), fontsize=fontsize)

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
savename = 'Ltrace_quad' + tag + '-' + str(iter1).zfill(8) + '_' + str(iter2).zfill(8) + '.png'

if clas0['saveplot']:
    print ('Saving the etrace plot at ' + plotdir + savename)
    plt.savefig(plotdir + savename, dpi=300)

# Show the plot
if clas0['showplot']:
    plt.show()
plt.close()