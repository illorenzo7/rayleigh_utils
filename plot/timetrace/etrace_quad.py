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

# See if magnetism is "on"
magnetism = get_parameter(dirname, 'magnetism')

# SPECIFIC ARGS for etrace:
kwargs_default = dict({'the_file': None, 'xminmax': None, 'xmin': None, 'xmax': None, 'minmax': None, 'min': None, 'max': None, 'coords': None, 'ntot': 500, 'xiter': False, 'log': False, 'growth': False, 'growthfrac': 0.5, 'xvals': np.array([]), 'inte': False, 'nquadr': None, 'nquadlat': None, 'etype': 'tot', 'legfrac': None, 'nomag': False, 'noke': False})

# make figure kwargs
nlines = get_num_lines(clas0.dirname_label)
lineplot_fig_dimensions['margin_top_inches'] = (nlines+2)*default_line_height
make_figure_kwargs_default.update(lineplot_fig_dimensions)
make_figure_kwargs_default['margin_top_inches'] += 2*default_line_height

kwargs_default.update(make_figure_kwargs_default)

# plots two more columns with energies in CZ and RZ separately 
# update these defaults from command-line
kw = update_dict(kwargs_default, clas)
kw_make_figure = update_dict(make_figure_kwargs_default, clas)


fontsize = default_titlesize
the_file = kw.the_file
xminmax = kw.xminmax
xmin = kw.xmin
xmax = kw.xmax
minmax = kw.minmax
ymin = kw.min
ymax = kw.max
coords = kw.coords
ntot = kw.ntot
xiter = kw.xiter
logscale = kw.log
growth = kw.growth
growthfrac = kw.growthfrac
xvals = make_array(kw.xvals)
plot_inte = kw.inte
nquadlat = kw.nquadlat
nquadr = kw.nquadr
etype = kw.etype
legfrac = kw.legfrac
noke = kw.noke
nomag = kw.nomag

# deal with coords (if user wants minmax to only apply to certain subplots)
if not coords is None:
    numpanels = len(coords)//2
    acopy = np.copy(coords)
    coords = []
    for i in range(numpanels):
        coords.append((acopy[2*i], acopy[2*i + 1]))

# get desired data file
dataname = 'G_Avgs_trace'
if the_file is None:
    if not nquadlat is None:
        dataname += '_nquadlat%i' %nquadlat
    if not nquadr is None:
        dataname += '_nquadr%i' %nquadr
    the_file = get_widest_range_file(clas0['datadir'], dataname)

print ('Getting data from ' + the_file)
di = get_dict(the_file)
vals = di['vals']
rvals = di['rvals']
latvals = di['latvals']
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
ntimes, nq, nquadlat, nquadr = np.shape(vals)
nplots = nquadlat*nquadr

# now finally get the shape of the "vals" array
ntimes, nq, nquadlat, nquadr = np.shape(vals)
ntimes = len(xaxis)
nplots = nquadlat*nquadr

# create the figure dimensions
kw_make_figure.nplots = nplots
kw_make_figure.ncol = nquadr
fig, axs, fpar = make_figure(**kw_make_figure)

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
        ax = axs[ilat, ir]

        all_e = []
        # choose etype: tot, fluc, or mean of energies

        # KINETIC ENERGY
        if etype == 'tot':
            rke = vals_loc[:, lut[402]]
            tke = vals_loc[:, lut[403]]
            pke = vals_loc[:, lut[404]]
            label_pre = ''
            label_app = ''
        if etype == 'fluc':
            rke = vals_loc[:, lut[410]]
            tke = vals_loc[:, lut[411]]
            pke = vals_loc[:, lut[412]]
            label_pre = ''
            label_app = '\''
        if etype == 'mean':
            rke = vals_loc[:, lut[402]] - vals_loc[:, lut[410]]
            tke = vals_loc[:, lut[403]] - vals_loc[:, lut[411]]
            pke = vals_loc[:, lut[404]] - vals_loc[:, lut[412]]
            label_pre = '<'
            label_app = '>'
        ke = rke + tke + pke

        # INTERNAL ENERGY
        if plot_inte:
            if etype == 'fluc': # no inte for fluctuating S'
                inte = np.zeros(len(xaxis))
            else:
                inte = vals_loc[:, lut[701]]
        
        # MAGNETIC ENERGY
        if magnetism:
            if etype == 'tot':
                rme = vals_loc[:, lut[1102]]
                tme = vals_loc[:, lut[1103]]
                pme = vals_loc[:, lut[1104]]
            if etype == 'fluc':
                rme = vals_loc[:, lut[1110]]
                tme = vals_loc[:, lut[1111]]
                pme = vals_loc[:, lut[1112]]
            if etype == 'mean':
                rme = vals_loc[:, lut[1102]] - vals_loc[:, lut[1110]]
                tme = vals_loc[:, lut[1103]] - vals_loc[:, lut[1111]]
                pme = vals_loc[:, lut[1104]] - vals_loc[:, lut[1112]]
            me = rme + tme + pme

        # make line plots

        # KINETIC
        # collect all the total energies together for min/max vals

        # see if we should plot the growth phase
        if logscale: # only need to worry about this for log-scale
            if growth:
                itcut = 0
            else:
                tcut = tmin + growthfrac*(tmax - tmin)
                itcut = np.argmin(np.abs(times - tcut))
        else:
            itcut = 0

        if not noke:
            all_e += [rke[itcut:], tke[itcut:], pke[itcut:], ke[itcut:]]
            
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
            all_e += [inte[itcut:]]
            ax.plot(xaxis, inte, color_order[4], linewidth=lw_ke,\
                    label='INTE')

        # MAGNETIC
        if not nomag:
            if magnetism:
                all_e += [rme[itcut:], tme[itcut:], pme[itcut:], me[itcut:]]

                ax.plot(xaxis, me, color=color_order[0], linestyle='--',\
                        linewidth=lw, label=r'$\rm{ME_{tot}}$')
                ax.plot(xaxis, rme, color=color_order[1], linestyle='--',\
                        linewidth=lw, label=r'$\rm{ME_r}$')
                ax.plot(xaxis, tme, color=color_order[2], linestyle='--',\
                        linewidth=lw, label=r'$\rm{ME_\theta}$')
                ax.plot(xaxis, pme, color=color_order[3], linestyle='--',\
                        linewidth=lw, label=r'$\rm{ME_\phi}$')

        if ilat == 0 and ir == 0: # put a legend on the upper left axis
            plotleg = True
            ax.legend(loc='lower left', ncol=3, fontsize=0.7*fontsize, columnspacing=1)
        else:
            plotleg = False

        # set the y limits
        minmax_loc = minmax
        if not coords is None:
            if not (ilat, ir) in coords: # reset minmax_loc to None
                # (will become default) if not in desired coordinates
                minmax_loc = None
        if minmax_loc is None:
            minmax_loc = lineplot_minmax(xaxis, all_e, logscale=logscale, legfrac=legfrac, plotleg=plotleg)
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
    r1 = rvals[ir]
    r2 = rvals[ir+1]
    title = 'rad. range =\n [%1.2e, %1.2e]' %(r1, r2)
    axs[0, ir].set_title(title, fontsize=fontsize)

# y labels / side titles
for it in range(nquadlat):
    lat1 = latvals[it]
    lat2 = latvals[it+1]
    axs[it, 0].set_ylabel('lat. range = [%.1f, %.1f]' %(lat1, lat2), fontsize=fontsize)

# overall title 
iter1, iter2 = get_iters_from_file(the_file)
time_string = get_time_string(dirname, iter1, iter2) 
the_title = clas0.dirname_label + '\n' +  'energy trace (' + etype + ')' + '\n' + time_string
margin_x = fpar['margin_left'] + fpar['sub_margin_left']
margin_y = default_margin/fpar['height_inches']
fig.text(margin_x, 1 - margin_y, the_title,\
         ha='left', va='top', fontsize=default_titlesize)


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
tag = clas0['tag']
if xiter and tag == '':
    tag = '_xiter'
plotdir = my_mkdir(clas0['plotdir']) 
basename = dataname.replace('G_Avgs_trace', 'etrace' + '_' + etype)
savename = basename + tag + '-' + str(iter1).zfill(8) + '_' + str(iter2).zfill(8) + '.png'
if clas0.prepend:
    savename = dirname_stripped + '_' + savename

if clas0['saveplot']:
    print ('Saving the etrace plot at ' + plotdir + savename)
    plt.savefig(plotdir + savename, dpi=300)

# Show the plot
if clas0['showplot']:
    plt.show()
plt.close()
