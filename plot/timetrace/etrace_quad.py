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
kw_default = dict({'the_file': None, 'minmax': None, 'xmin': None, 'xmax': None, 'minmax': None, 'min': None, 'max': None, 'coords': None, 'ntot': 500, 'xiter': False, 'log': False, 'growth': False, 'growthfrac': 0.5, 'xvals': np.array([]), 'inte': False, 'nquadr': 1, 'nquadlat': 1, 'type': 'tot', 'legfrac': None, 'nomag': False, 'noke': False, 'tkappa': False, 'legloc': 'lower left'})

# make figure kw
nlines = get_num_lines(clas0.dirname_label)
lineplot_fig_dimensions['margin_top_inches'] = (nlines+2)*default_line_height
kw_make_figure_default.update(lineplot_fig_dimensions)
kw_make_figure_default['margin_top_inches'] += 2*default_line_height

kw_default.update(kw_make_figure_default)

# plots two more columns with energies in CZ and RZ separately 
# update these defaults from command-line
kw = update_dict(kw_default, clas)
kw_make_figure = update_dict(kw_make_figure_default, clas)
fontsize = default_titlesize

# deal with kw.coords (if user wants kw.minmax to only apply to certain subplots)
if not kw.coords is None:
    numpanels = len(kw.coords)//2
    acopy = np.copy(kw.coords)
    kw.coords = []
    for i in range(numpanels):
        kw.coords.append((acopy[2*i], acopy[2*i + 1]))

# get desired data file
dataname = 'G_Avgs_trace'
if kw.the_file is None:
    # need to get widest range file
    if kw.nquadlat > 1:
        dataname += '_nquadlat%i' %kw.nquadlat
    if kw.nquadr > 1:
        dataname += '_nquadr%i' %kw.nquadr
    kw.the_file = get_widest_range_file(clas0['datadir'], dataname)

print ('Getting data from ' + kw.the_file)
di = get_dict(kw.the_file)
vals = di['vals']
rvals = di['rvals']
latvals = di['latvals']
lut = di['lut']
qv = di['qv']
times = di['times']
iters = di['iters']

# get the x axis
time_unit, time_label, rotation, simple_label = get_time_unit(dirname, tkappa=kw.tkappa)
if not kw.xiter:
    xaxis = times/time_unit
else:
    xaxis = iters

# set kw.xminmax if not set by user
if kw.xminmax is None:
    # set kw.xmin possibly
    if kw.xmin is None:
        kw.xmin = xaxis[0]
    # set kw.xmax possibly
    if kw.xmax is None:
        kw.xmax = xaxis[-1]
    kw.xminmax = kw.xmin, kw.xmax

ixmin = np.argmin(np.abs(xaxis - kw.xminmax[0]))
ixmax = np.argmin(np.abs(xaxis - kw.xminmax[1]))

# Now shorten all the "x" arrays
xaxis = xaxis[ixmin:ixmax+1]
times = times[ixmin:ixmax+1]
iters = iters[ixmin:ixmax+1]
tmin, tmax = times[0], times[-1]
vals = vals[ixmin:ixmax+1, :]

# deal with x axis, maybe thinning data
if np.all(kw.ntot == 'full'):
    print ('kw.ntot = full')
    kw.ntot = len(times)
print ("kw.ntot = %i" %kw.ntot)
print ("before thin_data: len(xaxis) = %i" %len(xaxis))
xaxis = thin_data(xaxis, kw.ntot)
times = thin_data(times, kw.ntot)
iters = thin_data(iters, kw.ntot)
vals = thin_data(vals, kw.ntot)
print ("after thin_data: len(xaxis) = %i" %len(xaxis))

# now finally get the shape of the "vals" array
ntimes, nq, nquadlat, nquadr = np.shape(vals)
nplots = nquadlat*nquadr

# create the figure dimensions
kw_make_figure.nplots = nplots
kw_make_figure.ncol = nquadr
fig, axs, fpar = make_figure(**kw_make_figure)

# Make thin lines to see structure of variation for ME
lw = 0.5
lw_ke = 1. # bit thicker for KE to differentiate between ME

# See if y-axis should be on log scale (must do this before setting limits)
# Make all axes use scientific notation (except for y if kw.log=True)
if kw.log:
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

        all_e = [] # always need this for max val
        # choose kw.type: tot, fluc, or mean of energies

        # KINETIC ENERGY
        if kw.type == 'tot':
            rke = vals_loc[:, lut[402]]
            tke = vals_loc[:, lut[403]]
            pke = vals_loc[:, lut[404]]
            label_pre = ''
            label_app = ''
        if kw.type == 'fluc':
            rke = vals_loc[:, lut[410]]
            tke = vals_loc[:, lut[411]]
            pke = vals_loc[:, lut[412]]
            label_pre = ''
            label_app = '\''
        if kw.type == 'mean':
            rke = vals_loc[:, lut[402]] - vals_loc[:, lut[410]]
            tke = vals_loc[:, lut[403]] - vals_loc[:, lut[411]]
            pke = vals_loc[:, lut[404]] - vals_loc[:, lut[412]]
            label_pre = '<'
            label_app = '>'
        ke = rke + tke + pke

        # INTERNAL ENERGY
        if kw.inte:
            if kw.type == 'fluc': # no inte for fluctuating S'
                inte = np.zeros(len(xaxis))
            else:
                inte = vals_loc[:, lut[701]]
        
        # MAGNETIC ENERGY
        # get what we can from component energies,
        # otherwise try to get total energies directly (for Lydia's stuff)
        if magnetism:
            if kw.type == 'tot':
                if issubset([1102,1103,1104], qv): # get components
                    have_comp = True
                    rme = vals_loc[:, lut[1102]]
                    tme = vals_loc[:, lut[1103]]
                    pme = vals_loc[:, lut[1104]]
                else:
                    have_comp = False
                    me = vals_loc[:, lut[1101]]
            if kw.type == 'fluc':
                if issubset([1110,1111,1112], qv): # get components
                    have_comp = True
                    rme = vals_loc[:, lut[1110]]
                    tme = vals_loc[:, lut[1111]]
                    pme = vals_loc[:, lut[1112]]
                else:
                    have_comp = False
                    me = vals_loc[:, lut[1109]]
            if kw.type == 'mean':
                if issubset([1102,1103,1104,1110,1111,1112], qv): # get components
                    have_comp = True
                    rme = vals_loc[:, lut[1102]] - vals_loc[:, lut[1110]]
                    tme = vals_loc[:, lut[1103]] - vals_loc[:, lut[1111]]
                    pme = vals_loc[:, lut[1104]] - vals_loc[:, lut[1112]]
                else:
                    have_comp = False
                    me = vals_loc[:, lut[1105]]
            if have_comp:
                me = rme + tme + pme

        # make line plots

        # KINETIC
        # collect all the total energies together for min/max vals

        # see if we should plot the kw.growth phase
        if kw.log: # only need to worry about this for log-scale
            if kw.growth:
                itcut = 0
            else:
                tcut = tmin + kw.growthfrac*(tmax - tmin)
                itcut = np.argmin(np.abs(times - tcut))
        else:
            itcut = 0

        if not kw.noke:
            all_e += [rke, tke, pke, ke]
            
            ax.plot(xaxis, ke, color_order[0],\
                    linewidth=lw_ke, label=r'$\rm{KE_{tot}}$')
            ax.plot(xaxis, rke, color_order[1],\
                    linewidth=lw_ke, label=r'$\rm{KE_r}$')
            ax.plot(xaxis, tke, color_order[2],\
                    linewidth=lw_ke, label=r'$\rm{KE_\theta}$')
            ax.plot(xaxis, pke, color_order[3],\
                    linewidth=lw_ke, label=r'$\rm{KE_\phi}$')

        # INTERNAL
        if kw.inte:
            all_e += [inte]
            ax.plot(xaxis, inte, color_order[4], linewidth=lw_ke,\
                    label='INTE')

        # MAGNETIC
        if not kw.nomag:
            if magnetism:
                all_e += [me]
                if have_comp:
                    all_e += [rme, tme, pme]

                ax.plot(xaxis, me, color=color_order[0], linestyle='--',\
                        linewidth=lw, label=r'$\rm{ME_{tot}}$')
                if have_comp:
                    ax.plot(xaxis, rme, color=color_order[1], linestyle='--',\
                            linewidth=lw, label=r'$\rm{ME_r}$')
                    ax.plot(xaxis, tme, color=color_order[2], linestyle='--',\
                            linewidth=lw, label=r'$\rm{ME_\theta}$')
                    ax.plot(xaxis, pme, color=color_order[3], linestyle='--',\
                            linewidth=lw, label=r'$\rm{ME_\phi}$')

        if ilat == 0 and ir == 0: # put a legend on the upper left axis
            plotleg = True
            ax.legend(loc=kw.legloc, ncol=3, fontsize=0.7*fontsize, columnspacing=1)
        else:
            plotleg = False

        # set the y limits
        minmax_loc = kw.minmax
        if not kw.coords is None:
            if not (ilat, ir) in kw.coords: # reset minmax_loc to None
                # (will become default) if not in desired coordinates
                minmax_loc = None
        if minmax_loc is None:
            minmax_loc = lineplot_minmax(xaxis, all_e, log=kw.log, legfrac=kw.legfrac, plotleg=plotleg, ixcut=itcut)
        if not kw.min is None:
            minmax_loc = kw.min, minmax_loc[1]
        if not kw.max is None:
            minmax_loc = minmax_loc[0], kw.max
        ax.set_ylim((minmax_loc[0], minmax_loc[1]))

        ax.set_xlim((kw.xminmax[0], kw.xminmax[1]))
if kw.xiter:
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
iter1, iter2 = get_iters_from_file(kw.the_file)
time_string = get_time_string(dirname, iter1, iter2) 
the_title = clas0.dirname_label + '\n' +  'energy trace (' + kw.type + ')' + '\n' + time_string
margin_x = fpar['margin_left'] + fpar['sub_margin_left']
margin_y = default_margin/fpar['height_inches']
fig.text(margin_x, 1 - margin_y, the_title,\
         ha='left', va='top', fontsize=default_titlesize)

# mark times if desired
for ax in axs.flatten():
    y1, y2 = ax.get_ylim()
    yvals = np.linspace(y1, y2, 100)
    for time in kw.xvals:
        ax.plot(time + np.zeros(100), yvals,'k--')

# Get ticks everywhere
for ax in axs.flatten():
    plt.sca(ax)
    plt.minorticks_on()
    plt.tick_params(top=True, right=True, direction='in', which='both')

# Save the plot
tag = clas0['tag']
if kw.xiter and tag == '':
    tag = '_kw.xiter'
plotdir = my_mkdir(clas0['plotdir']) 
basename = dataname.replace('G_Avgs_trace', 'etrace' + '_' + kw.type)
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
