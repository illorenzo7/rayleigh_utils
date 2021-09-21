# Author: Loren Matilsky
# Created: 08/17/2021
# plot the trace in different quadrants
# This script plots the quantitities specified by --qvals or --groupname
# default is --groupname v

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
kw_default = dict({'the_file': None, 'xminmax': None, 'xmin': None, 'xmax': None, 'minmax': None, 'min': None, 'max': None, 'coords': None, 'ntot': 500, 'xiter': False, 'log': False, 'xvals': np.array([]), 'nquadr': None, 'nquadlat': None, 'qvals': None, 'groupname': 'b', 'totsig': None, 'titles': None, 'justtot': False})

# make figure kwargs
lineplot_fig_dimensions['margin_top_inches'] = 3/4
make_figure_kwargs_default.update(lineplot_fig_dimensions)
kw_default.update(make_figure_kwargs_default)

# lineplot kwargs
kw_default.update(lineplot_kwargs_default)

# then override defaults
kw = update_dict(kw_default, clas)
kw_make_figure = update_dict(make_figure_kwargs_default, clas)
kw_lineplot = update_dict(lineplot_kwargs_default, clas)

# deal with desired quantities
if kw.qvals is None: # it's a quantity group
    qgroup = get_quantity_group(kw.groupname, magnetism)
    kw.qvals = qgroup['qvals']
    if kw.titles is None:
        kw.titles = qgroup['titles']
    if kw.totsig is None:
        kw.totsig = qgroup['totsig']
else:
    kw.titles = parse_quantities(kw.qvals)[1]
    kw.groupname = input("choose a groupname to save your plot\n to not save it, enter 'nosave': ")

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
xvals = make_array(kw.xvals)
nquadlat = kw.nquadlat
nquadr = kw.nquadr

# state what we're plotting
print ("plotting the following quantities:")
print ("qvals = ", kw.qvals)

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
rbounds = di['rbounds']
nquadr = len(rbounds) - 1
latbounds = di['latbounds']
nquadlat = len(latbounds) - 1
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
#times = times[ixmin:ixmax+1]
iters = iters[ixmin:ixmax+1]
#tmin, tmax = times[ixmin], times[ixmax]
#vals = vals[ixmin:ixmax+1, :] # keep the full vals array for now

# deal with x axis, maybe thinning data
print ("ntot = %i" %ntot)
print ("before thin_data: len(xaxis) = %i" %len(xaxis))
xaxis = thin_data(xaxis, ntot)
#times = thin_data(times, ntot)
iters = thin_data(iters, ntot)
#vals = thin_data(vals, ntot) # don't thin the data quite yet for this guy
print ("after thin_data: len(xaxis) = %i" %len(xaxis))

# now finally get the shape of the "vals" array
ntimes, nq, nquadlat, nquadr = np.shape(vals)
ntimes = len(xaxis)
nplots = nquadlat*nquadr

# create the figure dimensions
kw_make_figure.nplots = nplots
kw_make_figure.ncol = nquadr
fig, axs, fpar = make_figure(**kw_make_figure)

# loop over different domains
for ilat in range(nquadlat):
    for ir in range(nquadr):
        vals_loc = vals[:, :, ilat, ir]
        ax = axs[ilat, ir]

        # get terms we want and plot them
        terms = []
        kw_lineplot.labels = []
        nterms = len(kw.qvals)
        for iterm in range(nterms):
            qval = kw.qvals[iterm]
            terms.append(vals_loc[:, lut[int(qval)]])
            kw_lineplot.labels.append(kw.titles[iterm])

        # might also need the total of these terms (with some signature: totsig)

        if not kw.totsig is None:
            tot_term = np.zeros_like(terms[0])
            for iterm in range(nterms):
                tot_term += terms[iterm]*kw.totsig[iterm]
            terms.append(tot_term)
            nterms += 1
            kw_lineplot.labels.append('sum')

        equation_sets = ['torque', 'teq', 'forcer', 'forcet', 'forcep', 'indr', 'indt', 'indp']
        if kw.groupname in equation_sets:
            # replace term with its time derivative
            terms[0] = drad(terms[0], times)
            kw_lineplot.labels[0] = 'd/dt'
            
            kw_lineplot.legfrac = 0.6

        if kw.justtot:
            diff = terms[-1] - terms[0]
            terms = [terms[0], terms[-1], diff]
            nterms = 3
            kw_lineplot.labels = ['d/dt (LHS)', 'sum (RHS)', 'difference']
            kw_lineplot.linestyles = ['-', '--', ':']
            kw_lineplot.colors = ['k', 'r', 'g']
            
        # now thin the data on the terms #and times
        #times = thin_data(times[ixmin:ixmax+1], ntot)
        for iterm in range(nterms):
            terms[iterm] = thin_data(terms[iterm][ixmin:ixmax+1], ntot)

        # now plot the terms
        if ilat == 0 and ir == 0:
            kw_lineplot.plotleg = True
        lineplot(xaxis, terms, ax, **kw_lineplot)

        if kw.groupname in equation_sets:
            # label the derivative ampltitude and sum amplitude
            if kw.justtot:
                tot = terms[-2]
            else:
                tot = terms[-1]
                diff = tot - terms[0]
            quant = 'rms(d/dt) = %1.3e' %rms(terms[0])
            quant += '\n' + 'rms(sum) = %1.3e' %rms(tot)
            quant += '\n' + 'rms(d/dt - sum) = %1.3e' %rms(diff)
            quant += '\n' + 'err = %1.3e' %(rms(diff)/rms(terms[0]))
            xmin, xmax = ax.get_xlim()
            dx = xmax - xmin
            ymin, ymax = ax.get_ylim()
            dy = ymax - ymin
            ax.text(xmin + 0.1*dx, ymin + 0.3*dy, quant, va='bottom', ha='left')
        

# Set some parameters defining all subplots
# x limits and label
axs[0, 0].set_xlim((xminmax[0], xminmax[1]))
if xiter:
    axs[-1, 0].set_xlabel('iteration #')
else:
    axs[-1, 0].set_xlabel('time [' + time_label + ']')

# x titles
for ir in range(nquadr):
    r1 = rbounds[ir]
    r2 = rbounds[ir+1]
    title = 'rad. range = [%.3f, %.3f]' %(r1, r2)
    axs[0, ir].set_title(title, fontsize=fontsize)

# y labels
for it in range(nquadlat):
    lat1 = latbounds[it]
    lat2 = latbounds[it+1]
    axs[it, 0].set_ylabel('lat. range = [%.1f, %.1f]' %(lat1, lat2), fontsize=fontsize)

# main title
maintitle = dirname_stripped
if not kw.groupname is None:
    maintitle += '\n' + 'groupname = %s' %kw.groupname
maintitle += '\n' + 'qvals = ' +(arr_to_str(kw.qvals, '%i'))
# Put the main title in upper left
fig.text(fpar['margin_left'] + fpar['sub_margin_left'], 1.0 - default_margin/fpar['height_inches'], maintitle, ha='left', va='top', fontsize=default_titlesize)

# save the figure if tag (or qgroup) was specified
if len(clas0['tag']) > 0 or not kw.groupname is None:
    basename = dataname.replace('G_Avgs_trace', 'timetrace')

    if not kw.groupname is None:
        basename += '_' + kw.groupname
    basename += clas0['tag']
    iter1, iter2 = get_iters_from_file(the_file)
    savename = basename + '-' + str(iter1).zfill(8) + '_' + str(iter2).zfill(8) + '.png'

    plotdir = my_mkdir(clas0['plotdir'] + 'timetrace/')

    iter1, iter2 = get_iters_from_file(the_file)
    savefile = plotdir + basename + '-' + str(iter1).zfill(8) + '_' + str(iter2).zfill(8) + '.png'
    print ('saving figure at ' + savefile)
    fig.savefig(savefile, dpi=300)

# Show the plot
if clas0['showplot']:
    plt.show()
plt.close()
