# Author: Loren Matilsky
# Created: 08/17/2021
# plot the trace in different quadrants
# This script plots the quantitities specified by --qvals or --groupname
# default is --groupname v

import matplotlib.pyplot as plt
import numpy as np
from scipy import integrate 
import sys, os
sys.path.append(os.environ['raco'])
sys.path.append(os.environ['rapp'])
sys.path.append(os.environ['rapl'])
from common import *
from plotcommon import *
from cla_util import *
from rayleigh_diagnostics import GridInfo
cumtrap = integrate.cumulative_trapezoid

# Get the run directory on which to perform the analysis
args = sys.argv
clas0, clas = read_clas(args)
dirname = clas0['dirname']
dirname_stripped = strip_dirname(dirname)

# See if magnetism is "on"
magnetism = get_parameter(dirname, 'magnetism')

# SPECIFIC ARGS for etrace:
kw_default = dict({'the_file': None, 'xminmax': None, 'xmin': None, 'xmax': None, 'minmax': None, 'min': None, 'max': None, 'coords': None, 'ntot': 500, 'xiter': False, 'log': False, 'xvals': np.array([]), 'kw.nquadr': None, 'nquadlat': None, 'qvals': None, 'groupname': 'v', 'totsig': None, 'titles': None, 'justtot': False, 'notot': False, 'printerr': False, 'dpi': 300, 'vol': False, 'tint': False, 'tkappa': False, 'tavg': None})

# make figure kw
lineplot_fig_dimensions['margin_top_inches'] = 1.
kw_make_figure_default.update(lineplot_fig_dimensions)
kw_default.update(kw_make_figure_default)

# lineplot kw
kw_default.update(kw_lineplot_default)

# then override defaults
kw = update_dict(kw_default, clas)
kw_make_figure = update_dict(kw_make_figure_default, clas)
kw_lineplot = update_dict(kw_lineplot_default, clas)

# deal with desired quantities
if kw.qvals is None: # it's a quantity group
    qgroup = get_quantity_group(kw.groupname, magnetism)
    kw.qvals = qgroup['qvals']
    if kw.titles is None:
        kw.titles = qgroup['titles']
    if kw.totsig is None:
        kw.totsig = qgroup['totsig']
else:
    kw.qvals = make_array(kw.qvals)
    kw.titles = parse_quantities(kw.qvals)[1]
    kw.groupname = input("choose a groupname to save your plot\n to not save it, enter 'nosave': ")

fontsize = default_titlesize

# state what we're plotting
print ("plotting the following quantities:")
print ("qvals = ", kw.qvals)
print ("totsig = ", kw.totsig)

# deal with coords (if user wants minmax to only apply to certain subplots)
if not kw.coords is None:
    numpanels = len(kw.coords)//2
    acopy = np.copy(kw.coords)
    kw.coords = []
    for i in range(numpanels):
        kw.coords.append((acopy[2*i], acopy[2*i + 1]))

# get desired data file
dataname = 'G_Avgs_trace'
if kw.the_file is None:
    if not kw.nquadlat is None:
        dataname += '_kw.nquadlat%i' %nquadlat
    if not kw.nquadr is None:
        dataname += '_kw.nquadr%i' %nquadr
    kw.the_file = get_widest_range_file(clas0['datadir'], dataname)

print ('Getting data from ' + kw.the_file)
di = get_dict(kw.the_file)
vals = di['vals']
if kw.vol:
    vols = di['volumes']
    shapevols = list(np.shape(vols))
    newshape = [1,1] + shapevols
    vals = vals*vols.reshape(newshape)
rvals = di['rvals']
kw.nquadr = len(rvals) - 1
latvals = di['latvals']
kw.nquadlat = len(latvals) - 1
lut = di['lut']
times = di['times']
iters = di['iters']

# get the x axis
time_unit, time_label, rotation, simple_label = get_time_unit(dirname, tkappa=kw.tkappa)
if not kw.xiter:
    xaxis = np.copy(times)/time_unit
else:
    xaxis = np.copy(iters)

# set xminmax if not set by user
if kw.xminmax is None:
    # set xmin possibly
    if kw.xmin is None:
        kw.xmin = xaxis[0]
    # set xmax possibly
    if kw.xmax is None:
        kw.xmax = xaxis[-1]
    kw.xminmax = kw.xmin, kw.xmax
ixmin = np.argmin(np.abs(xaxis - kw.xminmax[0]))
ixmax = np.argmin(np.abs(xaxis - kw.xminmax[1]))

# Now shorten all the "x" arrays
xaxis = xaxis[ixmin:ixmax+1]
iters = iters[ixmin:ixmax+1]

# possibly time average data
if kw.tavg is None:
    print (buff_line)
    print("No time average: tavg = None")
else:
    vals, intervals = sliding_average_fancy(vals, xaxis, kw.tavg)
    print (buff_line)
    print ("Performing time average")
    print ("mean(tavg) = %1.3e" %np.mean(intervals))
    print ("std(tavg) = %1.3e" %np.std(intervals))

# deal with x axis, maybe thinning data
if np.all(kw.ntot == 'full'):
    print ('ntot = full')
    kw.ntot = len(times)
print ("ntot = %i" %kw.ntot)
print ("before thin_data: len(xaxis) = %i" %len(xaxis))
xaxis = thin_data(xaxis, kw.ntot)
#times = thin_data(times, kw.ntot)
iters = thin_data(iters, kw.ntot)
#vals = thin_data(vals, kw.ntot) # don't thin the data quite yet for this guy
print ("after thin_data: len(xaxis) = %i" %len(xaxis))

# now finally get the shape of the "vals" array
ntimes, nq, kw.nquadlat, nquadr = np.shape(vals)
ntimes = len(xaxis)
nplots = kw.nquadlat*nquadr

# create the figure dimensions
kw_make_figure.nplots = nplots
kw_make_figure.ncol = kw.nquadr
fig, axs, fpar = make_figure(**kw_make_figure)

# loop over different domains
for ilat in range(kw.nquadlat):
    for ir in range(kw.nquadr):
        vals_loc = vals[:, :, ilat, ir]
        ax = axs[ilat, ir]

        # get terms we want and plot them
        terms = []
        kw_lineplot.labels = []
        nterms = len(kw.qvals)
        for iterm in range(nterms):
            qval = kw.qvals[iterm]
            terms.append(vals_loc[:, lut[int(qval)]])
            #plt.plot(terms[-1])
            #plt.show()
            kw_lineplot.labels.append(kw.titles[iterm])

        # might also need the total of these terms (with some signature: totsig)

        if not kw.totsig is None:
            tot_term = np.zeros_like(terms[0])
            for iterm in range(nterms):
                tot_term += terms[iterm]*kw.totsig[iterm]
                # if totsig != 0,
                # also update the term itself with the proper sign
                if kw.totsig[iterm] != 0:
                    terms[iterm] *= kw.totsig[iterm]
            terms.append(tot_term)
            nterms += 1
            kw_lineplot.labels.append('sum')

        equation_sets = ['torque', 'teq', 'forcer', 'forcet', 'forcep', 'indr', 'indt', 'indp', 'meprodtotr', 'meprodtott', 'meprodtotp']
        if kw.groupname in equation_sets :
            if not kw.tint:
                # replace term with its time derivative
                terms[0] = drad(terms[0], times)
                kw_lineplot.labels[0] = '(d/dt)' + kw_lineplot.labels[0]

        if kw.tint: # integrate each term from the beginning
            terms[0] -= terms[0][0] # subtract off the initial value
            # so everything starts at 0
            for iterm in range(1,nterms): # now integrate the source terms in time
                #t0 = times[0]
                #terms[iterm] = indefinite_integral(terms[iterm], times, t0)
                terms[iterm] = cumtrap(terms[iterm], times, initial=0)
                kw_lineplot.labels[iterm] = r'$\int$' + kw_lineplot.labels[iterm]


        # name these in case we quantify their rms later
        ddt = terms[0]
        tot = terms[-1]
        difference = ddt - tot

        if kw.justtot:
            if 'meprod' in kw.groupname: # add induction + diffusion
                ddt *= (4*np.pi)
                induct = terms[1]
                difference = tot - ddt # need to recompute this now
                diffusion = terms[-2]
            if 'ind' in kw.groupname:
                induct = terms[-3]
                diffusion = terms[-2]
            if 'meprod' in kw.groupname or 'ind' in kw.groupname:
                terms = [ddt, induct, diffusion]
                kw_lineplot.labels = ['d/dt (LHS)', 'induction', 'diffusion']
            else:
                terms = [ddt, tot, difference]
                #kw_lineplot.labels = ['RHS - LHS', 'sum (RHS)', 'd/dt (LHS)']
                kw_lineplot.labels = [kw_lineplot.labels[0], kw_lineplot.labels[-1], kw_lineplot.labels[0] + ' - ' + kw_lineplot.labels[-1]]
            kw_lineplot.linestyles = ['-', '--', ':']
            kw_lineplot.colors = ['k', 'r', 'g']
            nterms = len(terms) # we reduced nterms
        elif kw.notot:
            terms = terms[1:-1]
            kw_lineplot.labels = kw_lineplot.labels[1:-1]
            nterms = len(terms) # we reduced nterms
        else: # make sure the last curve is dashed 
            #kw_lineplot.colors = ['k'] + color_order[1:nterms-1] + ['r']
            #kw_lineplot.linestyles = ['-'] + style_order[1:nterms-1] + ['--']
            kw_lineplot.linestyles = style_order[:nterms-1] + ['--']

        # now thin the data on the terms and times
        #times = thin_data(times[ixmin:ixmax+1], ntot)
        for iterm in range(nterms):
            terms[iterm] = thin_data(terms[iterm][ixmin:ixmax+1], kw.ntot)

        # now plot the terms
        if ilat == 0 and ir == 0:
            kw_lineplot.plotleg = True
        lineplot(xaxis, terms, ax, **kw_lineplot)

        if kw.groupname in equation_sets and kw.printerr:
            kw_lineplot.legfrac = 0.6

            # label the derivative ampltitude and sum amplitude
            quant = 'rms(d/dt) = %1.3e' %rms(ddt)
            quant += '\n' + 'rms(sum) = %1.3e' %rms(tot)
            quant += '\n' + 'rms(d/dt - sum) = %1.3e' %rms(difference)
            quant += '\n' + 'err = %1.3e' %(rms(difference)/rms(ddt))
            xmin, xmax = ax.get_xlim()
            dx = xmax - xmin
            ymin, ymax = ax.get_ylim()
            dy = ymax - ymin
            ax.text(xmin + 0.1*dx, ymin + 0.3*dy, quant, va='bottom', ha='left')

# Set some parameters defining all subplots
# x limits and label
axs[0, 0].set_xlim((kw.xminmax[0], kw.xminmax[1]))
if kw.xiter:
    axs[-1, 0].set_xlabel('iteration #')
else:
    axs[-1, 0].set_xlabel('time [' + time_label + ']')

# x titles
for ir in range(kw.nquadr):
    r1 = rvals[ir]
    r2 = rvals[ir+1]
    title = 'rad. range = [%.3f, %.3f]' %(r1, r2)
    axs[0, ir].set_title(title, fontsize=fontsize)

# y labels
for it in range(kw.nquadlat):
    lat1 = latvals[it]
    lat2 = latvals[it+1]
    axs[it, 0].set_ylabel('lat. range = [%.1f, %.1f]' %(lat1, lat2), fontsize=fontsize)

# main title
maintitle = dirname_stripped + '\nquadrant traces'
if kw.vol:
    maintitle += ' (volume-integrated'
else:
    maintitle += ' (volume-averaged'
if kw.tint:
    maintitle += ', time-integrated'
maintitle += ')'
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
    if kw.justtot:
        basename += '_justtot'
    elif kw.notot:
        basename += '_notot'
    if kw.tint:
        basename += '_tint'
    basename += clas0['tag']
    iter1, iter2 = get_iters_from_file(kw.the_file)
    savename = basename + '-' + str(iter1).zfill(8) + '_' + str(iter2).zfill(8) + '.png'

    plotdir = my_mkdir(clas0['plotdir'] + 'timetrace/')

    iter1, iter2 = get_iters_from_file(kw.the_file)
    savefile = plotdir + basename + '-' + str(iter1).zfill(8) + '_' + str(iter2).zfill(8) + '.png'
    if clas0['saveplot']:
        print ('saving figure at ' + savefile)
        fig.savefig(savefile, dpi=kw.dpi)

# Show the plot
if clas0['showplot']:
    plt.show()
plt.close()
