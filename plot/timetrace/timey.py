# Author: Loren Matilsky
# Date created: 03/02/2019
import matplotlib.pyplot as plt
import numpy as np
import sys, os
sys.path.append(os.environ['raco'])
sys.path.append(os.environ['raco'] + '/quantities_util')
sys.path.append(os.environ['rapl'])
sys.path.append(os.environ['rapl'] + '/timetrace')
from common import *
from cla_util import *
from plotcommon import *
from timey_util import *
from varprops import get_quantity_group

# Set fontsize
fontsize = default_titlesize

# Read command-line arguments (CLAs)
args = sys.argv
clas0, clas = read_clas(args)
dirname = clas0['dirname']
dirname_stripped = strip_dirname(dirname)
# See if magnetism is "on"
magnetism = clas0['magnetism']

# defaults
kwargs_default = dict({'rad': False, 'groupname': 'v', 'sampletag': '', 'the_file': None, 'isamplevals': np.array([0]), 'samplevals': None, 'rvals': None, 'qvals': 'all', 'ntot': 500, 'tavg': None, 'prepend': False, 'sub': False, 'xminmax': None, 'xmin': None, 'xmax': None, 'ntheta': None, 'tdt': False})

# also need make figure kwargs
make_figure_kwargs_default.update(timey_fig_dimensions)
kwargs_default.update(make_figure_kwargs_default)

# of course, also need plot_timey kwargs
kwargs_default.update(plot_timey_kwargs_default)

# check for bad keys
find_bad_keys(kwargs_default, clas, clas0['routinename'], justwarn=True)

# overwrite defaults
kw = update_dict(kwargs_default, clas)
kw_plot_timey = update_dict(plot_timey_kwargs_default, clas)
kw_make_figure = update_dict(make_figure_kwargs_default, clas)

# might need two colorbars
if not kw.ycut is None:  # need room for two colorbars
    kw_make_figure.sub_margin_right_inches *= 2
    kw_make_figure.margin_top_inches += 1/4

# baseline time unit
time_unit, time_label, rotation, simple_label = get_time_unit(dirname, tdt=kw.tdt)

# get grid info
di_grid = get_grid_info(dirname, ntheta=kw.ntheta)

if kw.rad:
    datatype = 'timerad'
    plotlabel = 'time-radius trace'
    yaxis = di_grid['rr']
    axislabel = 'radius'
    samplefmt = lat_fmt
    samplename = 'latval'
else:
    datatype = 'timelat'
    plotlabel = 'time-latitude trace'
    yaxis = di_grid['tt_lat']
    axislabel = 'latitude (deg)'
    samplefmt = '%1.3e'
    samplename = 'rval'

dataname = datatype + '_' + kw.groupname
if len(kw.sampletag) > 0:
    dataname += '_' + kw.sampletag

# get data
if kw.the_file is None:
    kw.the_file = get_widest_range_file(clas0['datadir'] +\
            datatype + '/', dataname)

# Read in the data
print ('reading ' + kw.the_file)
di = get_dict(kw.the_file)
vals = di['vals']
times = di['times']
iters = di['iters']
qvals_avail = np.array(di['qvals'])
samplevals_avail = di['samplevals']

# time unit
times /= time_unit

# set xminmax if not set by user
if kw.xminmax is None:
    # set xmin possibly
    if kw.xmin is None:
        kw.xmin = np.min(times)
    # set xmax possibly
    if kw.xmax is None:
        kw.xmax = np.max(times)
    kw.xminmax = kw.xmin, kw.xmax

ixmin = np.argmin(np.abs(times - kw.xminmax[0]))
ixmax = np.argmin(np.abs(times - kw.xminmax[1]))

# Now shorten all the "x" arrays
times = times[ixmin:ixmax + 1]
iters = iters[ixmin:ixmax + 1]
vals = vals[ixmin:ixmax + 1]

# determine desired levels to plot

# can control samplevals with rvals for time-latitude traces
if not kw.rad and not kw.rvals is None:
    kw.samplevals = kw.rvals

if not kw.samplevals is None: # isamplevals being set indirectly
    # check for special 'all' option
    if isall(kw.samplevals):
        kw.isamplevals = np.arange(len(samplevals_avail))
    else:
        kw.isamplevals = inds_from_vals(samplevals_avail, kw.samplevals)

# determine desired quantities to plot
if isall(kw.qvals):
    kw.qvals = qvals_avail

terms = []
for qval in kw.qvals:
    qind = np.argmin(np.abs(qvals_avail - qval))
    terms.append(vals[:, :, :, qind])

# determine titles
if kw.titles is None:
    kw.titles = parse_quantities(kw.qvals)[1]
    if kw.sub:
        for iplot in range(len(kw.titles)):
            kw.titles[iplot] += ' (sub. temp. mean)'

# maybe add some more quantities, dependent on qgroup
qgroup = get_quantity_group(kw.groupname, magnetism)
totsig = qgroup.totsig

if kw.groupname in ['torque']: # add two more terms at the top
    # time derivative of "field term"
    # and sum of force terms --- they should basically match up
    # but may not depending on sampling frequency
    d_dt = np.gradient(terms[0], times, axis=0)/time_unit
    # remember "times" is now normalized by time_unit
    the_sum = np.zeros_like(terms[1])
    nterms = len(terms)
    for iterm in range(1,nterms):
        the_sum += terms[iterm]*totsig[iterm]
        if totsig[iterm] < 0:
            terms[iterm] *= -1 # might as well also just plot the negative
            kw.titles[iterm] = '-' + kw.titles[iterm]
        
    if kw.sub:
        titletag = ' (sub. temp. mean)'
    else:
        titletag = ''
    title1 = 'd' + kw.titles[0] + '/dt' + titletag
    title2 = 'sum of force terms' + titletag

    # now insert these terms (and two titles)
    terms.insert(1, the_sum) # the_sum is now the second term
    terms.insert(1, d_dt) # d_dt is now the second term; the_sum is the third
    kw.titles.insert(1, title2)
    kw.titles.insert(1, title1)

# Will probably thin data --- if we do, need to keep track of 
# "thinned" times and iters
if kw.ntot == 'full':
    print(buff_line)
    print ("ntot = %i (full time series)" %kw.ntot)
    times_thin = times
    iters_thin = iters
else:
    print(buff_line)
    print ("ntot = %i" %kw.ntot)
    print ("before thin_data: len(times) = %i" %len(times))
    times_thin = thin_data(times, kw.ntot)
    iters_thin = thin_data(iters, kw.ntot)
    print ("after thin_data: len(times) = %i" %len(times))


# Loop over the desired levels and save plots
kw.isamplevals = make_array(kw.isamplevals) # needs to be array
firsttime = True # keep track if we've plotted any panels yet
# for printing purposes
for isampleval in kw.isamplevals:
    sampleval = samplevals_avail[isampleval]

    # set some labels 
    samplelabel = samplename + ' = ' + (samplefmt %sampleval)
    if kw.rad: # label things by colat (not lat) to be in sequential order
        position_tag = ('_colatval' + lon_fmt) %(90.0 - sampleval)
    else:
        position_tag = '_' + samplename + (samplefmt %sampleval)

    # say we are about to plot
    print(buff_line)
    print('plotting ' + samplelabel + ' (i = %02i)' %isampleval)
   
    # make plot
    nplots = kw_make_figure.nplots = len(terms)
    kw_make_figure.ncol = 1
    fig, axs, fpar = make_figure(**kw_make_figure)


    for iplot in range(nplots):
        ax = axs[iplot, 0]
        if kw.rad:
            field = terms[iplot][:, isampleval, :]
        else:
            field = terms[iplot][:, :, isampleval]

        # possibly time average data
        if kw.tavg is None:
            if firsttime:
                print (buff_line)
                print("No time average: tavg = None")
        else:
            field, intervals = sliding_average(field, times, kw.tavg)
            if firsttime:
                print (buff_line)
                print ("Performing time average, tavg = %.2f Prot" %np.mean(intervals))
                print ("sigma(tavg) = %.3f Prot" %np.std(intervals))

        # possibly subtract temporal mean
        if kw.sub: # full Omega (no subtraction)
            nx, ny = np.shape(field)
            tempmean = np.mean(field, axis=0).reshape((1, ny))
            field -= tempmean
            if firsttime:
                print ("Subtracting temporal mean from data")
        else:
            if firsttime:
                print ("Full data: no subtraction of temporal mean")

        # probably thin data
        if not kw.ntot == 'full':
            field = thin_data(field, kw.ntot)

        print ("plotting panel %02i of %02i" %(iplot, nplots))
        plot_timey(field, times_thin, yaxis, fig, ax, **kw_plot_timey)
                
        #  title the plot
        ax.set_title(kw.titles[iplot], fontsize=fontsize)

        # Turn the x tick labels off for the top strips
        #if iplot < nplots - 1:
        #    ax.set_xticklabels([])
        # Put time label on bottom strip        
        if iplot == nplots - 1:
            ax.set_xlabel('time (' + time_label + ')', fontsize=fontsize)
        # Put ylabel on middle strip
        if iplot == nplots//2:
            ax.set_ylabel(axislabel, fontsize=fontsize)

        # Only do the print messages once
        firsttime = False

    # Put some useful information on the title
    maintitle = dirname_stripped + '\n' +\
            plotlabel + '\n' +\
            'groupname = ' + kw.groupname + '\n' +\
            samplelabel
    if kw.tavg is None:
        maintitle += '\ntavg = none'
    else:
        maintitle += '\n' + ("tavg = %.2f Prot, sigma(tavg) = %.3f Prot"\
                %(np.mean(intervals), np.std(intervals)))

    maintitle += '\nm=0 (lon. avg.)'
    if not kw.ycut is None:
        maintitle += '\nycut = %1.3e' %kw.ycut

    margin_x = fpar['margin_left'] + fpar['sub_margin_left']
    margin_y = default_margin/fpar['height_inches']
    fig.text(margin_x, 1 - margin_y, maintitle,\
             ha='left', va='top', fontsize=default_titlesize)

    # Save the plot
    if clas0['saveplot']:
        # Make appropriate file name to save

        # save the figure
        iter1, iter2 = get_iters_from_file(kw.the_file)
        basename = dataname + '-%08i_%08i' %(iter1, iter2)
        plotdir = my_mkdir(clas0['plotdir'] + '/' + datatype + clas0['tag'])
        if kw.lon and not om is None:
            basename += '_om%.0f' %om
        savename = basename + position_tag + '.png'
        if kw.prepend:
            savename = dirname_stripped + '_' + savename
        print ("saving", plotdir + '/' + savename)
        plt.savefig(plotdir + '/' + savename, dpi=200)

    # Show the plot if only plotting at one latitude
    if clas0['showplot'] and len(kw.isamplevals) == 1:
        plt.show()
    else:
        plt.close()
    print ("=======================================")
