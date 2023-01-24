# Author: Loren Matilsky
# Date created: 03/02/2019
import matplotlib.pyplot as plt
import numpy as np
import sys, os
sys.path.append(os.environ['raco'])
sys.path.append(os.environ['rapl'])
sys.path.append(os.environ['rapl'] + '/timetrace')
from common import *
from cla_util import *
from plotcommon import *
from timey_util import *

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
kwargs_default = dict({'the_file': None, 'ntot': 500, 'rad': False, 'isamplevals': np.array([0]), 'samplevals': None, 'groupname': 'b'})

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

# add in groupname keys
kw.update(get_quantity_group(kw.groupname, magnetism))

# user may have wanted to change some groupname keys
kw = update_dict(kw, clas)

# baseline time unit
time_unit, time_label, rotation, simple_label = get_time_unit(dirname)

# get grid info
di_grid = get_grid_info(dirname)

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

# time range
times /= time_unit

# maybe thin data
if not kw.ntot == 'full':
    print ("ntot = %i" %kw.ntot)
    print ("before thin_data: len(times) = %i" %len(times))
    times = thin_data(times, kw.ntot)
    iters = thin_data(iters, kw.ntot)
    vals = thin_data(vals, kw.ntot)
    print ("after thin_data: len(times) = %i" %len(times))

# get raw traces of desired variables
terms = []
for qval in kw.qvals:
    qind = np.argmin(np.abs(qvals_avail - qval))
    terms.append(vals[:, :, :, qind])

# determine desired levels to plot
if not kw.samplevals is None: # isamplevals being set indirectly
    # check for special 'all' option
    if kw.samplevals == 'all':
        kw.isamplevals = np.arange(len(samplevals_avail))
    else:
        kw.samplevals = make_array(kw.samplevals)
        kw.isamplevals = np.zeros_like(kw.samplevals, dtype='int')
        for i in range(len(kw.samplevals)):
            kw.isamplevals[i] = np.argmin(np.abs(samplevals_avail - kw.samplevals[i]))

# Loop over the desired levels and save plots
for isampleval in kw.isamplevals:
    sampleval = samplevals_avail[isampleval]

    # set some labels 
    samplelabel = samplename + ' = ' + (samplefmt %sampleval)
    if kw.rad: # label things by colat (not lat) to be in sequential order
        position_tag = ('_colatval' + lon_fmt) %(90.0 - sampleval)
    else:
        position_tag = '_' + samplename + (samplefmt %sampleval)


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
        plot_timey(field, times, yaxis, fig, ax, **kw_plot_timey)
                
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

    # Put some useful information on the title
    maintitle = dirname_stripped + '\n' +\
            plotlabel + '\n' +\
            'groupname = ' + kw.groupname + '\n' +\
            samplelabel
    if kw.navg is None:
        maintitle += '\nt_avg = none'
    else:
        averaging_time = (times[-1] - times[0])/len(times)*kw.navg
        maintitle += '\n' + ('t_avg = %.1f Prot' %averaging_time)

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
        print ("saving", plotdir + '/' + savename)
        plt.savefig(plotdir + '/' + savename, dpi=200)

    # Show the plot if only plotting at one latitude
    if clas0['showplot'] and len(kw.isamplevals) == 1:
        plt.show()
    else:
        plt.close()
    print ("=======================================")
