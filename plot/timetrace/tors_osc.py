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
kwargs_default = dict({'the_file': None, 'ntot': 500, 'tavg': None, 'clat': 10, 'dlat': 0, 'om': None, 'rad': False, 'lon': False, 'isamplevals': np.array([0]), 'samplevals': None, 'rvals': None, 'groupname': 'v', 'sub': False, 'prepend': False, 'ntheta': None, 'tdt': False})

# also need make figure kwargs
timey_fig_dimensions['margin_top_inches'] = 1.0
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
time_unit, time_label, rotation, simple_label = get_time_unit(dirname, tdt=kw.tdt)

# get grid info
di_grid = get_grid_info(dirname, ntheta=kw.ntheta)
sint = di_grid['sint']
rr = di_grid['rr']
nr = di_grid['nr']
nt = di_grid['nt']

if kw.rad:
    datatype = 'timerad'
    plotlabel = 'torsional oscillation: time-radius'
    yaxis = di_grid['rr']
    axislabel = 'radius'
    samplefmt = lat_fmt
    samplename = 'latval'
else:
    datatype = 'timelat'
    plotlabel = 'torsional oscillation: time-latitude'
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
iter1, iter2 = get_iters_from_file(kw.the_file)
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
times = times[ixmin:ixmax+1]
iters = iters[ixmin:ixmax+1]
vals = vals[ixmin:ixmax+1, :]

# maybe thin data
if not kw.ntot == 'full':
    print ("ntot = %i" %kw.ntot)
    print ("before thin_data: len(times) = %i" %len(times))
    times = thin_data(times, kw.ntot)
    iters = thin_data(iters, kw.ntot)
    vals = thin_data(vals, kw.ntot)
    print ("after thin_data: len(times) = %i" %len(times))

# get differential rotation
qind = np.argmin(np.abs(qvals_avail - 3))
vp = vals[:, :, :, qind]
if kw.rad:
    sampletheta = np.pi/2.0 - samplevals_avail*np.pi/180.0
    xx = rr.reshape((1, 1, nr))*\
            np.sin(sampletheta).reshape((1, len(sampletheta), 1))
else:
    xx = samplevals_avail.reshape((1, 1, len(samplevals_avail)))*\
            sint.reshape((1, nt, 1))

# frame rate
eq = get_eq(dirname)
Om0 = 2*np.pi/eq.trot

# differential rotation in the rotating frame. 
Om = vp/xx

if not kw.sub: # full Omega (no subtraction)
    subplottitle = 'full rotation: ' + r'$(\Omega - \Omega_0)/\Omega_0$'
else:
    dummy, n1, n2 = np.shape(Om)
    tempmean = np.mean(Om, axis=0).reshape((1, n1, n2))
    Om -= tempmean
    subplottitle = 'residual rotation: ' + r'$(\Omega - \langle\Omega\rangle_t)/\Omega_0$'

# determine desired levels to plot

# can control samplevals with rvals for time-latitude traces
if not kw.rad and not kw.rvals is None:
    kw.samplevals = kw.rvals

if not kw.samplevals is None: # isamplevals being set indirectly
    # check for special 'all' option
    if isall(kw.samplevals):
        kw.isamplevals = np.arange(len(samplevals_avail))
    else:
        kw.samplevals = make_array(kw.samplevals)
        kw.isamplevals = np.zeros_like(kw.samplevals, dtype='int')
        for i in range(len(kw.samplevals)):
            kw.isamplevals[i] = np.argmin(np.abs(samplevals_avail - kw.samplevals[i]))


# Loop over the desired levels and save plots
firsttime = True # keep track if we've plotted any panels yet
# for printing purposes
kw.isamplevals = make_array(kw.isamplevals) # needs to be array
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

    # set up figure
    fig, axs, fpar = make_figure(**kw_make_figure)
    ax = axs[0, 0]

    # get desired "field" to plot
    if kw.rad:
        field = Om[:, isampleval, :]/Om0
    else:
        field = Om[:, :, isampleval]/Om0

    # possibly time average data
    if kw.tavg is None:
        if firsttime:
            print (buff_line)
            print("No time average: tavg = None")
    else:
        field, intervals = sliding_average(field, times, kw.tavg)
        if firsttime:
            print (buff_line)
            print ("Performing time average, tavg = %.2f t_rot" %np.mean(intervals))
            print ("sigma(tavg) = %.3f t_rot" %np.std(intervals))

    plot_timey(field, times, yaxis, fig, ax, **kw_plot_timey)
            
    #  title the plot
    ax.set_title(subplottitle, fontsize=fontsize)

    # time label
    ax.set_xlabel('time (' + time_label + ')', fontsize=fontsize)
    # ylabel on middle strip
    ax.set_ylabel(axislabel, fontsize=fontsize)

    # Put some useful information on the title
    maintitle = dirname_stripped + '\n' +\
            plotlabel + '\n' +\
            samplelabel

    if kw.tavg is None:
        maintitle += '\ntavg = none'
    else:
        maintitle += '\n' + ("tavg = %.2f t_rot, sigma(tavg) = %.3f t_rot"\
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
        basename = datatype + '_torsosc-%08i_%08i' %(iter1, iter2)
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
    print (buff_line)
