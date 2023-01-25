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
from rayleigh_diagnostics import Shell_Slices

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
kwargs_default = dict({'irvals': np.array([0]), 'rvals': None, 'ntot': 500, 'groupname': 'b', 'mmax': 10, 'mval': 1, 'imag': False, 'abs': False, 'qvals': None, 'rad': False})

# also need make figure kwargs
timey_fig_dimensions['margin_top_inches'] += 1/4
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

# check if we want the real or imaginary vals
if kw.imag:
    part = 'imag'
elif kw.abs:
    part = 'abs'
    kw_plot_timey.posdef = True
else:
    part = 'real'

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

# mval
mval = kw.mval

# get the rvals we want
slice_info = get_sliceinfo(dirname)
rvals_avail = slice_info.samplevals

irvals = kw.irvals
if not kw.rvals is None: # irvals haven't been set directly
    if isall(kw.rvals):
        irvals = np.arange(slice_info.nsamplevals)
    else:
        kw.rvals = make_array(kw.rvals)
        irvals = inds_from_vals(rvals_avail, kw.rvals)

# and the qvals
qvals = kw.qvals

# maybe qvals can be all (everything available)
if isall(qvals): # but probably won't use this option here ... would 
    # make too many panels
    qvals = np.sort(slice_info.qv)

if qvals is None: # it's a quantity group
    groupname = kw.groupname
    qgroup = get_quantity_group(groupname, magnetism)
    qvals = qgroup['qvals']
    titles = qgroup['titles']
else:
    titles = []
    for qval in qvals:
        titles.append(str(qval))
    groupname = input("choose a groupname to save your plot: ")

irvals = make_array(irvals)
#qvals = make_array(qvals)

# mmax (if needed)
mmax = kw.mmax


print (buff_line)
print ("plotting time-latitude trace for mval = %03i" %mval)
print ("irvals = " + arr_to_str(irvals, "%i"))
print ("rvals  = " + arr_to_str(slice_info.samplevals[irvals], "%1.3e"))
print ("qvals  = " + arr_to_str(qvals, "%i"))
print (buff_line)

# Loop over the desired levels and save plots
firstplot = True
for irval in irvals:
    # for each plot, collect the terms (qvals) we want
    terms = []

    # must read in each data file separately
    count = 0
    for qval in qvals:
        dataname = ('mtrace_qval%04i_irval%02i' %(qval, irval))
        print ("dataname = ", dataname)

        # get data
        the_file = get_widest_range_file(clas0['datadir'] +\
                    '/mtrace_mmax%03i/' %mmax, dataname)

        # Read in the data
        print ('reading ' + the_file)
        di = get_dict(the_file)
        if part == 'imag':
            vals = np.imag(di['vals'][:, mval, :])
        elif part == 'abs':
            vals = np.abs(di['vals'][:, mval, :])
        else:
            vals = np.real(di['vals'][:, mval, :])
        times = di['times']
        iters = di['iters']

        # time range
        times /= time_unit

        # maybe thin data
        if not kw.ntot == 'full':
            if count == len(qvals) - 1: # last one
                print (buff_line)
                print ("ntot = %i" %kw.ntot)
                print ("before thin_data: len(times) = %i" %len(times))

            times = thin_data(times, kw.ntot)
            iters = thin_data(iters, kw.ntot)
            vals = thin_data(vals, kw.ntot)

            if count == len(qvals) - 1: # last one
                print ("after thin_data: len(times) = %i" %len(times))
                print (buff_line)
        
        terms.append(vals)
        count += 1

    # set some labels 
    samplelabel = samplename + ' = ' + (samplefmt %sampleval)
    if kw.rad: # label things by colat (not lat) to be in sequential order
        position_tag = ('_colatval' + lon_fmt) %(90.0 - sampleval)
    else:
        position_tag = '_' + samplename + (samplefmt %sampleval)

    axislabel = 'latitude (deg)'
    sampleval = slice_info.samplevals[irval]

    # Put some useful information on the title
    maintitle = dirname_stripped + '\n' +\
            plotlabel + '\n' +\
            'groupname = ' + kw.groupname + '\n' +\
            samplelabel

    maintitle += '\nmval=%03i' %mval
    if part == 'imag':
        maintitle += '\nimaginary part'
    elif part == 'abs':
        maintitle += '\nabsolute magnitude'
    else:
        maintitle += '\nreal part'
    if kw.navg is None:
        maintitle += '\nt_avg = none'
    else:
        averaging_time = (times[-1] - times[0])/len(times)*kw.navg
        maintitle += '\n' + ('t_avg = %.1f Prot' %averaging_time)

    print('plotting rval = %1.3e (i = %02i)' %(sampleval, irval))

    # make plot
    nplots = kw_make_figure.nplots = len(terms)
    kw_make_figure.ncol = 1
    fig, axs, fpar = make_figure(**kw_make_figure)

    for iplot in range(nplots):
        ax = axs[iplot, 0]
        field = terms[iplot]
        plot_timey(field, times, yaxis, fig, ax, **kw_plot_timey)
                
        #  title the plot
        ax.set_title(titles[iplot], fontsize=fontsize)

        # Put time label on bottom strip        
        if iplot == nplots - 1:
            ax.set_xlabel('time (' + time_label + ')', fontsize=fontsize)
        # Put ylabel on middle strip
        if iplot == nplots//2:
            ax.set_ylabel(axislabel, fontsize=fontsize)

    fig.text(fpar['margin_left'], 1 - fpar['margin_top'], maintitle, fontsize=fontsize, ha='left', va='bottom')

    # Save the plot
    if clas0['saveplot']:
        # Make appropriate file name to save

        # save the figure
        iter1, iter2 = get_iters_from_file(the_file)
        basename = 'mtracetimelat_' + groupname + '-%08i_%08i' %(iter1, iter2)
        plotdir = my_mkdir(clas0['plotdir'] +\
                '/mtracetimelat_mval%03i' %mval + clas0['tag'])
        savename = basename + position_tag + '_' + part + '.png'
        print ("saving", plotdir + '/' + savename)
        plt.savefig(plotdir + '/' + savename, dpi=200)

    # Show the plot if only plotting at one latitude
    if clas0['showplot'] and len(irvals) == 1:
        plt.show()
    else:
        plt.close()
    print (buff_line)
