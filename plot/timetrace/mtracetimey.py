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

# Get the Rayleigh data directory
radatadir = dirname + '/Shell_Slices/'

# Get all the file names in datadir and their integer counterparts
file_list, int_file_list, nfiles = get_file_lists(radatadir, args)

# read first file for some metadata
a0 = Shell_Slices(radatadir + file_list[0], '')

# defaults
kwargs_default = dict({'ntot': 500, 'groupname': 'b', 'irvals': np.array([0]), 'rvals': None, 'mmax': 10, 'mval': 1})

kwargs_default.update(plot_timey_kwargs_default)

# check for bad keys
find_bad_keys(kwargs_default, clas, clas0['routinename'], justwarn=True)

# overwrite defaults
kw = update_dict(kwargs_default, clas)
# add in groupname keys
kw.update(get_quantity_group(kw.groupname, magnetism))

# user may have wanted to change some groupname keys
kw = update_dict(kw, clas)
kw_plot_timey = update_dict(plot_timey_kwargs_default, clas)

# baseline time unit
time_unit, time_label, rotation, simple_label = get_time_unit(dirname)

# get grid info
di_grid = get_grid_info(dirname)

sampleaxis = di_grid['tt_lat']
mval = kw.mval

# get the rvals we want
irvals = kw.irvals
if not kw.rvals is None: # irvals haven't been set directly
    if isall(kw.rvals):
        irvals = np.arange(a0.nr)
    else:
        irvals = np.zeros_like(kw.rvals, dtype='int')
        for i in range(len(kw.rvals)):
            irvals[i] = np.argmin(np.abs(a0.radius/rsun - kw.rvals[i]))

qvals = kw.qvals
# everything must be array
irvals = make_array(irvals)
#qvals = make_array(qvals)

# mmax (if needed)
mmax = kw.mmax

# set figure dimensions
sub_width_inches = 7.5
sub_height_inches = 2.0
margin_bottom_inches = 1/2 # space for x-axis label
margin_top_inches = 1
if kw.lon:
    margin_top_inches =  1 + 1/4
margin_left_inches = 5/8 # space for latitude label
margin_right_inches = 7/8 # space for colorbar
if 'ycut' in clas:
    margin_right_inches *= 2

nplots = len(qvals)

print (buff_line)
print ("plotting time-latitude trace for mval = %03i" %mval)
print ("irvals = ", irvals)
print ("r/rsun = ", a0.radius[irvals]/rsun)
print ("qvals = ", qvals)
print (buff_line)

# Loop over the desired levels and save plots
firstplot = True
for irval in irvals:
    # for each plot, collect the terms (qvals) we want
    terms = []

    # must read in each data file separately
    for qval in qvals:
        dataname = ('mtrace_qval%04i_irval%02i' %(qval, irval)) +\
                clas0['tag']

        # get data
        the_file = get_widest_range_file(clas0['datadir'] +\
                    '/mtrace_mmax%03i/' %mmax, dataname)

        # Read in the data
        print ('reading ' + the_file)
        di = get_dict(the_file)
        vals = np.real(di['vals'][:, mval, :])
        times = di['times']
        iters = di['iters']

        # time range
        iter1, iter2 = get_iters_from_file(the_file)
        times /= time_unit

        # maybe thin data
        if not kw.ntot == 'full':
            if firstplot:
                print ("ntot = %i" %kw.ntot)
                print ("before thin_data: len(times) = %i" %len(times))

            times = thin_data(times, kw.ntot)
            iters = thin_data(iters, kw.ntot)
            vals = thin_data(vals, kw.ntot)

            if firstplot:
                print ("after thin_data: len(times) = %i" %len(times))
                firstplot = False
        
        terms.append(vals)

    # set some labels 
    axislabel = 'latitude (deg)'
    rval = a0.radius[irval]/rsun
    samplelabel =  r'$r/R_\odot$' + ' = %.3f' %rval
    position_tag = '_rval%.3f' %rval

    # Put some useful information on the title
    maintitle = dirname_stripped 
    maintitle += '\n' + samplelabel
    maintitle += '\nmval=%03i' %mval
    if kw.navg is None:
        maintitle += '\nt_avg = none'
    else:
        averaging_time = (times[-1] - times[0])/len(times)*kw.navg
        maintitle += '\n' + ('t_avg = %.1f Prot' %averaging_time)

    print('plotting rval = %0.3f (i = %02i)' %(rval, irval))

    # make plot
    fig, axs, fpar = make_figure(nplots=nplots, ncol=1, sub_width_inches=sub_width_inches, sub_height_inches=sub_height_inches, margin_left_inches=margin_left_inches, margin_right_inches=margin_right_inches, margin_top_inches=margin_top_inches, margin_bottom_inches=margin_bottom_inches)

    for iplot in range(nplots):
        ax = axs[iplot, 0]
        field = terms[iplot]
        plot_timey(field, times, sampleaxis, fig, ax, **kw_plot_timey)
                
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

    fig.text(fpar['margin_left'], 1 - fpar['margin_top'], maintitle, fontsize=fontsize, ha='left', va='bottom')

    # Save the plot
    if clas0['saveplot']:
        # Make appropriate file name to save

        # save the figure
        if kw.groupname is None:
            basename = dataname
        else:
            basename = 'mtrace_' + kw.groupname + '_irval%02i' %irval
        basename += '-%08i_%08i' %(iter1, iter2)
        plotdir = my_mkdir(clas0['plotdir'] + '/timelat_mval%03i' %mval)
        savename = basename + clas0['tag'] + position_tag + '.png'
        print ("saving", plotdir + '/' + savename)
        plt.savefig(plotdir + '/' + savename, dpi=200)

    # Show the plot if only plotting at one latitude
    if clas0['showplot'] and len(irvals) == 1:
        plt.show()
    else:
        plt.close()
    print ("=======================================")
