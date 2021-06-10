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
from timey_util import plot_timey

# Get CLAs
args = sys.argv
if not '--qvals' in args:
    args += ['--qvals', 'b'] # make default qvals = B field
clas0, clas = read_clas(args)
clas['plotcbar'] = False
dirname = clas0['dirname']
dirname_stripped = strip_dirname(dirname)

# get grid info
di_grid = get_grid_info(dirname)

datatype = 'timelat'
sampleaxis = di_grid['tt_lat']
axislabel = 'latitude (deg)'
samplelabel =  r'$r/R_\odot$' + ' = %.3f'
rad = False
lon = False
if 'rad' in clas:
    rad = True
    datatype = 'timerad'
    sampleaxis = di['rr']/rsun
    axislabel = r'$r/R_\odot$'
    samplelabel = 'lat = %.0f'
elif 'lon' in clas:
    lon = True
    clat = clas['clat']
    datatype = 'timelon_clatN%02i' %clat
    sampleaxis = di_grid['lons']
    axislabel = 'longitude (deg)'
    samplelabel = 'lat = %.0f ' +  r'$r/R_\odot$' + ' = %.3f'

dataname = datatype + clas0['tag']

# get data
if 'the_file' in clas: 
    the_file = clas['the_file']
else:
    the_file = get_widest_range_file(clas0['datadir'], dataname)

# Read in the data
print ('reading ' + the_file)
di = get_dict(the_file)
vals = di['vals']
times = di['times']
iters = di['iters']
samplevals_avail = di['samplevals']
qvals_avail = np.array(di['qvals'])

# baseline time unit
iter1, iter2 = get_iters_from_file(the_file)
time_unit, time_label, rotation, simple_label = get_time_unit(dirname)
times /= time_unit

# get raw traces of desired variables
terms = []
for qval in make_array(clas['qvals']):
    qind = np.argmin(np.abs(qvals_avail - qval))
    terms.append(vals[:, :, :, qind])

# set figure dimensions
sub_width_inches = 7.5
sub_height_inches = 2.0
margin_inches = 1./4.
sub_margin_bottom_inches = 1/2 # space for x-axis label + colorbar
margin_top_inches = 3/4
margin_left_inches = 1/2 # space for latitude label
nplots = len(terms)

# make plot
#fig, axs, fpar = make_figure(nplots=nplots, sub_width_inches=sub_width_inches, sub_height_inches=sub_height_inches, margin_left_inches=margin_left_inches, margin_top_inches=margin_top_inches, margin_bottom_inches=margin_bottom_inches)

# determine desired levels to plot
if not 'isamplevals' in clas:
    if not 'samplevals' in clas:
        isamplevals = np.array([0]) # just plot the top radius by default
    else: # get isamplevals from samplevals
        samplevals = clas['samplevals']
        if samplevals == 'all':
            isamplevals = np.arange(len(samplevals_avail))
        else:
            samplevals = make_array(samplevals)
            isamplevals = np.zeros_like(samplevals, dtype='int')
            for i in range(len(samplevals)):
                isamplevals[i] = np.argmin(np.abs(samplevals_avail - rsamplevals[i]))
else:
    isamplevals = make_array(clas['isamplevals'])

# Loop over the desired levels and save plots
for isampleval in isamplevals:
    sampleval = samplevals_avail[isampleval]
    print('plotting sampleval = %0.3f (i = %02i)' %(sampleval, isampleval))
   
    # Make appropriate file name to save
    savename = dataname + ('_%08i_%08i_' %(iter1, iter2)) +\
        ('sampleval%0.3f' %sampleval) + '.png'

    # make plot
    fig, axs, fpar = make_figure(nplots=nplots, ncol=1, sub_width_inches=sub_width_inches, sub_height_inches=sub_height_inches, margin_left_inches=margin_left_inches, margin_top_inches=margin_top_inches, sub_margin_bottom_inches=sub_margin_bottom_inches)

    for iplot in range(nplots):
        ax = axs[iplot, 0]
        field = terms[iplot][:, :, isampleval]
        plot_timey(field, times, sampleaxis, fig, ax, **clas)
                
        #  title the plot
        ax.set_title(clas['titles'][iplot], fontsize=default_titlesize)

    # Turn the x tick labels off for the top strips
    for iax in range(nplots):
        ax = axs[iax, 0]
        if iax < nplots - 1:
            ax.set_xticklabels([])
        if iax == nplots - 1:
            ax.set_xlabel('time (' + time_label + ')')
        if iax == nplots//2:
            # Label y-axis (latitude in degrees)
            ax.set_ylabel('latitude (deg)')

    # Put some useful information on the title
    if lon:
        sampleval = (clat, sampleval)
    maintitle = dirname_stripped + '\n' + (samplelabel %sampleval)
    if 'navg' in clas:
        navg = clas['navg']
        averaging_time = (times[-1] - times[0])/len(times)*navg
        maintitle += '\n' + ('t_avg = %.1f Prot' %averaging_time)
    else: 
        maintitle += '\nt_avg = none'
    fig.text(fpar['margin_left'], 1 - fpar['margin_top'], maintitle, fontsize=default_titlesize, ha='left', va='bottom')

    # Save the plot
    if clas0['saveplot']:
        # save the figure
        plotdir = my_mkdir(clas0['plotdir'] + 'timelat/')
        print ('Saving the time-latitude plot at ')
        print (plotdir + savename)
        print ("=======================================")
        #plt.savefig(plotdir + savename, dpi=200)

    # Show the plot if only plotting at one latitude
    #if clas0['showplot'] and len(isamplevals) == 1:
    plt.show()
    plt.close()
