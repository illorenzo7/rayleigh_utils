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

# Set fontsize
fontsize = default_titlesize

# Get CLAs
args = sys.argv
if not '--qvals' in args:
    args += ['--qvals', 'b'] # make default qvals = B field
clas0, clas = read_clas(args)
#clas['plotcbar'] = False
dirname = clas0['dirname']
dirname_stripped = strip_dirname(dirname)

if 'ntot' in clas:
    ntot = clas['ntot']
else:
    ntot = 2000

# baseline time unit
time_unit, time_label, rotation, simple_label = get_time_unit(dirname)

# get grid info
di_grid = get_grid_info(dirname)

datatype = 'timelat'
sampleaxis = di_grid['tt_lat']
rad = False
lon = False
shav = False
if 'rad' in clas:
    rad = True
    datatype = 'timerad'
    sampleaxis = di_grid['rr']/rsun
elif 'lon' in clas:
    lon = True
    if 'clat' in clas:
        clat = clas['clat']
    else:
        clat = 10
    if 'dlat' in clas:
        dlat = clas['dlat']
    else:
        dlat = 0
    if 'om' in clas: # plot in frame rotating at om [nHz]
        om = clas['om']
        om0 = 1/time_unit*1e9 # frame rate, nHz
    else:
        om = None

    datatype = 'timelon_clat' + lat_format(clat) + '_dlat%03.0f' %dlat
    sampleaxis = di_grid['lons']
elif 'shav' in clas:
    shav = True
    datatype = 'timeshav'
    sampleaxis = di_grid['rr']/rsun
    # just a loop placeholder...
    isamplevals = [0]

dataname = datatype + clas0['tag']
print ('dataname = ', dataname)
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
qvals_avail = np.array(di['qvals'])
if not shav:
    samplevals_avail = di['samplevals']

# time range
iter1, iter2 = get_iters_from_file(the_file)
times /= time_unit

# Subtract DR, if desired
if lon:
    if not om is None:
        print ("plotting in rotating frame om = %.1f nHz" %om)
        print ("compare this to frame rate    = %.1f nHz" %om0)
        phi_deflections = (times*(om - om0)) % 1 # between zero and one
        nphi = len(sampleaxis)
        for it in range(len(times)):
            phi_deflection = phi_deflections[it]
            nroll = int(phi_deflection*nphi)
            vals[it] = np.roll(vals[it], -nroll, axis=0)

# maybe thin data
if not ntot == 'full':
    print ("ntot = %i" %ntot)
    print ("before thin_data: len(times) = %i" %len(times))
    times = thin_data(times, ntot)
    iters = thin_data(iters, ntot)
    vals = thin_data(vals, ntot)
    print ("after thin_data: len(times) = %i" %len(times))

# get raw traces of desired variables
terms = []
for qval in make_array(clas['qvals']):
    qind = np.argmin(np.abs(qvals_avail - qval))
    terms.append(vals[:, :, :, qind])

# set figure dimensions
sub_width_inches = 7.5
sub_height_inches = 2.0
margin_bottom_inches = 1/2 # space for x-axis label
margin_top_inches = 3/4
if lon:
    margin_top_inches =  1 + 1/4
margin_left_inches = 5/8 # space for latitude label
margin_right_inches = 7/8 # space for colorbar
if 'ycut' in clas:
    margin_right_inches *= 2
nplots = len(terms)

# determine desired levels to plot
if not shav:
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
                    isamplevals[i] = np.argmin(np.abs(samplevals_avail - samplevals[i]))
    else:
        isamplevals = make_array(clas['isamplevals'])

# Loop over the desired levels and save plots
for isampleval in isamplevals:
    if not shav:
        sampleval = samplevals_avail[isampleval]

    # set some labels 
    axislabel = 'latitude (deg)'
    samplelabel =  r'$r/R_\odot$' + ' = %.3f' %sampleval
    position_tag = '_rval%.3f' %sampleval
    if rad:
        axislabel = r'$r/R_\odot$'
        samplelabel = 'lat = ' + lat_format(sampleval)
        position_tag = '_lat' + lat_format(sampleval)
    elif lon:
        axislabel = 'longitude (deg)'
        samplelabel = 'clat = ' + lat_format(clat) + '\n' +  r'$r/R_\odot$' + ' = %.3f' %sampleval
        if not om is None:
            samplelabel += '\n' + (r'$\Omega_{\rm{frame}}$' + ' = %.1f nHz ' + '\n' + r'$\Omega_{\rm{frame}} - \Omega_0$' + ' = %.2f nHz') %(om, om - om0)
        else:
            samplelabel += '\n' + r'$\Omega_{\rm{frame}} = \Omega_0$'

        position_tag = '_clat' + lat_format(clat) + '_rval%.3f' %sampleval
    elif shav:
        axislabel = r'$r/R_\odot$'
        samplelabel = ''
        position_tag = ''

    # Put some useful information on the title
    maintitle = dirname_stripped 
    if not shav:
        maintitle += '\n' + samplelabel
    if 'navg' in clas:
        navg = clas['navg']
        averaging_time = (times[-1] - times[0])/len(times)*navg
        maintitle += '\n' + ('t_avg = %.1f Prot' %averaging_time)
    else: 
        maintitle += '\nt_avg = none'

    if not shav:
        print('plotting sampleval = %0.3f (i = %02i)' %(sampleval, isampleval))
   
    # make plot
    fig, axs, fpar = make_figure(nplots=nplots, ncol=1, sub_width_inches=sub_width_inches, sub_height_inches=sub_height_inches, margin_left_inches=margin_left_inches, margin_right_inches=margin_right_inches, margin_top_inches=margin_top_inches, margin_bottom_inches=margin_bottom_inches)

    for iplot in range(nplots):
        ax = axs[iplot, 0]
        if rad:
            field = terms[iplot][:, isampleval, :]
        else:
            field = terms[iplot][:, :, isampleval]
        plot_timey(field, times, sampleaxis, fig, ax, fontsize=fontsize, **clas)
                
        #  title the plot
        ax.set_title(clas['titles'][iplot], fontsize=fontsize)

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
        if lon:
            plotdir = my_mkdir(clas0['plotdir'] + 'timelon/')
            basename = 'timelon_%08i_%08i' %(iter1, iter2)
            if not om is None:
                basename += '_om%.0f' %om
        else:
            plotdir = my_mkdir(clas0['plotdir'] + datatype + '/')
            basename = datatype + '_%08i_%08i' %(iter1, iter2)
        savename = basename + clas0['tag'] + position_tag + '.png'
        print ("saving", plotdir + savename)
        print ("=======================================")
        plt.savefig(plotdir + savename, dpi=200)

    # Show the plot if only plotting at one latitude
    if clas0['showplot'] and len(isamplevals) == 1:
        plt.show()
    else:
        plt.close()
