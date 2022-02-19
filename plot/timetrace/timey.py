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
kwargs_default = dict({'the_file': None, 'ntot': 2000, 'clat': 10, 'dlat': 0, 'om': None, 'rad': False, 'lon': False, 'shav': False, 'isamplevals': np.array([0]), 'samplevals': None, 'rcut': None, 'groupname': 'b'})
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

datatype = 'timelat'
sampleaxis = di_grid['tt_lat']
if kw.rad:
    datatype = 'timerad'
    sampleaxis = di_grid['rr']/rsun
elif kw.lon:
    if not kw.om is None:
        om0 = 1/time_unit*1e9 # frame rate, nHz
    datatype = 'timelon_clat' + lat_format(clat) + '_dlat%03.0f' %dlat
    sampleaxis = di_grid['lons']
elif kw.shav:
    datatype = 'timeshav'
    sampleaxis = di_grid['rr']/rsun

dataname = datatype
if 'groupname' in kw:
    dataname += '_' + kw.groupname
if not kw.rcut is None:
    dataname += '_rcut%0.3f' %kw.rcut

dataname += clas0['tag']

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
if not kw.shav:
    samplevals_avail = di['samplevals']

# time range
iter1, iter2 = get_iters_from_file(kw.the_file)
times /= time_unit

# Subtract DR, if desired
if kw.lon:
    if not kw.om is None:
        print ("plotting in rotating frame om = %.1f nHz" %kw.om)
        print ("compare this to frame rate    = %.1f nHz" %om0)
        phi_deflections = (times*(kw.om - om0)) % 1 # between zero and one
        nphi = len(sampleaxis)
        for it in range(len(times)):
            phi_deflection = phi_deflections[it]
            nroll = int(phi_deflection*nphi)
            vals[it] = np.roll(vals[it], -nroll, axis=0)

# maybe thin data
if not kw.ntot == 'full':
    print ("ntot = %i" %kw.ntot)
    print ("before thin_data: len(times) = %i" %len(times))
    times = thin_data(times, kw.ntot)
    iters = thin_data(iters, kw.ntot)
    vals = thin_data(vals, kw.ntot)
    print ("after thin_data: len(times) = %i" %len(times))

# these all need to be arrays
kw.qvals = make_array(kw.qvals)
kw.isamplevals = make_array(kw.isamplevals)
kw.samplevals = make_array(kw.samplevals)

# get raw traces of desired variables
terms = []
for qval in kw.qvals:
    qind = np.argmin(np.abs(qvals_avail - qval))
    terms.append(vals[:, :, :, qind])

# set figure dimensions
sub_width_inches = 7.5
sub_height_inches = 2.0
margin_bottom_inches = 1/2 # space for x-axis label
margin_top_inches = 3/4
if kw.lon:
    margin_top_inches =  1 + 1/4
margin_left_inches = 5/8 # space for latitude label
margin_right_inches = 7/8 # space for colorbar
if 'ycut' in clas:
    margin_right_inches *= 2
nplots = len(terms)

# determine desired levels to plot
if not kw.shav: # kw.shav means integrated over latitude, so can't choose
    # specific latitude levels
    if not kw.samplevals is None: # isamplevals being set indirectly
        # check for special 'all' option
        itsall = False
        if len(kw.samplevals) == 1:
            if kw.samplevals[0] == 'all':
                itsall = True
        if itsall:
            kw.isamplevals = np.arange(len(samplevals_avail))
        else:
            kw.isamplevals = np.zeros_like(kw.samplevals, dtype='int')
            for i in range(len(kw.samplevals)):
                kw.isamplevals[i] = np.argmin(np.abs(samplevals_avail - kw.samplevals[i]))

# Loop over the desired levels and save plots
for isampleval in kw.isamplevals:
    if not kw.shav:
        sampleval = samplevals_avail[isampleval]

    # set some labels 
    axislabel = 'latitude (deg)'
    samplelabel =  r'$r/R_\odot$' + ' = %.3f' %sampleval
    position_tag = '_rval%.3f' %sampleval
    if kw.rad:
        axislabel = r'$r/R_\odot$'
        samplelabel = 'lat = ' + lat_format(sampleval)
        position_tag = '_lat' + lat_format(sampleval)
    elif kw.lon:
        axislabel = 'longitude (deg)'
        samplelabel = 'clat = ' + lat_format(clat) + '\n' +  r'$r/R_\odot$' + ' = %.3f' %sampleval
        if not kw.om is None:
            samplelabel += '\n' + (r'$\Omega_{\rm{frame}}$' + ' = %.1f nHz ' + '\n' + r'$\Omega_{\rm{frame}} - \Omega_0$' + ' = %.2f nHz') %(om, om - om0)
        else:
            samplelabel += '\n' + r'$\Omega_{\rm{frame}} = \Omega_0$'

        position_tag = '_clat' + lat_format(clat) + '_rval%.3f' %sampleval
    elif kw.shav:
        axislabel = r'$r/R_\odot$'
        samplelabel = ''
        position_tag = ''

    # Put some useful information on the title
    maintitle = dirname_stripped 
    if not kw.shav:
        maintitle += '\n' + samplelabel
    if kw.navg is None:
        maintitle += '\nt_avg = none'
    else:
        averaging_time = (times[-1] - times[0])/len(times)*kw.navg
        maintitle += '\n' + ('t_avg = %.1f Prot' %averaging_time)

    if not kw.shav:
        print('plotting sampleval = %0.3f (i = %02i)' %(sampleval, isampleval))
   
    # make plot
    fig, axs, fpar = make_figure(nplots=nplots, ncol=1, sub_width_inches=sub_width_inches, sub_height_inches=sub_height_inches, margin_left_inches=margin_left_inches, margin_right_inches=margin_right_inches, margin_top_inches=margin_top_inches, margin_bottom_inches=margin_bottom_inches)

    for iplot in range(nplots):
        ax = axs[iplot, 0]
        if kw.rad:
            field = terms[iplot][:, isampleval, :]
        else:
            field = terms[iplot][:, :, isampleval]
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
        basename = dataname + '_%08i_%08i' %(iter1, iter2)
        plotdir = my_mkdir(clas0['plotdir'] + '/' + datatype)
        if kw.lon and not om is None:
            basename += '_om%.0f' %om
        savename = basename + clas0['tag'] + position_tag + '.png'
        print ("saving", plotdir + '/' + savename)
        plt.savefig(plotdir + '/' + savename, dpi=200)

    # Show the plot if only plotting at one latitude
    if clas0['showplot'] and len(kw.isamplevals) == 1:
        plt.show()
    else:
        plt.close()
    print ("=======================================")
