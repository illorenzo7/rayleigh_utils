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
kwargs_default = dict({'the_file': None, 'ntot': 500, 'clat': 10, 'dlat': 0, 'om': None, 'irvals': np.array([0]), 'rvals': None, 'qvals': np.array([1])})
kwargs_default.update(plot_timey_kwargs_default)

# check for bad keys
find_bad_keys(kwargs_default, clas, clas0['routinename'], justwarn=True)

# overwrite defaults
kw = update_dict(kwargs_default, clas)
kw_plot_timey = update_dict(plot_timey_kwargs_default, clas)

# get the rvals we want
radlevs = get_slice_levels(dirname)
irvals = kw.irvals

rvals = kw.rvals
irvals = kw.irvals
if not kw.rvals is None: # irvals haven't been set directly
    if np.all(kw.rvals == 'all'):
        irvals = np.arange(radlevs.nr)
    else:
        irvals = np.zeros_like(kw.rvals, dtype='int')
        for i in range(len(kw.rvals)):
            irvals[i] = np.argmin(np.abs(radlevs.radius/rsun - kw.rvals[i]))

# and the qvals
qvals = make_array(kw.qvals)

# everything must be array
irvals = make_array(irvals)

# baseline time unit
time_unit, time_label, rotation, simple_label = get_time_unit(dirname)

# get grid info
di_grid = get_grid_info(dirname)
nphi = di_grid['nphi']
lons = di_grid['lons']

om0 = 1/time_unit*1e9 # frame rate, nHz

# set figure dimensions
sub_width_inches = 3.0
sub_height_inches = 9.0
sub_margin_bottom_inches = 1/2 # space for x-axis and label
margin_top_inches =  1 + 1/4
if not kw.om is None:
    margin_top_inches += 1/4
sub_margin_left_inches = 3/4 # space for time label
sub_margin_right_inches = 7/8 # space for colorbar

# loop over data and make plots
for irval in irvals:
    rval = radlevs.radius[irval]/rsun
    for qval in qvals:
        dataname = 'timelon_clat' + lat_format(kw.clat) + '_dlat%03.0f' %kw.dlat + ('_qval%04i_irval%02i' %(qval, irval)) + clas0['tag']

        # get data
        if kw.the_file is None:
            kw.the_file = get_widest_range_file(clas0['datadir'] +\
                'timelon/', dataname)

        # Read in the data
        print ('reading ' + kw.the_file)
        di = get_dict(kw.the_file)
        vals = di['vals']
        times = di['times']
        iters = di['iters']

        # time range
        iter1, iter2 = get_iters_from_file(kw.the_file)

        # Subtract DR, if desired
        if not kw.om is None:
            #print ("plotting in rotating frame om = %.1f nHz" %kw.om)
            #print ("compare this to frame rate    = %.1f nHz" %om0)
            rate_wrt_frame = (kw.om - om0)/1e9
            #phi_deflections = ((times - times[0])*rate_wrt_frame) % 1 
            t0 = times[0]
            # between zero and one
            for it in range(len(times)):
                phi_deflection = ((times[it] - t0)*rate_wrt_frame) % 1
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

        # set some labels 
        samplelabel = 'clat = ' + lat_format(kw.clat) + '\n' +  r'$r/R_\odot$' + ' = %.3f' %rval
        if not kw.om is None:
            samplelabel += '\n' + (r'$\Omega_{\rm{frame}}$' + ' = %.1f nHz ' + '\n' + r'$\Omega_{\rm{frame}} - \Omega_0$' + ' = %.2f nHz') %(kw.om, kw.om - om0)
        else:
            samplelabel += '\n' + r'$\Omega_{\rm{frame}} = \Omega_0 =$' +  '%.2f nHz' %om0

        # Put some useful information on the title
        maintitle = dirname_stripped + '\nqval = %i' %qval
        if not kw.shav:
            maintitle += '\n' + samplelabel

        # Display at terminal what we are plotting
        savename = dataname + '-' + str(iter1).zfill(8) + '_' + str(iter2).zfill(8) + '.png'
  
        # make plot
        fig, axs, fpar = make_figure(sub_width_inches=sub_width_inches, sub_height_inches=sub_height_inches, sub_margin_left_inches=sub_margin_left_inches, sub_margin_right_inches=sub_margin_right_inches, margin_top_inches=margin_top_inches, sub_margin_bottom_inches=sub_margin_bottom_inches)
        ax = axs[0, 0]

        # plot the colormesh
        plot_timey(vals.T, lons, times/time_unit, fig, ax, **kw_plot_timey)
        # make time go downward
        ax.invert_yaxis()

        # title plot
        fig.text(fpar['sub_margin_left'] + fpar['margin_left'], 1 - fpar['margin_top'], maintitle, fontsize=fontsize, ha='left', va='bottom')
        #ax.set_title(maintitle, fontsize=fontsize, ha='left', va='bottom')

        # Put lon label on bottom
        ax.set_xlabel('longitude (deg)', fontsize=fontsize)

        # put time label on side
        ax.set_ylabel('time (' + time_label + ')', fontsize=fontsize)

        # Save the plot
        if clas0['saveplot']:
            # Make appropriate file name to save

            # save the figure
            basename = dataname + '_%08i_%08i' %(iter1, iter2)
            if not kw.om is None:
                plotdir = my_mkdir(clas0['plotdir'] + '/timelon_om%.1f' %kw.om)
            else:
                plotdir = my_mkdir(clas0['plotdir'] + '/timelon')

        if clas0['saveplot']:
            print ("saving", plotdir + '/' + savename)
            plt.savefig(plotdir + '/' + savename, dpi=200)

        # Show the plot if only plotting at one latitude
        if clas0['showplot'] and len(irvals) == 1:
            plt.show()
        else:
            plt.close()
        print ("=======================================")
