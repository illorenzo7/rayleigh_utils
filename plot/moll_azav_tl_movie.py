# Created: 03/29/2019
# Author: Loren Matilsky
#
# Purpose: create a string of figures with 3 subplot each: for Mollweide 
# slices, AZ_Avgs, and a time-latitude diagram for a specified interval

# Note that the desired quantity must be available in both Shell_Slices and
# AZ_Avgs; this will generally be one of [1, 2, 3, 501, 502, 801, 802, 803]
# FOR LATER: Change to be compatiable with derivative thermo variables rho, T
import matplotlib as mpl
mpl.use('TkAgg')
from matplotlib import colors
import matplotlib.pyplot as plt
plt.rcParams['mathtext.fontset'] = 'dejavuserif'
csfont = {'fontname':'DejaVu Serif'}
import numpy as np
import sys, os
sys.path.append(os.environ['raco'])
sys.path.append(os.environ['rapp'])
from common import *
        get_dict, get_widest_range_file, sci_format, saturate_array,\
        get_symlog_params, get_exp, append_logfile, allthrees_start,\
        sci_format
from sslice_util import plot_moll, get_satvals
from get_sslice import get_sslice
from azav_util import plot_azav
from rayleigh_diagnostics import Shell_Slices, AZ_Avgs
from get_parameter import get_parameter
from varprops import texlabels, var_indices, texunits

# Get command line arguments, mainly the directory on which to perform the
# plotting

# Analysis directory
dirname = sys.argv[1]

# Rayleigh data dirs
slicedatadir = dirname + '/Shell_Slices/'
azdatadir = dirname + '/AZ_Avgs/'
file_list, int_file_list, nfiles = get_file_lists(slicedatadir)
# Depending on when/how the simulation ended, data may differ in output number 
# by at most 1
# Thus ...
file_list = file_list[:-1]; 
int_file_list = int_file_list[:-1]
nfiles -= 1

# Rotation rate and period
Omega0 = get_parameter(dirname, 'angular_velocity')
Prot = 2*np.pi/Omega0

# Other defaults
minmax = None
linscale = None # may be used for symlog stuff
linthresh = None
symlog = False
restart = False
log_progress = False # user can demand using a logfile via -logfile
varname = 'bp' # by default plot the azimuthal field
ir = 0 # by default plot just below the surface
rval = None # can also find ir by finding the closest desired radial value
            # (normalized by Rsun)
clon = 0.
count = 0 # start the plot count at 0000 by default
tlabel_string_width = None
alpha = 0.05 # how greyed out some of the TL is
rbcz = None

# Get desired data range to plot
args = sys.argv[2:]
nargs = len(args)
default_range = True
for i in range(nargs):
    arg = args[i]
    if arg in range_options:
        default_range = False

if default_range:
    index_first, index_last = 0, nfiles - 1  
    # By default plot all the shell slices
else:
    index_first, index_last = get_desired_range(int_file_list, args)

# Get the file names (zero-filled iteration strings) associated with desired
# plotting range
fnames = file_list[index_first:index_last+1] 

# Change other defaults
for i in range(nargs):
    arg = args[i]
    if arg == '-minmax':
        minmax = float(args[i+1]), float(args[i+2]),\
                float(args[i+3]), float(args[i+4]),\
                float(args[i+5]), float(args[i+6])
    elif arg == '-rbcz':
        rbcz = float(args[i+1])
    elif arg == '-symlog':
        symlog = True
    elif arg == '-linscale':
        linscale = float(args[i+1]), float(args[i+2]), float(args[i+3])    
    elif arg == '-linthresh':
        linthresh = float(args[i+1]), float(args[i+2]), float(args[i+3])    
    elif arg == '-var':
        varname = args[i+1]
    elif arg == '-ir':
        ir = int(args[i+1])
    elif arg == '-rval':
        rval = float(args[i+1])
    elif arg == '-clon':
        clon = float(args[i+1])
    elif arg == '-start':
        count = int(args[i+1])
    elif arg == '-restart':
        restart = True # attempt to complete a partially filled
                        # plot directory
    elif arg == '-tlabel':
        tlabel_string_width = int(args[i+1])
    elif arg == '-alpha':
        alpha = float(args[i+1])
    elif arg == '-logfile':
        log_progress = True

# varlabel and var_index will depend on the variable being plotted
varlabel = texlabels[varname]
var_index = var_indices[varname]

# Get some useful info from the first sslice
a0 = Shell_Slices(slicedatadir + fnames[0], '')

# Figure out ir if user specified rval
if not rval is None:
    ir = np.argmin(np.abs(a0.radius/rsun - rval))
# Replace desired rval (or "None") by the actual rval
rval = a0.radius[ir]/rsun

# Start time
time0 = a0.time[0]

# Create the save directory if it doesn't already exist
plotdir = dirname + '/plots/moll_azav_tl_movies/' + varname + '/' +\
        ('rval%0.3f' %rval) + '/'
if not os.path.isdir(plotdir):
    os.makedirs(plotdir)

# Name of logfile to catalog progress, if desired
logfile = dirname + '/plots/zz_logfile_mollazavtlmovie_' + varname +\
    ('_rval%0.3f' %rval) + '.txt'
if log_progress:
    append_logfile(logfile, '==================================\n')
    append_logfile(logfile, '========== Begin Log File ==========\n')
    
# Get grid information from last AZ_Avgs file
az0 = AZ_Avgs(azdatadir + file_list[index_last], '')
rr = az0.radius
ri, ro = np.min(rr), np.max(rr)
d = ro - ri
rr_depth = (ro - rr)/d
rr_height = (rr - ri)/d
sint = az0.sintheta
cost = az0.costheta
tt = np.arccos(cost)
tt_lat = (np.pi/2 - tt)*180/np.pi
nr = az0.nr
nt = az0.ntheta

# compute some derivative quantities for the grid
tt_2d, rr_2d = np.meshgrid(tt, rr, indexing='ij')
sint_2d = np.sin(tt_2d); cost_2d = np.cos(tt_2d)

# Get max time, in rotations
max_time = az0.time[-1]/Prot
max_time_exp = int(np.floor(np.log10(max_time)))
# Set the tlabel width unless the user wants to do so manually
if tlabel_string_width is None:
    tlabel_string_width = max_time_exp + 3 # Make room for full whole number
        # width, a decimal point, and one tenths place after decimal

# Read in the time-latitude data (dictionary form)
datadir = dirname + '/data/'
time_latitude_file = get_widest_range_file(datadir, 'time-latitude')
print ('Getting time-latitude trace from ' + datadir +\
       time_latitude_file + ' ...')
if log_progress:
    append_logfile (logfile, 'Getting time-latitude trace from ' +\
            datadir + time_latitude_file + ' ...\n')
di = get_dict(datadir + time_latitude_file)

vals = di['vals']
times = di['times']/Prot - allthrees_start
rinds = di['rinds'] # radial locations sampled for the trace
rvals_sampled = rr[rinds]/rsun
ir_tl = np.argmin(np.abs(rvals_sampled - rval))
qvals = np.array(di['qvals'])
iq_tl = np.argmin(np.abs(qvals - var_index))

# Same values used in each time-latitude plot
tl_vals = vals[:, :, ir_tl, iq_tl]

# Get saturation values for Shell_Slice and AZ_Avgs of variable in question
# possibly also get linscale/linthresh if -symlog was specified
field_slice = get_sslice(a0, varname, dirname=dirname)[:, :, ir]
field_az = az0.vals[:, :, az0.lut[var_index], 0]
if minmax is None: # Set saturation values by the field from
    # sslice max/min from first slice
    min_slice, max_slice = get_satvals(field_slice, symlog=symlog)
    min_az, max_az = get_satvals(field_az, symlog=symlog)
    min_tl, max_tl = get_satvals(tl_vals, symlog=symlog)
else:
    min_slice, max_slice = minmax[0], minmax[1]
    min_az, max_az = minmax[2], minmax[3]
    min_tl, max_tl = minmax[4], minmax[5]
    
# May need linthresh/linscale for each plot if -symlog was True
    
# get the default values first
linthresh_slice, linscale_slice =\
    get_symlog_params(field_slice, field_max=max_slice)
linthresh_az, linscale_az =\
    get_symlog_params(field_az, field_max=max_az)
linthresh_tl, linscale_tl =\
    get_symlog_params(tl_vals, field_max=max_tl)
if not linthresh is None: # ... then possibly overwrite them
    linthresh_slice, linthresh_az, linthresh_tl = linthresh
if not linscale is None:
    linscale_slice, linscale_az, linscale_tl = linscale  

# Factor out the exponent in tl_vals (if symlog is not True)    
# This was in the for-loop before, which messed up the exponent (first plot
# would have, e.g., 10^3, all following would have 10^0))
if not symlog: # don't yet need to worry about logscale = True, since I have
    # only been using moll_azav_tl_movie for the raw (signed) fluid variables
    maxabs_tl = max(abs(min_tl), abs(max_tl))
    tl_exp = get_exp(maxabs_tl)
    tl_vals /= 10**tl_exp
    min_tl /= 10**tl_exp
    max_tl /= 10**tl_exp

# Saturate the array at the extremal values (necessary when using
# contourf for the color plot)
saturate_array(tl_vals, min_tl, max_tl)

# Must set the norm if using symlog
# default norm_tl = None
norm_tl = None
if symlog:
    norm_tl = colors.SymLogNorm(linthresh=linthresh_tl,\
            linscale=linscale_tl, vmin=min_tl, vmax=max_tl)
    nlevs_per_interval = 100
    log_thresh = np.log10(linthresh_tl)
    log_max = np.log10(max_tl)
    levels_neg = -np.logspace(log_max, log_thresh, nlevs_per_interval,\
            endpoint=False)
    levels_mid = np.linspace(-linthresh_tl, linthresh_tl,\
            nlevs_per_interval, endpoint=False)
    levels_pos = np.logspace(log_thresh, log_max, nlevs_per_interval)
    levels_tl = np.hstack((levels_neg, levels_mid, levels_pos))
else:
    levels_tl = np.linspace(min_tl, max_tl, 150)

# General parameters for figure and axes
moll_height_inches = 3.9
moll_width_inches = 2*moll_height_inches
azav_height_inches = moll_height_inches
azav_width_inches = 0.5*azav_height_inches
tl_height_inches = 1.5
margin_inches = 1/8
margin_top_inches = 1/2
margin_bottom_inches = 0.75*(2 - (rbcz is None)) 
    # larger bottom margin to make room for colorbar(s)
hspace_inches = 3/4

# Global figure dimensions
fig_height_inches = margin_bottom_inches + tl_height_inches +\
        hspace_inches + moll_height_inches + margin_top_inches
fig_width_inches = margin_inches + moll_width_inches +\
        margin_inches + azav_width_inches + margin_inches
fig_aspect = fig_height_inches/fig_width_inches

# "Non-dimensional" figure parameters
margin_x = margin_inches/fig_width_inches
margin_y = margin_inches/fig_height_inches
hspace = hspace_inches/fig_height_inches
margin_top = margin_top_inches/fig_height_inches
margin_bottom = margin_bottom_inches/fig_height_inches

moll_width = moll_width_inches/fig_width_inches
moll_height = moll_height_inches/fig_height_inches
azav_width = azav_width_inches/fig_width_inches
azav_height = azav_height_inches/fig_height_inches

tl_height = 1. - margin_top - moll_height - hspace - margin_bottom
tl_width = 0.5*moll_width
#tl_center = margin_x + moll_width + 0.5*margin_x
tl_center = margin_x + 0.5*(1 - 2*margin_x)
tl_left = tl_center - 0.5*tl_width

# If trying to restart, see where the last run got to:
if restart:
    already_plotted = os.listdir(plotdir)
    n_plotted = len(already_plotted)
    # Make a new list of the remaining plots to make;
    # Make the last plot again, since it may be f**ked up
    count = n_plotted - 1
    fnames = fnames[n_plotted - 1:]
    print ("------------------------")
    print ("Restart desired:")
    print ("Remaking img%04i.png" %count)
    print ("Will plot %s through %s" %(fnames[0], fnames[-1]))
    print ("as img%04i.png through img%04i.png"\
            %(count, count + len(fnames) - 1))
    print ("------------------------")
else: # or else announce we are plotting everything
    print ("Plotting slice/AZ_Avg " + fnames[0] + " through " + fnames[-1])
    print ("Or img%04i.png through img%04i.png"\
            %((count, count + len(fnames) - 1)))

# Now do the actual plotting
for fname in fnames:
    # Read in desired shell slice
    a = Shell_Slices(slicedatadir + fname, '')
    az = AZ_Avgs(azdatadir + fname, '')

    # Loop over times and make plots
    for j in range(a.niter):
        time = a.time[j]/Prot - allthrees_start
            
        # Create the plot using subplot axes
        savename = 'img' + str(count).zfill(4) + '.png'
        print('Plotting moll_azav_tl: ' + varname +\
                (', rval = %0.3f, ' %rval) + 'iter ' + fname + ' ...')
        if log_progress:
            append_logfile(logfile, 'Plotting moll_azav_tl: ' + varname +\
                (', rval = %0.3f, ' %rval) + 'iter ' + fname + ' ...\n')        
        fig = plt.figure(figsize=(fig_width_inches, fig_height_inches))
        ax_moll = fig.add_axes([margin_x,\
                margin_bottom + tl_height + hspace,\
                moll_width, moll_height])
        ax_azav = fig.add_axes([margin_x + moll_width + margin_x,\
                margin_bottom + tl_height + hspace,\
                azav_width, azav_height])
        ax_tl = fig.add_axes([tl_left, margin_bottom, tl_width, tl_height])

        # Colorbar axes to accompany the time-latitude plot (the other 
        # colorbars come pre-packaged with plot_moll and plot_azav)
        tl_cbar_left = tl_left + tl_width + 2*margin_x 
        tl_cbar_bottom = margin_bottom
        tl_cbar_height = 0.9*tl_height
        tl_cbar_width = 1/20*tl_cbar_height
        cax_tl = fig.add_axes([tl_cbar_left, tl_cbar_bottom,\
                tl_cbar_width, tl_cbar_height])
       
        # Make the Mollweide plot
        vals = get_sslice(a, varname, dirname=dirname)
        field = vals[:, :, ir]
        plot_moll(field, cost, fig=fig, ax=ax_moll, varname=varname,\
                minmax=(min_slice, max_slice), clon=clon, symlog=symlog,\
                linscale=linscale_slice, linthresh=linthresh_slice) 
        time_string = '%.1f' %time
        time_string = time_string.zfill(tlabel_string_width)
        title = varlabel + '     ' +\
                (r'$r/R_\odot\ =\ %0.3f$' %rval) +\
                '     ' + (r'$t = %s\ P_{\rm{rot}}$' %time_string)
        fig.text(margin_x + 0.5*moll_width, 1. - 0.5*margin_top,\
                title, ha='center', va='center', **csfont, fontsize=14)

        # Make the AZ_Avgs plot:
        field = az.vals[:, :, az.lut[var_index], 0]
        plot_azav (field, rr, cost, fig=fig, ax=ax_azav,\
               units=texunits[varname], minmax=(min_az, max_az),\
               plotcontours=False, plotlatlines=True, cbar_fs=10,\
               rvals=[rval], symlog=symlog, linthresh=linthresh_az,\
               linscale=linscale_az)

        # Make the time-latitude plot underneath everything

        # Arrays for plotting time-latitude diagram (use two, with the one
        # after the current time greyed out)
        #times2, lats2 = np.meshgrid(times, tt_lat, indexing='ij')
        it_current = np.argmin(np.abs(times - time))
        times_1 = times[:it_current+2]
            # Have some overlap to prevent getting
            # empty arrays at the first and last times
        times_2 = times[it_current:]
        times_2d_1, lats_2d_1 = np.meshgrid(times_1, tt_lat, indexing='ij')
        times_2d_2, lats_2d_2 = np.meshgrid(times_2, tt_lat, indexing='ij')
        tl_vals_1 = tl_vals[:it_current+2]
        tl_vals_2 = tl_vals[it_current:]
            
        im_tl = ax_tl.contourf(times_2d_1, lats_2d_1, tl_vals_1,\
                cmap='RdYlBu_r', vmin=min_tl, vmax=max_tl, norm=norm_tl,\
                levels=levels_tl) # set colorbar by un-greyed-out image
        ax_tl.contourf(times_2d_2, lats_2d_2, tl_vals_2,\
                cmap='RdYlBu_r', vmin=min_tl, vmax=max_tl, alpha=alpha,\
                norm=norm_tl, levels=levels_tl) 
                # ...grey out uncovered time with alpha

        # Set up the time-latitude colorbar
        cbar = plt.colorbar(im_tl, cax=cax_tl)    
            
        # Set custom ticks, dependent on True/False-ness of symlog
        if symlog:
            title = r'$\rm{G}$'
            nlin = 5
            nlog = 6
            lin_ticks = np.linspace(-linthresh_tl, linthresh_tl, nlin)
            log_ticks1 = np.linspace(min_tl, -linthresh_tl, nlog,\
                    endpoint=False)
            log_ticks2 = -log_ticks1[::-1]
            ticks = np.hstack((log_ticks1, lin_ticks, log_ticks2))
            nticks = nlin + 2*nlog
            cbar.set_ticks(ticks)
            ticklabels = []
            for i in range(nticks):
                ticklabels.append(r'')
            ticklabels[0] = sci_format(min_tl)
            ticklabels[nlog] = sci_format(-linthresh_tl)
            ticklabels[nticks//2] = r'$0$'
            ticklabels[nlog + nlin - 1] = sci_format(linthresh_tl)
            ticklabels[nticks - 1] = sci_format(max_tl)
            cbar.set_ticklabels(ticklabels)
        else:
            title = r'$\times10^{%i}\ \rm{G}$' %tl_exp
            cbar.set_ticks([min_tl, 0, max_tl])
            cbar.set_ticklabels(['%1.1f' %min_tl, '0', '%1.1f' %max_tl])

        cax_tl.set_title(title, **csfont,\
        fontsize=10)

        # Label the x-axis (time)
        xlabel = 'time (' + r'$P_{\rm{rot}}$' + ')'
        ax_tl.set_xlabel(xlabel, **csfont)

        # Get ticks everywhere for TL plot
        plt.sca(ax_tl)
        plt.minorticks_on()
        plt.tick_params(top=True, bottom=True,\
                left=True, right=True, direction='in', which='both')

        # Label y-axis (latitude in degrees)
        ax_tl.set_ylabel('latitude (deg.)', **csfont)
        ax_tl.set_yticks(np.arange(-90, 90, 30))
        ax_tl.set_ylim(-90, 90)
        
        # Save the plot
        print('Saving plot: ' + plotdir + savename + ' ...')
        if log_progress:
            append_logfile(logfile, 'Saving plot: ' + plotdir + savename +\
                       ' ...\n')
        plt.savefig(plotdir + savename, dpi=150)
        plt.close()
        count += 1
if log_progress:
    append_logfile(logfile, '========== End Log File ==========\n')
    append_logfile(logfile, '==================================\n')
