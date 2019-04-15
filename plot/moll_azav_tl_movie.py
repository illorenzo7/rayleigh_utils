# Created: 03/29/2019
# Author: Loren Matilsky
#
# Purpose: create a string Mollweide slices for every shell slice in a
# directory, to blend together in an ffmpeg movie
import matplotlib as mpl
mpl.use('TkAgg')
import matplotlib.pyplot as plt
plt.rcParams['mathtext.fontset'] = 'dejavuserif'
csfont = {'fontname':'DejaVu Serif'}
from binormalized_cbar import MidpointNormalize
import numpy as np
import sys, os
sys.path.append(os.environ['co'])
sys.path.append(os.environ['rapp'])
from common import get_file_lists, strip_dirname, rsun, get_desired_range,\
        get_dict, get_widest_range_file
from sslice_util import plot_moll, get_satvals, get_sslice
from azavg_util import plot_azav
from rayleigh_diagnostics import Shell_Slices, AZ_Avgs
from get_parameter import get_parameter
from varprops import texlabels, var_indices, texunits

# Get command line arguments
dirname = sys.argv[1]
slicedatadir = dirname + '/Shell_Slices/'
azdatadir = dirname + '/AZ_Avgs/'

file_list, int_file_list, nfiles = get_file_lists(slicedatadir)
# Different data times may differ in output number by at most 1
# Thus ...
file_list = file_list[:-1]; 
int_file_list = int_file_list[:-1]
nfiles -= 1

Omega0 = get_parameter(dirname, 'angular_velocity')
Prot = 2*np.pi/Omega0

minmax = None
restart = False
varname = 'bp' # by default plot the zonal field
ir = 0 # by default plot just below the surface
rval = None # can also find ir by finding the closest point
clon = 0
count = 0 # start the plot count at 0000 by default

args = sys.argv[2:]
nargs = len(args)
try:
    index_first, index_last = get_desired_range(int_file_list, args)
except:
    index_first, index_last = 0, nfiles - 1  
    # By default plot all the shell slices

fnames = file_list[index_first:index_last+1] 

for i in range(nargs):
    arg = args[i]
    if arg == '-minmax':
        minmax = float(args[i+1]), float(args[i+2]),\
                float(args[i+3]), float(args[i+4]),\
                float(args[i+5]), float(args[i+6])
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

posdef = False
if 'sq' in varname:
    posdef = True

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
rval = a0.radius[ir]/rsun
plotdir = dirname + '/plots/moll_azav_tl_movies/' + varname + '/' +\
        ('rval%0.3f' %rval) + '/'
if not os.path.isdir(plotdir):
    os.makedirs(plotdir)

# Get grid information from first AZ_Avgs file
az0 = AZ_Avgs(azdatadir + file_list[index_first], '')
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
xx = rr_2d*sint_2d
zz = rr_2d*cost_2d

# Get saturation values for Shell_Slice and AZ_Avgs of variable in question
if minmax is None: # Set saturation values by the field from
    # Sslice max/min from first slice
    field_slice = get_sslice(a0, varname, dirname=dirname)[:, :, ir]
    min_slice, max_slice = -3*np.std(field_slice), 3.*np.std(field_slice)
    if posdef:
        min_slice, max_slice = get_satvals(field_slice, posdef=True)

    field_az = az0.vals[:, :, az0.lut[var_index], 0]
    nstd = 5
    min_az, max_az = -nstd*np.std(field_az), nstd*np.std(field_az)
else:
    min_slice, max_slice = minmax[0], minmax[1]
    min_az, max_az = minmax[2], minmax[3]

# Read in the time-latitude data (dictionary form)
datadir = dirname + '/data/'
time_latitude_file = get_widest_range_file(datadir, 'time-latitude')
print ('Getting time-latitude trace from ' + datadir +\
       time_latitude_file + ' ...')
di = get_dict(datadir + time_latitude_file)

vals = di['vals']
times = di['times']/Prot
rinds = di['rinds'] # radial locations sampled for the trace
rvals_sampled = rr[rinds]/rsun
ir_tl = np.argmin(np.abs(rvals_sampled - rval))
qvals = np.array(di['qvals'])
iq_tl = np.argmin(np.abs(qvals - var_index))

# Arrays for plotting (every time!)
times2, lats2 = np.meshgrid(times, tt_lat, indexing='ij')
tl_vals = vals[:, :, ir_tl, iq_tl]
if minmax is None:
    min_tl, max_tl = -3.*np.std(tl_vals), 3.*np.std(tl_vals)
else:
    min_tl, max_tl = minmax[4], minmax[5]

# General parameters for main axis/color bar
moll_height_inches = 4.
moll_width_inches = 2*moll_height_inches
azav_height_inches = moll_height_inches
azav_width_inches = 0.5*azav_height_inches
tl_height_inches = 1.5
margin_inches = 1/8
margin_top_inches = 1/2
margin_bottom_inches = 1/2
margin_left_tl_inches = 3/4
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
margin_left_tl = margin_left_tl_inches/fig_width_inches

moll_width = moll_width_inches/fig_width_inches
moll_height = moll_height_inches/fig_height_inches
azav_width = azav_width_inches/fig_width_inches
azav_height = azav_height_inches/fig_height_inches

tl_height = 1. - margin_top - moll_height - hspace - margin_bottom
tl_width = 0.5*moll_width
#tl_width = 1. - margin_left_tl - margin_x

print ("Plotting slice/AZ_Avg " + fnames[0] + " through " + fnames[-1])
print ("Or img%04i.png through img%04i.png"\
        %((count, count + len(fnames) - 1)))

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

firstplot = True
for fname in fnames:
    # Read in desired shell slice
    if firstplot and not restart:
        a = a0
        az = az0
        firstplot = False
    else:
        a = Shell_Slices(slicedatadir + fname, '')
        az = AZ_Avgs(azdatadir + fname, '')

    for j in range(a.niter):
        time = a.time[j]
            
        # Create the plot using subplot axes

        # Loop over times and make plots
        savename = 'img' + str(count).zfill(4) + '.png'
        print('Plotting moll_azav_tl: ' + varname +\
                (', rval = %0.3f, ' %rval) + 'iter ' + fname + ' ...')
        print('Saving plot: ' + plotdir + savename + ' ...')
        fig = plt.figure(figsize=(fig_width_inches, fig_height_inches))
        ax_moll = fig.add_axes([margin_x,\
                margin_bottom + tl_height + hspace,\
                moll_width, moll_height])
        ax_azav = fig.add_axes([margin_x + moll_width + margin_x,\
                margin_bottom + tl_height + hspace,\
                azav_width, azav_height])
        ax_tl = fig.add_axes([margin_left_tl, margin_bottom,\
                tl_width, tl_height])

        # Colorbar to accompany the time-latitude plot
        tl_cbar_left = margin_left_tl + tl_width + 2*margin_x 
        tl_cbar_bottom = margin_bottom
        tl_cbar_height = 0.9*tl_height
        tl_cbar_width = 1/20*tl_cbar_height
        cax_tl = fig.add_axes([tl_cbar_left, tl_cbar_bottom,\
                tl_cbar_width, tl_cbar_height])
       
        # Make the Mollweide plot
        plot_moll(fig, ax_moll, a, dirname, varname, ir=ir,\
                minmax=(min_slice, max_slice), clon=clon, plot_title=False) 
        title = varlabel + '     ' +\
                (r'$r/R_\odot\ =\ %0.3f$' %rval) +\
                '     ' + (r'$t = %02.1f\ P_{\rm{rot}}$' %(time/Prot))
        fig.text(margin_x + 0.5*moll_width, 1. - 0.5*margin_top,\
                title, ha='center', va='center', **csfont, fontsize=14)
#        ax_moll.axis('off')

        # Make the AZ_Avgs plot:
        var_az = az.vals[:, :, az.lut[var_index], 0]
        plot_azav (fig, ax_azav, var_az, rr, cost, sint,\
               units = texunits[varname], boundstype = 'manual',\
               caller_minmax = (min_az, max_az),\
               norm=MidpointNormalize(0), plotcontours=False,\
               plotlatlines=True, fsize=10)

        # Mark the line of the rvalue we're plotting
        r_over_ro = rval/(ro/rsun)
        linex = r_over_ro*sint
        linez = r_over_ro*cost
        ax_azav.plot(linex, linez, 'k--', linewidth=0.5)


        # Make the time-latitude plot underneath everything
        # Factor out the exponent
        maxabs_tl = max(abs(min_tl), abs(max_tl))
        tl_exp = int(np.floor(np.log10(maxabs_tl)))
        tl_vals /= 10**tl_exp
        min_tl /= 10**tl_exp
        max_tl /= 10**tl_exp
        im_tl = ax_tl.pcolormesh(times2, lats2, tl_vals,\
                cmap='RdYlBu_r', vmin=min_tl, vmax=max_tl)

        # Set up the time-latitude colorbar
        cax_tl.set_title(r'$\times10^{%i}\ \rm{G}$' %tl_exp, **csfont,\
                fontsize=10)
        cbar = plt.colorbar(im_tl, cax=cax_tl)
        cbar.set_ticks([min_tl, 0, max_tl])
        cbar.set_ticklabels(['%1.1f' %min_tl, '0', '%1.1f' %max_tl])

        # Put a vertical line at current time
        linex = np.zeros(100) + time/Prot
        liney = np.linspace(-90, 90, 100)
        ax_tl.plot(linex, liney, 'k--', linewidth=0.8)

        # Label the x-axis (time)
        xlabel = 'time (' + r'$P_{\rm{rot}}$' + ')'
        ax_tl.set_xlabel(xlabel, **csfont)

        # Get ticks everywhere
        plt.sca(ax_tl)
        plt.minorticks_on()
        plt.tick_params(top=True, bottom=True,\
                left=True, right=True, direction='in', which='both')

        # Label y-axis (radius in units of rsun)
        ax_tl.set_ylabel('latitude (deg.)', **csfont)
        ax_tl.set_yticks(np.arange(-90, 90, 30))
        ax_tl.set_ylim(-90, 90)

        plt.savefig(plotdir + savename, dpi=200)
        count += 1
        plt.close()
