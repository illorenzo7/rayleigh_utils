# Author: Loren Matilsky
# Created: 07/02/2020
# This script plots the DR and angular momentum in the meridional plane
# for many different times (to make up the frames of a movie)
# Saves plots in subdirectory
# plots/diffrot_amom_times_sample/
import numpy as np
import pickle
import matplotlib as mpl
mpl.use('TkAgg')
import matplotlib.pyplot as plt
plt.rcParams['mathtext.fontset'] = 'dejavuserif'
csfont = {'fontname':'DejaVu Serif'}
import sys, os
sys.path.append(os.environ['rapp'])
sys.path.append(os.environ['raco'])
sys.path.append(os.environ['rapl'])
from rayleigh_diagnostics import AZ_Avgs, GridInfo
from azav_util import plot_azav, streamfunction
from common import get_widest_range_file, strip_dirname, get_dict,\
        get_file_lists, get_desired_range, sci_format

from get_parameter import get_parameter
from get_eq import get_eq
from time_scales import compute_Prot, compute_tdt
from translate_times import translate_times

# Get directory name and stripped_dirname for plotting purposes
dirname = sys.argv[1]
dirname_stripped = strip_dirname(dirname)

# Directory with AZ_Avgs
radatadir = dirname + '/AZ_Avgs/'

# Get all the file names in datadir and their integer counterparts
file_list, int_file_list, nfiles = get_file_lists(radatadir)

# Directory with data and plots, make the plotting directory if it doesn't
# already exist    
datadir = dirname + '/data/'
plotdir = dirname + '/plots/dr_amom_mc_times_sample/'
if not os.path.isdir(plotdir):
    os.makedirs(plotdir)

# Read command-line arguments (CLAs)
plotcontours = True
plotlatlines = False
minmaxdr = None
minmaxamom = None
minmaxmc = None
AZ_Avgs_file = get_widest_range_file(datadir, 'AZ_Avgs')
rvals = []
navg = 1 # by default average over 1 AZ_Avgs instance (no average)
# for navg > 1, a "sliding average" will be used.
nlevs = 20
plotboundary = True
nskip = 1 # by default don't skip anything
ntot = None 

args = sys.argv[2:]
nargs = len(args)

the_tuple = get_desired_range(int_file_list, args)
if the_tuple is None: # By default average over the first 50 files
    index_first, index_last = 0, 49
else:
    index_first, index_last = the_tuple

# Change defaults
for i in range(nargs):
    arg = args[i]
    if arg == '-minmaxdr' or arg == '-minmax':
        minmaxdr = float(args[i+1]), float(args[i+2])
    elif arg == '-minmaxamom':
        minmaxamom = float(args[i+1]), float(args[i+2])
    elif arg == '-minmaxmc':
        minmaxmc = float(args[i+1]), float(args[i+2])
    elif arg == '-nocontour':
        plotcontours = False
    elif arg == '-rvals':
        rvals_str = args[i+1].split()
        rvals = []
        for rval_str in rvals_str:
            rvals.append(float(rval_str))
    elif arg == '-lats':
        plotlatlines = True
    elif arg == '-nlevs':
        nlevs = int(args[i+1])
    elif arg == '-nobound':
        plotboundary = False
    elif arg == '-navg':
        navg = int(args[i+1])
        if navg % 2 == 0:
            print ("Please don't enter even values for navg!")
            print ("Replacing navg = %i with navg = %i" %(navg, navg + 1))
            navg += 1
    elif arg == '-nskip':
        nskip = int(args[i+1])
    elif arg == '-ntot':
        ntot = int(args[i+1])

# Check to make sure index_last didn't fall beyond the last possible index ...
if index_last > nfiles - navg:
    index_last = nfiles - navg
# compute number files in range
n_analyze = index_last - index_first + 1
# then set nskip if user specified ntot
if not ntot is None:
    nskip = n_analyze//ntot

# Get the baseline time unit
rotation = get_parameter(dirname, 'rotation')
if rotation:
    time_unit = compute_Prot(dirname)
    time_label = r'$\rm{P_{rot}}$'
else:
    time_unit = compute_tdt(dirname)
    time_label = r'$\rm{TDT}$'

# Read in first AZ_Avgs file
az0 = AZ_Avgs(radatadir + file_list[index_first], '')

# Get necessary grid info
rr = az0.radius
ri, ro = np.min(rr), np.max(rr)
cost = az0.costheta
sint = az0.sintheta
nr = az0.nr
nt = az0.ntheta
nq = az0.nq
xx = rr.reshape((1, nr))*sint.reshape((nt, 1))
tt_lat = (np.pi/2. - np.arccos(cost))*180./np.pi

# Get the density
eq = get_eq(dirname)
rho = (eq.density).reshape((1, nr))

# Get the radial/horizontal integration weights 
gi = GridInfo(dirname + '/grid_info', '')
rw = gi.rweights
tw = (gi.tweights).reshape((nt, 1))

# Get index for v_phi (hopefully you don't change quantity codes partway
# through!
lut = az0.lut
ind_vr = lut[1]
ind_vt = lut[2]
ind_vp = lut[3]

# Set up universal figure dimensions
fig_width_inches = 7. # TOTAL figure width, in inches
    # (i.e., 8x11.5 paper with 1/2-inch margins)
margin_inches = 1./8. # margin width in inches (for both x and y) and 
    # horizontally in between figures
margin_bottom_inches = 0.75
    # larger bottom margin to make room for colorbar(s)
margin_top_inches = 1 # wider top margin to accommodate subplot titles AND metadata
margin_subplot_top_inches = 1/4 # margin to accommodate just subplot titles
#nplots = 4 + 2*magnetism + 1*forced
ncol = 3 # put three plots per row
nrow = 1
#nrow = np.int(np.ceil(nplots/3))

subplot_width_inches = (fig_width_inches - (ncol + 1)*margin_inches)/ncol
    # Make the subplot width so that ncol subplots fit together side-by-side
    # with margins in between them and at the left and right.
subplot_height_inches = 2*subplot_width_inches # Each subplot should have an
    # aspect ratio of y/x = 2/1 to accommodate meridional planes. 
fig_height_inches = margin_top_inches + nrow*(subplot_height_inches +\
        margin_subplot_top_inches + margin_bottom_inches)
fig_aspect = fig_height_inches/fig_width_inches

# "Margin" in "figure units"; figure units extend from 0 to 1 in BOTH 
# directions, so unitless dimensions of margin will be different in x and y
# to force an equal physical margin
margin_x = margin_inches/fig_width_inches
margin_y = margin_inches/fig_height_inches
margin_top = margin_top_inches/fig_height_inches
margin_bottom = margin_bottom_inches/fig_height_inches
margin_subplot_top = margin_subplot_top_inches/fig_height_inches

# Subplot dimensions in figure units
subplot_width = subplot_width_inches/fig_width_inches
subplot_height = subplot_height_inches/fig_height_inches

# Plotting loop
print ('Plotting AZ_Avgs files %s through %s'\
       %(file_list[index_first], file_list[index_last]))
print('Saving plots to '  + plotdir)

# May have to do sliding average (for navg > 1)
# Perform an average of v_phi, then subtract the first AZ_Avg 
# and add the last AZ_Avg with each iteration
# This won't be great if you store multiple records in single AZ_Avg file...

# First time can be initialized here from az0
t1 = az0.time[0]
iter1 = az0.iters[0]

# Initialize the vals array with the average of first navg arrays
print ("Performing initial average over navg = %i files" %navg)
vals = np.zeros((nt, nr, nq))
for i in range(index_first, index_first + navg):
    az = AZ_Avgs(radatadir + file_list[i], '')

    for j in range(az.niter):
        vals += az.vals[:, :, :, j]/(navg*az.niter)
# Second time can be initialized now
t2 = az.time[-1]
iter2 = az.iters[-1]

# Now perform a sliding average
count = 1
for i in range(index_first, index_last + 1):
    if i > index_first: # only past the first point is it necessary to do anything
        print ("Reading AZ_Avgs/", file_list[i])
        print ("Reading AZ_Avgs/", file_list[i + navg - 1])
        az1 = AZ_Avgs(radatadir + file_list[i], '')
        az2 = AZ_Avgs(radatadir + file_list[i + navg - 1], '')
        print ("Subtracting AZ_Avgs/", file_list[i])
        for j in range(az1.niter):
            vals -= az1.vals[:, :, :, j]/(navg*az1.niter)
        print ("Adding AZ_Avgs/", file_list[i + navg - 1])
        for j in range(az2.niter):
            vals += az2.vals[:, :, :, j]/(navg*az2.niter)

        t1 = az1.time[0]
        t2 = az2.time[-1]
        iter1 = az1.iters[0]
        iter2 = az2.iters[-1]
    if i % nskip == 0:
        print ("Plot number %03i" %count)
        count += 1
        # Make the savename like for Mollweide times sample
        savename = 'dr_amom_mc_iter' + str(iter1).zfill(8) + '.png'
        print('Plotting: ' + savename)

        # Get average velocity
        vr_av = vals[:, :, ind_vr]
        vt_av = vals[:, :, ind_vt]
        vp_av = vals[:, :, ind_vp]

        # Get differential rotation in the rotating frame. 
        Om = vp_av/xx
        diffrot = Om*1.0e9/2/np.pi # rad/s --> nHz

        # Get angular momentum in the rotating frame. 
        amom = rho*xx*vp_av

        #  Get DR contrast from 0 to 60 deg lat
        it0, it60 = np.argmin(np.abs(tt_lat)), np.argmin(np.abs(tt_lat - 60))
        Delta_Om = diffrot[it0, 0] - diffrot[it60, 0]

        # Compute the mass flux
        rhovm = rho*np.sqrt(vr_av**2 + vt_av**2)

        # Compute the streamfunction
        psi = streamfunction(rho*vr_av, rho*vt_av, rr, cost)

        # Make CCW negative and CW positive
        rhovm *= np.sign(psi)

        # Get the total integrated absolute amom
        amom_tot = np.sum(np.abs(amom)*tw, axis=0)
        amom_tot = np.sum(amom_tot*rw)
        shell_volume = 4./3.*np.pi*(ro**3. - ri**3.)
        amom_tot *= shell_volume

        # Generate the actual figure of the correct dimensions
        fig = plt.figure(figsize=(fig_width_inches, fig_height_inches))

        # Make 3 subplot axes
        ax1 = fig.add_axes((margin_x, margin_bottom,\
                subplot_width, subplot_height))
        ax2 = fig.add_axes((2*margin_x + subplot_width, margin_bottom,\
                subplot_width, subplot_height))
        ax3 = fig.add_axes((3*margin_x + 2*subplot_width, margin_bottom,\
                subplot_width, subplot_height))

        # Plot the DR with metadata at the top
        plot_azav (diffrot, rr, cost, fig=fig, ax=ax1, units='nHz',\
                nlevs=nlevs, minmax=minmaxdr, rvals=rvals,\
                plotlatlines=plotlatlines, plotboundary=plotboundary)

        # Put directory name in center
        fsize = 12.
        line_height = 1./4./fig_height_inches
        fig.text(margin_x + 0.5*(1 - 2*margin_x), 1 - margin_y, dirname_stripped,\
                 ha='center', va='top', fontsize=fsize, **csfont)

        # Make time label in center
        t_c = (t1 + t2)/2.
        Dt = t2 - t1
        if rotation:
            time_string = ('t = %.1f ' %(t_c/time_unit))\
                    + time_label + (r'$\ (\Delta t = %.2f\ $'\
                    %(Dt/time_unit)) + time_label + ')'
        else:
            time_string = ('t = %.3f' %(t_c/time_unit))\
                    + time_label + (r'$\ \Delta t = %.4f\ $'\
                    %(Dt/time_unit)) + time_label

        fig.text(margin_x + 0.5*(1 - 2*margin_x), 1 - margin_y - line_height,\
                time_string, ha='center', va='top', fontsize=fsize, **csfont)

        # Label DR stuff
        fig.text(margin_x, 1 - margin_y - 2*line_height, r'$\Omega - \Omega_0$',\
                 ha='left', va='top', fontsize=fsize, **csfont)
        fig.text(margin_x, 1 - margin_y - 3*line_height,\
                 r'$\Delta\Omega_{\rm{tot}} = %.1f\ nHz$' %Delta_Om,\
                 ha='left', va='top', fontsize=fsize, **csfont)

        # Make the angular momentum plot
        plot_azav (amom, rr, cost, fig=fig, ax=ax2,\
                units=r'$\rm{g\ cm^{-1}\ s^{-1}}$', nlevs=nlevs, minmax=minmaxamom,\
                rvals=rvals, plotlatlines=plotlatlines, plotboundary=plotboundary)
        # Label amom stuff
        amom_label = r'$\mathcal{L}\equiv \overline{\rho}r\sin\theta\langle v_\phi\rangle$'
        fig.text(2*margin_x + subplot_width, 1 - margin_y - 2*line_height,\
                amom_label, ha='left', va='top', fontsize=fsize, **csfont)
        fig.text(2*margin_x + subplot_width, 1 - margin_y - 3*line_height,\
                 r'$\int|\mathcal{L}|dV = $' + sci_format(amom_tot, 3),\
                 ha='left', va='top', fontsize=fsize, **csfont)

        # Plot mass flux
        plot_azav (rhovm, rr, cost, fig=fig, ax=ax3,\
            units=r'$\rm{g}\ \rm{cm}^{-2}\ \rm{s}^{-1}$', plotcontours=False,\
            minmax=minmaxmc, plotlatlines=plotlatlines, rvals=rvals,\
            plotboundary=plotboundary)

        # Plot streamfunction contours, if desired
        if plotcontours:
            lilbit = 0.01
            maxabs = np.max(np.abs(psi))
            levels = (-maxabs/2., -maxabs/4., -lilbit*maxabs, 0., lilbit*maxabs,\
                    maxabs/4., maxabs/2.)
            plot_azav (psi, rr, cost, fig=fig, ax=ax3, plotfield=False,\
                levels=levels, plotlatlines=plotlatlines,\
                plotboundary=plotboundary)

        # Label MC stuff
        fig.text(3*margin_x + 2*subplot_width, 1 - margin_y - 2*line_height,\
                r'$|\langle\overline{\rho}\mathbf{v}_m\rangle|$',\
                 ha='left', va='top', fontsize=fsize, **csfont)

        # Save figure
        savefile = plotdir + savename
        plt.savefig(savefile, dpi=300)
        plt.close()
        print('----------------------------------')
