# Author: Loren Matilsky
# Created: 01/28/2019
# This script plots the axial torques in the meridional plane (viscous, 
# Meridional Circ., Reynolds stress, and Maxwell torques (mean and turbulent) 
# if applicablein the meridional plane 
# ...for the Rayleigh run directory indicated by [dirname]. To use an AZ_Avgs file
# different than the one associated with the longest averaging range, use
# -usefile [complete name of desired AZ_Avgs file]
# Saves plot in
# [dirname]_torque_[first iter]_[last iter].png

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
from rayleigh_diagnostics import AZ_Avgs
from azav_util import plot_azav
from common import get_widest_range_file, strip_dirname, get_dict, get_file_lists, get_desired_range

from get_parameter import get_parameter
from read_inner_vp import read_inner_vp
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
plotdir = dirname + '/plots/torque_times_sample/'
if not os.path.isdir(plotdir):
    os.makedirs(plotdir)

# Read command-line arguments (CLAs)
showplot = True
saveplot = True
plotcontours = True
plotlatlines = True
minmax = None
linthresh = None
linscale = None
minmaxrz = None
linthreshrz = None
linscalerz = None
AZ_Avgs_file = get_widest_range_file(datadir, 'AZ_Avgs')
forced = False
rvals = []
rbcz = None
symlog = False
navg = 1 # by default average over 1 AZ_Avgs instance (no average)
# for navg > 1, a "sliding average" will be used.

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
    if arg == '-minmax':
        minmax = float(args[i+1]), float(args[i+2])
    elif arg == '-minmaxrz':
        minmaxrz = float(args[i+1]), float(args[i+2])
    elif arg == '-rbcz':
        rbcz = float(args[i+1])
    elif arg == '-noshow':
        showplot = False
    elif arg == '-nosave':
        saveplot = False
    elif arg == '-nocontour':
        plotcontours = False
    elif arg == '-usefile':
        AZ_Avgs_file = args[i+1]
        AZ_Avgs_file = AZ_Avgs_file.split('/')[-1]
    elif arg == '-forced':
        forced = True
    elif arg == '-rvals':
        rvals_str = args[i+1].split()
        rvals = []
        for rval_str in rvals_str:
            rvals.append(float(rval_str))
    elif arg == '-symlog':
        symlog = True
    elif arg == '-linthresh':
        linthresh = float(args[i+1])
    elif arg == '-linscale':
        linscale = float(args[i+1])
    elif arg == '-linthreshrz':
        linthreshrz = float(args[i+1])
    elif arg == '-linscalerz':
        linscalerz = float(args[i+1])
    elif arg == '-nolats':
        plotlatlines = False
    elif arg == '-navg':
        navg = int(args[i+1])
        if navg % 2 == 0:
            print ("Please don't enter even values for navg!")
            print ("Replacing navg = %i with navg = %i" %(navg, navg + 1))
            navg += 1

# Check to make sure index_last didn't fall beyond the last possible index ...
if index_last > nfiles - navg:
    index_last = nfiles - navg

# See if magnetism is "on"
try:
    magnetism = get_parameter(dirname, 'magnetism')
except:
    magnetism = False # if magnetism wasn't specified, it must be "off"

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
cost = az0.costheta
sint = az0.sintheta
nr = az0.nr
nt = az0.ntheta
nq = az0.nq
xx = rr.reshape((1, nr))*sint.reshape((nt, 1))

# Get indices for torques (hopefully you don't change quantity codes partway thru!
lut = az0.lut

# Get torque indices
ind_pp = lut[1801]
ind_mm = lut[1802]
ind_cor = lut[1803]
ind_visc = lut[1804]


# Set up universal figure dimensions
fig_width_inches = 7. # TOTAL figure width, in inches
    # (i.e., 8x11.5 paper with 1/2-inch margins)
margin_inches = 1./8. # margin width in inches (for both x and y) and 
    # horizontally in between figures
margin_bottom_inches = 0.75*(2 - (rbcz is None)) 
    # larger bottom margin to make room for colorbar(s)
margin_top_inches = 1 # wider top margin to accommodate subplot titles AND metadata
margin_subplot_top_inches = 1/4 # margin to accommodate just subplot titles
nplots = 4 + 2*magnetism + 1*forced
ncol = 3 # put three plots per row
nrow = np.int(np.ceil(nplots/3))

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
# Perform an average of the torques, then subtract the first AZ_Avg 
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
for i in range(index_first, index_last + 1):
    print ("Plot number %03i" %(i - index_first + 1))
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

    # Make the savename like for Mollweide times sample
    savename = 'torque_iter' + str(iter1).zfill(8) + '.png'
    print('Plotting: ' + savename)

    # Get torques
    torque_rs, torque_mc, torque_visc = -vals[:, :, ind_pp],\
            -vals[:, :, ind_mm] + vals[:, :, ind_cor], vals[:, :, ind_visc]
    torque_tot = torque_rs + torque_mc + torque_visc

    if magnetism:
        ind_Maxwell_mean = lut[1805]
        ind_Maxwell_rs = lut[1806]
       
        torque_Maxwell_mean = vals[:, :, ind_Maxwell_mean]
        torque_Maxwell_rs = vals[:, :, ind_Maxwell_rs]
        
        torque_tot += torque_Maxwell_mean + torque_Maxwell_rs

    if forced: # compute torque associated with forcing function
        mean_vp = vals[:, :, lut[3]]
        tacho_r = get_parameter(dirname, 'tacho_r')
        tacho_dr = get_parameter(dirname, 'tacho_dr')
        tacho_tau = get_parameter(dirname, 'tacho_tau')
        forcing = np.zeros((nt, nr))
        eq = get_eq(dirname)
        rho = eq.rho
        forcing_coeff = -rho/tacho_tau*0.5*(1.0 - np.tanh((rr - tacho_r)/(tacho_dr*rr[0])))
        forcing = forcing_coeff.reshape((1, nr))*mean_vp

        # convert forcing function into a torque
        torque_forcing = xx*forcing
        torque_tot += torque_forcing

    # Set up lists to generate plot
    torques = [torque_rs, torque_mc, torque_visc, torque_tot]
    titles = [r'$\tau_{\rm{rs}}$', r'$\tau_{\rm{mc}}$', r'$\tau_{\rm{v}}$',\
              r'$\tau_{\rm{tot}}$']
    units = r'$\rm{g}\ \rm{cm}^{-1}\ \rm{s}^{-2}$'

    if magnetism:
        torques.insert(3, torque_Maxwell_mean)
        torques.insert(4, torque_Maxwell_rs)
        titles.insert(3, r'$\tau_{\rm{mm}}$')
        titles.insert(4, r'$\tau_{\rm{ms}}$')

    if forced and magnetism:
        torques.insert(5, torque_forcing)
        titles.insert(5, r'$\tau_{\rm{forcing}}$')
    elif forced:
        torques.insert(3, torque_forcing)
        titles.insert(3, r'$\tau_{\rm{forcing}}$')

    # Generate the actual figure of the correct dimensions
    fig = plt.figure(figsize=(fig_width_inches, fig_height_inches))

    for iplot in range(nplots):
        ax_left = margin_x + (iplot%ncol)*(subplot_width + margin_x)
        ax_bottom = 1 - margin_top - subplot_height - margin_subplot_top -\
                (iplot//ncol)*(subplot_height + margin_subplot_top +\
                margin_bottom)
        ax = fig.add_axes((ax_left, ax_bottom, subplot_width, subplot_height))
        plot_azav (torques[iplot], rr, cost, fig=fig, ax=ax, units=units,\
               minmax=minmax, plotcontours=plotcontours, rvals=rvals,\
               minmaxrz=minmaxrz, rbcz=rbcz, symlog=symlog,\
        linthresh=linthresh, linscale=linscale, linthreshrz=linthreshrz,\
        linscalerz=linscalerz, plotlatlines=plotlatlines)

        ax.set_title(titles[iplot], verticalalignment='bottom', **csfont)
    # Make title
    t_c = (t1 + t2)/2.
    Dt = t2 - t1
    if rotation:
        time_string = ('t = %.1f ' %(t_c/time_unit))\
                + time_label + '\n' + (r'$\ (\Delta t = %.2f\ $'\
                %(Dt/time_unit)) + time_label + ')'
    else:
        time_string = ('t = %.3f' %(t_c/time_unit))\
                + time_label + (r'$\ \Delta t = %.4f\ $'\
                %(Dt/time_unit)) + time_label

    # Put some metadata in upper left
    fsize = 12
    fig.text(margin_x, 1 - 0.1*margin_top, dirname_stripped,\
             ha='left', va='top', fontsize=fsize, **csfont)
    fig.text(margin_x, 1 - 0.3*margin_top, 'Torque balance (zonally averaged)',\
             ha='left', va='top', fontsize=fsize, **csfont)
    fig.text(margin_x, 1 - 0.5*margin_top, time_string,\
             ha='left', va='top', fontsize=fsize, **csfont)
    # Save figure
    savefile = plotdir + savename
    plt.savefig(savefile, dpi=300)
    plt.close()
    print('----------------------------------')
