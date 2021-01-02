# Author: Loren Matilsky
# Created: 07/02/2020
# Modified + Parallelized: 01/02/2021
# This script plots the DR and angular momentum in the meridional plane
# for many different times (to make up the frames of a movie)
# Saves plots in subdirectory
# plots/diffrot_amom_times_sample/

# initialize communication
from mpi4py import MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()

# Start timing immediately
comm.Barrier()
if rank == 0:
    # timing module
    import time
    # info for print messages
    import sys, os
    sys.path.append(os.environ['raco'])
    from common import fill_str
    from mpi_util import opt_workload
    lent = 50
    char = '.'
    nproc = comm.Get_size()
    t1_glob = time.time()
    t1 = t1_glob + 0.0
    if nproc > 1:
        print ('processing in parallel with %i ranks' %nproc)
        print ('communication initialized')
    else:
        print ('processing in serial with 1 rank')
    print(fill_str('processes importing necessary modules', lent, char),\
            end='')

# modules needed by everyone
import numpy as np
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
from azav_util import plot_azav, streamfunction
from common import sci_format
from rayleigh_diagnostics import AZ_Avgs
reading_func = AZ_Avgs
dataname = 'AZ_Avgs'

if rank == 0:
    # modules needed only by proc 0 
    from common import strip_dirname, get_file_lists, get_desired_range
    from get_parameter import get_parameter
    from time_scales import compute_Prot, compute_tdt
    from rayleigh_diagnostics import AZ_Avgs, GridInfo

# Checkpoint and time
comm.Barrier()
if rank == 0:
    t2 = time.time()
    print ('%8.2e s' %(t2 - t1))
    print(fill_str('proc 0 distributing the file lists', lent, char),\
            end='')
    t1 = time.time()

if rank == 0:
    # Get directory name and stripped_dirname for plotting purposes
    dirname = sys.argv[1]
    dirname_stripped = strip_dirname(dirname)

    # Directory with AZ_Avgs
    radatadir = dirname + '/' + dataname + '/'

    # Get all the file names in datadir and their integer counterparts
    file_list, int_file_list, nfiles = get_file_lists(radatadir)

    # read command-line arguments
    # defaults
    args = sys.argv[2:]
    nargs = len(args)
    plotcontours = True
    plotlatlines = False
    rvals = []
    rbcz = None
    minmax = None
    minmaxrz = None
    symlog = False
    linthresh = None
    linscale = None
    linthreshrz = None
    linscalerz = None
    navg = 1 # by default average over 1 AZ_Avgs instance (no average)
    # for navg > 1, a "sliding average" will be used.
    tag = '' # optional way to tag save directory
    nlevs = 20
    plotboundary = True
    nskip = 1 # by default don't skip anything
    ntot = None 

    # Change defaults
    for i in range(nargs):
        arg = args[i]
        if arg == '-minmax':
            minmax = float(args[i+1]), float(args[i+2])
        elif arg == '-minmaxrz':
            minmaxrz = float(args[i+1]), float(args[i+2])
        elif arg == '-rbcz':
            rbcz = float(args[i+1])
        elif arg == '-nocontour':
            plotcontours = False
        elif arg == '-rvals':
            rvals_str = args[i+1].split()
            rvals = []
            for rval_str in rvals_str:
                rvals.append(float(rval_str))
        elif arg == '-lats':
            plotlatlines = True
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
        elif arg == '-tag':
            tag = '_' + args[i+1]

    # Get files we want
    the_tuple = get_desired_range(int_file_list, args)
    if the_tuple is None: # By default average over the first 50 files
        index_first, index_last = 0, 49
    else:
        index_first, index_last = the_tuple
    # Check to make sure index_last didn't fall beyond the last possible
    # index
    if index_last > nfiles - navg:
        index_last = nfiles - navg

    # Remove parts of file lists we don't need
    #file_list = file_list[index_first:index_last + 1]
    file_list = file_list[index_first:index_last + navg]
    #int_file_list = int_file_list[index_first:index_last + 1]
    int_file_list = int_file_list[index_first:index_last + navg]
    nfiles = index_last - index_first + 1

    # Get the problem size
    nproc_min, nproc_max, n_per_proc_min, n_per_proc_max =\
            opt_workload(nfiles, nproc)

    # Create plotdir if doesn't already exist
    plotdir = dirname +  '/torque_times_sample' + tag + '/'
    if not os.path.isdir(plotdir):
        os.makedirs(plotdir)

    # Distribute file_list to each process
    for k in range(nproc - 1, -1, -1):
        # distribute the partial file list to other procs 
        if k >= nproc_min: # last processes analyzes more files
            my_nfiles = np.copy(n_per_proc_max)
            istart = nproc_min*n_per_proc_min + (k - nproc_min)*my_nfiles
            iend = istart + my_nfiles + navg - 1
        else: # first processes analyze fewer files
            my_nfiles = np.copy(n_per_proc_min)
            istart = k*my_nfiles
            iend = istart + my_nfiles + navg - 1

        # Get the file list portion for rank k
        my_files = np.copy(int_file_list[istart:iend])
        
        # remainder of the start index (for nskip)
        first_rem = (istart - index_first) % nskip

        # send  my_files, my_nfiles, first_rem if nproc > 1
        if k >= 1:
            comm.send([my_files, my_nfiles, first_rem], dest=k)
else: # recieve my_files, my_nfiles
    my_files, my_nfiles, first_rem = comm.recv(source=0)

# collect and broadcast metadata
if rank == 0:
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
    # See if magnetism is "on"
    magnetism = get_parameter(dirname, 'magnetism')

    # Will need the first data file for a number of things
    a0 = reading_func(radatadir + file_list[0], '')

    # Get necessary grid info
    rr = a0.radius
    ri, ro = np.min(rr), np.max(rr)
    shell_volume = 4./3.*np.pi*(ro**3. - ri**3.)
    cost = a0.costheta
    sint = a0.sintheta
    nr = a0.nr
    nt = a0.ntheta
    nq = a0.nq
    xx = rr.reshape((1, nr))*sint.reshape((nt, 1))
    tt_lat = (np.pi/2. - np.arccos(cost))*180./np.pi

    # Get the radial/horizontal integration weights 
    gi = GridInfo(dirname + '/grid_info', '')
    rw = gi.rweights
    tw = (gi.tweights).reshape((nt, 1))
    meta = [dirname, radatadir, plotdir, nt, nr, nq, first_rem, nskip, xx,\
            rr, rbcz, cost, shell_volume, rw, tw, tt_lat, time_unit,\
            time_label, navg, plotcontours, plotlatlines, plotboundary,\
            rvals, nlevs, minmaxrz, minmax, dirname_stripped, rotation, magnetism,\
            symlog, linthresh, linthreshrz, linscale, linscalerz]
else:
    meta = None
dirname, radatadir, plotdir, nt, nr, nq, first_rem, nskip, xx,\
        rr, rbcz, cost, shell_volume, rw, tw, tt_lat, time_unit,\
        time_label, navg, plotcontours, plotlatlines, plotboundary,\
        rvals, nlevs, minmaxrz, minmax, dirname_stripped, rotation, magnetism,\
        symlog, linthresh, linthreshrz, linscale, linscalerz =\
            comm.bcast(meta, root=0)

# Checkpoint and time
comm.Barrier()
if rank == 0:
    t2 = time.time()
    print ('%8.2e s' %(t2 - t1))
    print(fill_str('plotting', lent, char), end='\r')
    t1 = time.time()

# make plots
# Set up figure dimensions
fig_width_inches = 7. # TOTAL figure width, in inches
    # (i.e., 8x11.5 paper with 1/2-inch margins)
margin_inches = 1./8. # margin width in inches (for both x and y) and 
    # horizontally in between figures
margin_bottom_inches = 0.75*(2 - (rbcz is None)) 
    # larger bottom margin to make room for colorbar(s)
margin_top_inches = 1 # wider top margin to accommodate subplot titles AND metadata
margin_subplot_top_inches = 1/4 # margin to accommodate just subplot titles
nplots = 4 + 2*magnetism 
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

# May have to do sliding average (for navg > 1)
# Perform an average then subtract the first AZ_Avg 
# and add the next AZ_Avg with each iteration
# This won't be great if you store multiple records in single AZ_Avg file
# (the average will slide more)

# Initialize the vals array with the average of first navg arrays
vals = np.zeros((nt, nr, nq))
for i in range(navg):
    a = reading_func(radatadir + str(my_files[i]).zfill(8), '')
    # initialize first time
    if i == 0:
        t1 = a.time[0]
        iter1 = a.iters[0]
    for j in range(a.niter):
        vals += a.vals[:, :, :, j]/(navg*a.niter)
# initialize second time here
t2 = a.time[-1]
iter2 = a.iters[-1]

# Now perform a sliding average
for i in range(my_nfiles):
    if i > 0: # only past the first point is it necessary to do anything
        a = reading_func(radatadir + str(my_files[i-1]).zfill(8), '')
        a2 = reading_func(radatadir + str(my_files[i+navg-1]).zfill(8), '')
        for j in range(a.niter):
            vals -= a.vals[:, :, :, j]/(navg*a.niter)
        for j in range(a2.niter):
            vals += a2.vals[:, :, :, j]/(navg*a2.niter)

        t1 = a.time[0]
        t2 = a2.time[-1]
        iter1 = a.iters[0]
        iter2 = a2.iters[-1]
    if (i - first_rem) % nskip == 0:
        # Make the savename like for Mollweide times sample
        savename = 'torque_iter' + str(my_files[i]).zfill(8) + '.png'

        # Get torques
        ind_pp = a.lut[1801]
        ind_mm = a.lut[1802]
        ind_cor = a.lut[1803]
        ind_visc = a.lut[1804]
        torque_rs, torque_mc, torque_visc = -vals[:, :, ind_pp],\
                -vals[:, :, ind_mm] + vals[:, :, ind_cor],\
                vals[:, :, ind_visc]
        torque_tot = torque_rs + torque_mc + torque_visc

        if magnetism:
            ind_Maxwell_mean = a.lut[1805]
            ind_Maxwell_rs = a.lut[1806]
           
            torque_Maxwell_mean = vals[:, :, ind_Maxwell_mean]
            torque_Maxwell_rs = vals[:, :, ind_Maxwell_rs]
            
            torque_tot += torque_Maxwell_mean + torque_Maxwell_rs

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
        fig.text(margin_x, 1 - 0.3*margin_top, 'Torque balance (zonal average)',
                 ha='left', va='top', fontsize=fsize, **csfont)
        fig.text(margin_x, 1 - 0.5*margin_top, time_string,\
                 ha='left', va='top', fontsize=fsize, **csfont)
        # Save figure
        savefile = plotdir + savename
        plt.savefig(savefile, dpi=300)
        plt.close()
    if rank == 0:
        pcnt_done = i/my_nfiles*100.
        print(fill_str('computing', lent, char) +\
            ('rank 0 %5.1f%% done' %pcnt_done), end='\r')

# Checkpoint and time
comm.Barrier()
if rank == 0:
    t2 = time.time()
    print(fill_str('\nplotting time', lent, char), end='')
    print ('%8.2e s                                 ' %(t2 - t1))
    print(fill_str('total time', lent, char), end='')
    print ('%8.2e s' %(t2 - t1_glob))
    print ('view frames using ')
    print ('eog ' + plotdir)
