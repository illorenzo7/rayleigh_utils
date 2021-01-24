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

import sys, os
sys.path.append(os.environ['rapp'])
sys.path.append(os.environ['raco'])
sys.path.append(os.environ['rapl'])

# Start timing immediately
comm.Barrier()
if rank == 0:
    # timing module
    import time
    # info for print messages
    from common import *
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
from sslice_util import plot_moll
from common import *
from plotcommon import axis_range
from sslice_util import plot_moll
from get_sslice import get_sslice
from varprops import texlabels
from rayleigh_diagnostics import Shell_Slices
reading_func = Shell_Slices
dataname = 'Shell_Slices'

if rank == 0:
    # modules needed only by proc 0 
    from rayleigh_diagnostics import GridInfo

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

    # Directory with data
    radatadir = dirname + '/' + dataname + '/'

    # Get all the file names in datadir and their integer counterparts
    file_list, int_file_list, nfiles = get_file_lists(radatadir)

    # Get specific range desired for plotting
    args = sys.argv[2:]
    nargs = len(args)

    the_tuple = get_desired_range(int_file_list, args)
    if the_tuple is None: # By default plot the last 10 Shell_Slices
        index_first, index_last = nfiles - 11, nfiles - 1  
    else:
        index_first, index_last = the_tuple
    nfiles = index_last - index_first + 1 # this is the number of files we 
    # will plot (modulo the nskip or ntot filters)

    # read command-line arguments
    # defaults
    ir = 0 # by default plot just below the surface
    rval = None # can also find ir by finding the closest point
                # to a local radius divided by rsun
    varlist = ['vr'] # by default plot the radial velocity
    clon = 0.
    ncol = 2 # columns of subplots in figure
    must_smooth = False
    minmax = None
    nskip = 1 # by default don't skip any slices in the range
    ntot = None 
    rootname = None

    # Change defaults
    for i in range(nargs):
        arg = args[i]
        if arg == '-clon':
            clon = float(args[i+1])
        elif arg == '-ir':
            ir = int(args[i+1])
        elif arg == '-rval':
            rval = float(args[i+1])
        elif arg == '-var' or arg == '-qval':
            st = args[i+1]
            if st == 'indr':
                varlist = ['801', '1601', '1602', '1603', '1601plus1602plus1603', '1605']
            elif st == 'indt':
                varlist = ['802', '1606', '1607', '1608', '1606plus1607plus1608', '1610']
            elif st == 'indp':
                varlist = ['803', '1611', '1612', '1613', '1611plus1612plus1613', '1615']
            else:
                varlist = st.split()
        elif arg == '-smooth':
            dlon = int(args[i+1])
            print ("smoothing nonfield vars over %i degrees in lon." %dlon)
            prepend = str(dlon).zfill(3) + 'smooth'
            must_smooth = True
        elif arg == '-minmax':
            minmax = float(args[i+1]), float(args[i+2])
        elif arg == '-nskip':
            nskip = int(args[i+1])
        elif arg == '-ntot':
            ntot = int(args[i+1])
            nskip = nfiles//ntot
        elif arg == '-ncol':
            ncol = int(args[i+1])
        elif arg == '-tag':
            rootname = args[i+1]

    if must_smooth:
        for i in range(len(varlist)):
            var = varlist[i]
            if not var in ['1', '2', '3', '801', '802', '803']:
                varlist[i] = prepend + var

    # get the root name to save plots
    if rootname is None:
        rootname = varlist[0]
        if len(varlist) > 1:
            print ("-tag not specified, so the rootname of the plots will be")
            print ("the first variable in varlist")
            print ("varlist[0] = ", varlist[0])
            print ("you may wish to rerun specifying -tag [helpful_name]")

    # Get the problem size
    nproc_min, nproc_max, n_per_proc_min, n_per_proc_max =\
            opt_workload(nfiles, nproc)

    # Make plot (sub-)directory if it doesn't already exist
    plotdir = dirname + '/plots/moll/times_sample/' + rootname + '/'
    if not os.path.isdir(plotdir):
        os.makedirs(plotdir)

    # Distribute file_list to each process
    for k in range(nproc - 1, -1, -1):
        # distribute the partial file list to other procs 
        if k >= nproc_min: # last processes analyzes more files
            my_nfiles = np.copy(n_per_proc_max)
            istart = nproc_min*n_per_proc_min + (k - nproc_min)*my_nfiles
        else: # first processes analyze fewer files
            my_nfiles = np.copy(n_per_proc_min)
            istart = k*my_nfiles
        iend = istart + my_nfiles

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

    meta = [dirname, radatadir, plotdir, rootname, first_rem, nskip, ir, rval, clon, time_unit, time_label, rotation, minmax, varlist, ncol, dirname_stripped]
else:
    meta = None

dirname, radatadir, plotdir, rootname, first_rem, nskip, ir, rval, clon, time_unit, time_label, rotation, minmax, varlist, ncol, dirname_stripped = comm.bcast(meta, root=0)

# figure dimensions
nplots = len(varlist)
if ncol > nplots:
    ncol = nplots # (no point in having empty columns)
nrow = np.int(np.ceil(nplots/ncol))
fig_width_inches = 6.*ncol

# General parameters for main axis/color bar
margin_bottom_inches = 1./2.
margin_top_inches = 3./4. # margin for big title
margin_subplot_top_inches = 1/4 # margin to accommodate just subplot titles
margin_inches = 1./8.

subplot_width_inches = (fig_width_inches - (ncol + 1)*margin_inches)\
    /ncol
subplot_height_inches = 0.5*subplot_width_inches
fig_height_inches = margin_top_inches + nrow*(margin_bottom_inches +\
    subplot_height_inches + margin_subplot_top_inches)
fig_aspect = fig_height_inches/fig_width_inches

# "Non-dimensional" figure parameters
margin_x = margin_inches/fig_width_inches
margin_y = margin_inches/fig_height_inches
margin_bottom = margin_bottom_inches/fig_height_inches
margin_top = margin_top_inches/fig_height_inches
margin_subplot_top = margin_subplot_top_inches/fig_height_inches

subplot_width = subplot_width_inches/fig_width_inches
subplot_height = subplot_height_inches/fig_height_inches

# Checkpoint and time
comm.Barrier()
if rank == 0:
    t2 = time.time()
    print ('%8.2e s' %(t2 - t1))
    print(fill_str('plotting', lent, char), end='\r')
    t1 = time.time()

# loop over data files and make plots
for i in range(my_nfiles):
    if (i - first_rem) % nskip == 0: # respect nskip!
        a = reading_func(radatadir + str(my_files[i-1]).zfill(8), '')

        # Find desired radius (by default ir=0--near outer surface)
        if i == 0: # only need to do this once
            if not rval is None:
                ir = np.argmin(np.abs(a.radius/rsun - rval))
            rval = a.radius[ir]/rsun 
            # in any case, this is the actual rvalue we get

        # Loop over the slices in the file 
        for j in range(a.niter):
            # Get local time (in seconds)
            t_loc = a.time[j]
            iter_loc = a.iters[j]

            # Savename
            savename = 'moll_' + rootname + ('_rval%0.3f' %rval) + '_iter' +\
                    str(iter_loc).zfill(8) + '.png'

            # make figure
            fig = plt.figure(figsize=(fig_width_inches, fig_height_inches))

            # loop over variables
            for iplot in range(nplots):
                # get the variable's associated field
                varname = varlist[iplot]
                vals = get_sslice(a, varname, dirname=dirname, j=j)
                field = vals[:, :, ir]
           
                # Make axes and plot the Mollweide projection
                ax_left = margin_x + (iplot%ncol)*(subplot_width + margin_x)
                ax_bottom = 1. - margin_top - margin_subplot_top -\
                        subplot_height - (iplot//ncol)*(subplot_height +\
                        margin_subplot_top + margin_bottom)
                ax = fig.add_axes((ax_left, ax_bottom, subplot_width, subplot_height))
                
                plot_moll(field, a.costheta, fig=fig, ax=ax, clon=clon,\
                        varname=varname, minmax=minmax)     

                # label the subplot
                varlabel = texlabels.get(varname, 0)
                if varlabel == 0:
                    varlabel = varname.replace('plus', ' + ')
                    varlabel = varlabel.replace('times', ' ' + r'$\times$' + ' ')
                    varlabel = varlabel.replace('smooth', '')
                    varlabel = varlabel[3:]
                    varlabel = 'qvals = ' + varlabel
                ax.set_title(varlabel, verticalalignment='bottom', **csfont)

                # Make title
                if iplot == 0:
                    ax_xmin, ax_xmax, ax_ymin, ax_ymax = axis_range(ax)
                    ax_delta_x = ax_xmax - ax_xmin
                    ax_delta_y = ax_ymax - ax_ymin
                    ax_center_x = ax_xmin + 0.5*ax_delta_x    

                    if rotation:
                        time_string = ('t = %.1f ' %(t_loc/time_unit)) + time_label +\
                                ' (1 ' + time_label + (' = %.2f days)'\
                                %(time_unit/86400.))
                    else:
                        time_string = ('t = %.3f ' %(t_loc/time_unit)) + time_label +\
                                ' (1 ' + time_label + (' = %.1f days)'\
                                %(time_unit/86400.))

                    title = dirname_stripped +\
                        '\n' + r'$\rm{Mollweide}$' + '     '  + time_string +\
                        '\n' + (r'$r/R_\odot\ =\ %0.3f$' %rval)
                    fig.text(ax_center_x, ax_ymax + 0.02*ax_delta_y +\
                            margin_subplot_top, title,\
                         verticalalignment='bottom', horizontalalignment='center',\
                         fontsize=10, **csfont)   
                
            plt.savefig(plotdir + savename, dpi=300)
            plt.close()

    if rank == 0:
        pcnt_done = i/my_nfiles*100.
        print(fill_str('plotting', lent, char) +\
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
