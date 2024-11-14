# Author: Loren Matilsky
# Date created: 11/29/2022
# Description: go-to slice plotting routine for everything

# initialize communication
from mpi4py import MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()

# import common
import sys, os
sys.path.append(os.environ['raco'])
from common import *

# start the clock
comm.Barrier()
if rank == 0:
    import time
    nproc = comm.Get_size()
    t1_glob = time.time()
    t1 = t1_glob + 0.0
    if nproc > 1:
        print ('processing in parallel with %i ranks' %nproc)
        print ('communication initialized')
    else:
        print ('processing in serial with 1 rank')
    print(fill_str('proc 0 preparing problem size'))

# additional modules needed
import matplotlib.pyplot as plt
import numpy as np
sys.path.append(os.environ['rapp'])
sys.path.append(os.environ['rapl'])
sys.path.append(os.environ['rapl'] + '/azav')
from plotcommon import *
from cla_util import *
from slice_util import *
from rayleigh_diagnostics import Equatorial_Slices, Meridional_Slices
from rayleigh_diagnostics_alt import Shell_Slices
from azav_util import plot_azav, plot_azav_kwargs_default, azav_fig_dimensions
#from get_slice import get_slice, get_label

# Get CLAs
args = sys.argv
clas0, clas = read_clas(args)
dirname = clas0.dirname
dirname_stripped = strip_dirname(dirname)

# SPECIFIC ARGS
kwargs_default = dotdict(dict({'type': None, 'isamplevals': 0, 'samplevals': None, 'rvals': None, 'varnames': 'vr', 'clons': None, 'clats': None, 'clonrange': None, 'clatrange': None, 't0': False, 'movie': False, 'moviesampleval': False, 'movieclon': False, 'movieclat': False, 'shrink': False, 'prepend': False, 'dpi': 300, 'lonav': False}))

# these guys need to update right away to choose fig dimensions
if clas.type is None:
    plottype = 'moll'
else:
    plottype = clas.type
if 'lonav'in clas:
    kwargs_default.lonav = True
lonav = kwargs_default.lonav

if plottype == 'moll':
    if lonav:
        fig_dimensions = molllonav_fig_dimensions
    else:
        fig_dimensions = moll_fig_dimensions
if plottype == 'ortho':
    fig_dimensions = ortho_fig_dimensions
if plottype in ['moll', 'ortho']:
    plotting_func = plot_moll_or_ortho
    plotting_func_kwargs_default = plot_moll_or_ortho_kwargs_default
    if plottype == 'ortho':
        plotting_func_kwargs_default['ortho'] = True
    dataname = 'Shell_Slices'
    reading_func = Shell_Slices
    samplelabel = 'rval'
    samplefmt = '%1.3e'
if plottype == 'eq':
    fig_dimensions = eq_fig_dimensions
    plotting_func = plot_eq
    plotting_func_kwargs_default = plot_eq_kwargs_default
    dataname = 'Equatorial_Slices'
    reading_func = Equatorial_Slices
if plottype == 'mer':
    fig_dimensions = azav_fig_dimensions
    plotting_func = plot_azav
    plotting_func_kwargs_default = plot_azav_kwargs_default
    plotting_func_kwargs_default['plotcontours'] = False
    dataname = 'Meridional_Slices'
    reading_func = Meridional_Slices
    samplelabel = 'lonval'
    samplefmt = lon_fmt

# Rayleigh data dir
radatadir = dirname + '/' + dataname + '/'

# now we can update the default kwargs
kwargs_default.update(plotting_func_kwargs_default)
make_figure_kwargs_default.update(fig_dimensions)
kwargs_default.update(make_figure_kwargs_default)
if rank == 0:
    for key in range_options: # add the range options key
        kwargs_default[key] = None
    print (buff_line)
    print ("plottype = " + plottype)
    find_bad_keys(kwargs_default, clas, 'plot/slice/raslice', justwarn=True)

# update relevant keyword args
kw = update_dict(kwargs_default, clas)
kw_plotting_func = update_dict(plotting_func_kwargs_default, clas)
kw_make_figure = update_dict(make_figure_kwargs_default, clas)

# figure out all the different plots we need
if rank == 0:
    # get desired file names in datadir and their integer counterparts
    # by default, read in last available file
    clas_mod = dict({'iters': 'last'})
    clas_mod.update(clas)
    file_list, int_file_list, nfiles = get_file_lists(radatadir, clas_mod)

    # figure out what clons and clats we need
    kw.clons = make_array(kw.clons)
    if not kw.clonrange is None:
        clonmin, clonmax, nclon = kw.clonrange
        kw.clons = np.linspace(clonmin, clonmax, nclon)
    if kw.clons is None:
        kw.clons = np.array([0.0])
    nclon = len(kw.clons)

    kw.clats = make_array(kw.clats)
    if not kw.clatrange is None:
        clatmin, clatmax, nclat = kw.clatrange
        kw.clats = np.linspace(clatmin, clatmax, nclat)
    if kw.clats is None:
        kw.clats = np.array([20.0])
    nclat = len(kw.clats)

    # get desired varnames
    sliceinfo = get_sliceinfo(dirname, dataname)
    if isall(kw.varnames):
        kw.varnames = array_of_strings(sliceinfo.qv)
    # no matter what, this needs to be an array of strings string
    kw.varnames = array_of_strings(make_array(kw.varnames))
    nq = len(kw.varnames)

    # get desired sampling locations
    # note: for equatorial slices it doesn't matter -- only one equator!

    # first the indices must be an array
    kw.isamplevals = make_array(kw.isamplevals)

    # can control samplevals with rvals for time-latitude traces
    if not kw.rvals is None:
        kw.samplevals = kw.rvals

    # get the samplevals we want (not needed for equatorial slices)
    if not plottype == 'eq':
        if not kw.samplevals is None: # samplevals have been set directly
            # need the available sampling locations
            if isall(kw.samplevals):
                kw.isamplevals = np.arange(sliceinfo.nsamplevals)
            else:
                kw.samplevals = make_array(kw.samplevals)
                kw.isamplevals = inds_from_vals(sliceinfo.samplevals, kw.samplevals)

    # these are the sample vals we end up with
    if plottype == 'eq':
        kw.samplevals = np.array([0.0])
    else:
        kw.samplevals = make_array(sliceinfo.samplevals[kw.isamplevals])
    nsamplevals = len(kw.samplevals)

    # say what we are plotting
    print (buff_line)

    # print file list
    # see if user specified a range or if it was default
    its_default = True
    for key in clas.keys():
        if key in range_options:
            its_default = False
    if 'iters' in clas.keys() or its_default:
        print ('plotting the following files:', arr_to_str(file_list, '%s'))
    else:
        print ('plotting %i %s files:\n%s through %s'\
            %(nfiles, dataname, file_list[0], file_list[-1]))

    # print varnames
    print (buff_line)
    print (("plotting %i variables:\nvarnames = " %nq) +\
            arr_to_str(kw.varnames, "%s"))

    # (possibly) clons and clats
    if not plottype == 'mer':
        # print clons
        print (buff_line)
        print (("plotting %i central longitudes:\nclons = " %nclon) +\
                arr_to_str(kw.clons, lon_fmt))
    if plottype == 'ortho':
        # print clats
        print (buff_line)
        print (("plotting %i central latitudes:\nclats = " %nclat) +\
                arr_to_str(kw.clats, lat_fmt))

    # (possibly) sampling locations
    if not plottype == 'eq':
        print (buff_line)
        print ("plotting %i sampling locations:" %nsamplevals)
        print ("i%ss = " %samplelabel + arr_to_str(kw.isamplevals, '%i'))
        print ("%ss = " %samplelabel + arr_to_str(kw.samplevals, samplefmt))

    # calculate total number of figures
    nfigures = nclat*nclon*nsamplevals*nq*nfiles
    print (buff_line)
    print ("nfigures = %i x %i x %i x %i x %i = %i"\
            %(nclat, nclon, nsamplevals, nq, nfiles, nfigures))
    print (buff_line)

    # prepare the epic loop!
    plotting_instructions = []
    if kw.movie:
        count = 0
    for fname in file_list:
        if kw.movie:
            count += 1

        for varname in kw.varnames:
            basic = is_basic(varname)
            if basic:
                varlabel = get_label(varname)
                simple_label = varname
            else:
                varlabel, simple_label = get_label(varname)

            if kw.movieclon:
                count = 0
            for clon in kw.clons:
                if kw.movieclon:
                    count += 1

                if kw.movieclat:
                    count = 0
                for clat in kw.clats:
                    if kw.movieclat:
                        count += 1

                    if kw.moviesampleval:
                        count = 0
                    for isampleval in kw.isamplevals:
                        if kw.moviesampleval:
                            count += 1
                        # get the sample val
                        if plottype == 'eq':
                            sampleval = 0. # (just a placeholder)
                        else:
                            sampleval = sliceinfo.samplevals[isampleval]

                        # now we're in the big loop
                        # get plot directory and image name to save the file
                        if kw.movie:
                            plotdir = 'movie_' + plottype + '_time/' + simple_label
                            if not plottype == 'mer':
                                plotdir += ('_clon' + lon_fmt) %clon
                            if plottype == 'ortho':
                                plotdir += ('_ccolat' + lon_fmt) %(90.0 - clat)
                            if not plottype == 'eq':
                                plotdir += ('_' + samplelabel + samplefmt) %sampleval
                            savename = 'img%04i.png' %count
                        elif kw.moviesampleval:
                            plotdir = 'movie_' + plottype + sample_label + '/' + fname + '_' + simple_label
                            if not plottype == 'mer':
                                plotdir += ('_clon' + lon_fmt) %clon
                            if plottype == 'ortho':
                                plotdir += ('_ccolat' + lon_fmt) %(90.0 - clat)
                            savename = 'img%04i.png' %count
                        elif kw.movieclon:
                            plotdir = 'movie_' + plottype + '_clon/' + fname + '_' + simple_label
                            if plottype == 'ortho':
                                plotdir += ('_ccolat' + lon_fmt) %(90.0 - clat)
                            if not plottype == 'eq':
                                plotdir += ('_' + samplelabel + samplefmt) %sampleval
                        elif kw.movieclat:
                            plotdir = 'movie_' + plottype + '_clat/' + fname + '_' + simple_label
                            if not plottype == 'mer':
                                plotdir += ('_clon' + lon_fmt) %clon
                            if not plottype == 'eq':
                                plotdir += ('_' + samplelabel + samplefmt) %sampleval
                            savename = 'img%04i.png' %count
                        else:
                            plotdir = clas0['plotdir'] + plottype + clas0['tag']
                            if lonav:
                                savename = plottype + 'lonav_' + fname + '_' + simple_label 
                            else:
                                savename = plottype + '_' + fname + '_' + simple_label 
                            if kw.prepend:
                                savename = dirname_stripped + '_' + savename
                            if not plottype == 'mer':
                                savename += ('_clon' + lon_fmt) %clon
                            if plottype == 'ortho':
                                savename += ('_ccolat' + lon_fmt) %(90.0 - clat)
                            if not plottype == 'eq':
                                savename += ('_' + samplelabel + samplefmt) %sampleval
                            savename += '.png'

                        savefile = plotdir + '/' + savename
                        plotdir = my_mkdir(plotdir)

                        plotting_instructions.append([fname,\
                                varname,\
                                clon,\
                                clat,\
                                isampleval,\
                                sampleval,\
                                savefile,\
                                varlabel])


# Checkpoint
comm.Barrier()
if rank == 0:
    print(fill_str('proc 0 distributing the plotting instructions'), end='')
                        
# distribute the plotting instructions
if rank == 0:
    # get the problem size
    nproc_min, nproc_max, n_per_proc_min, n_per_proc_max =\
            opt_workload(nfigures, nproc)

    # distribute plotting instructions to each process
    for k in range(nproc - 1, -1, -1):
        # distribute the partial file list to other procs 
        if k < nproc_max: # first processes analyzes more files
            my_nfigures = np.copy(n_per_proc_max)
            istart = k*my_nfigures
            iend = istart + my_nfigures
        else: # last processes analyze fewer files
            my_nfigures = np.copy(n_per_proc_min)
            istart = nproc_max*n_per_proc_max + (k - nproc_max)*my_nfigures
            iend = istart + my_nfigures

        # get the file list portion for rank k
        my_instructions = plotting_instructions[istart:iend]
        # send  my_files, my_nfigures if nproc > 1
        if k >= 1:
            comm.send([my_instructions, my_nfigures], dest=k)
else: # recieve my_files, my_nfigures
    my_instructions, my_nfigures = comm.recv(source=0)

# Checkpoint and time
comm.Barrier()
if rank == 0:
    t2 = time.time()
    print (format_time(t2 - t1))
    print(fill_str('beginning the plotting job'))
    t1 = time.time()

# now loop over and plot figures
for ifigure in range(my_nfigures):
    # local instructions for this plot
    fname, varname, clon, clat, isampleval, sampleval, savefile, varlabel =\
            my_instructions[ifigure]

    if not plottype == 'mer':
        kw_plotting_func.clon = clon        
    if not plottype in ['mer', 'eq']:
        kw_plotting_func.clat = clat

    # get the slice
    if plottype in ['moll', 'ortho']:
        varname_root, deriv, primevar, sphvar = get_varprops(varname)
        qvals = None # by default, but update
        if varname_root in var_indices:
            qvals = var_indices[varname_root]
        elif varname_root == 'omz':
            qvals = [301, 302]
        a = reading_func(radatadir + fname, '', irvals=isampleval, qvals=qvals)
    else:
        a = reading_func(radatadir + fname, '')

    # get the variable
    vals = get_slice(a, varname, dirname=dirname)

    # get sample location
    if dataname == 'Shell_Slices':
        field = vals[:, :, 0] # only read in one sample val
    elif dataname == 'Meridional_Slices':
        field = vals[isampleval, :, :]
    else: # equatorial slices
        field = vals

    # make plot
    fig, axs, fpar = make_figure(**kw_make_figure)
    ax = axs[0, 0]
    if plottype in ['moll', 'ortho']:
        plotting_args = field, a.costheta, fig, ax
    elif plottype == 'mer':
        plotting_args = field, a.radius, a.costheta, fig, ax
    elif plottype == 'eq':
        plotting_args = field, a.radius, fig, ax

    plotting_func(*plotting_args, **kw_plotting_func)

    # possibly add lon avg.
    if lonav:
        fig_width_inches, fig_height_inches = fig.get_size_inches()
        fig_aspect = fig_height_inches/fig_width_inches
        # get ax dimensions
        ax_left, ax_right, ax_bottom, ax_top = axis_range(ax)
        ax_width = ax_right - ax_left
        ax_height = ax_top - ax_bottom

        ax_line_left = ax_right + 1/8/fig_width_inches
        ax_line_bottom = ax_bottom
        ax_line_width = 1 - 5/8/fig_width_inches - ax_line_left
        ax_line_height = ax_height
        ax_line = fig.add_axes([ax_line_left, ax_line_bottom, ax_line_width, ax_line_height])

        # make line plot
        field_lonav = np.mean(field, axis=0)
        tt_lat = (np.pi/2. - np.arccos(a.costheta))*180./np.pi
        ax_line.plot(field_lonav, tt_lat)

        # set ylim
        ax_line.set_ylim(-90., 90.)

        # mark zero lines
        mark_axis_vals(ax_line, 'x')
        mark_axis_vals(ax_line, 'y')

        # move ticks to right
        ax_line.yaxis.tick_right()
        ax_line.yaxis.set_label_position("right")

        # Get ticks everywhere
        plt.sca(ax_line)
        plt.minorticks_on()
        plt.tick_params(top=True, right=True, direction='in', which='both')

        # label the axes
        ax_line.set_title('lon. avg.')
        ax_line.set_ylabel('latitude (degrees)')


    # make title
    if kw.t0:
        iter0 = int(fname)
        floatwidth = 6
    else:
        iter0 = None
        floatwidth = None

    time_string = get_time_string(dirname, a.iters[0], iter0=iter0, floatwidth=floatwidth)

    if plottype == 'moll':
        location_and_perspective =\
            (samplelabel + ' = ' + samplefmt) %sampleval +\
                ('\n' + r'$\phi_0$' + ' = ' + lon_fmt_tex) %clon
    elif plottype == 'ortho':
        location_and_perspective =\
            (samplelabel + ' = ' + samplefmt) %sampleval +\
                ('\n' + r'$\phi_0$' + ' = ' + lon_fmt_tex) %clon +\
                ('    ' + r'$\lambda_0$' + ' = ' + lat_fmt_tex) %clat
    elif plottype == 'eq':
        location_and_perspective = (r'$\phi_0$' + ' = ' + lon_fmt_tex) %clon
    elif plottype == 'mer':
        location_and_perspective = (samplelabel + ' = ' + samplefmt) %sampleval

    title = dirname_stripped + '\n' +\
            varlabel + '\n' +\
            location_and_perspective + '\n' +\
            time_string
    ax.set_title(title, va='bottom', fontsize=default_titlesize)

    # save by default
    if clas0['saveplot']:
        plt.savefig(savefile, dpi=kw.dpi)
    if rank == 0:
        pcnt_done = (ifigure+1)/my_nfigures*100.
        # print what we saved and how far along we are
        if clas0['saveplot']:
            print(fill_str("saved " + savefile + " (dpi=%i)" %kw.dpi) + '\n' +\
                    ('rank 0 %5.1f%% done' %pcnt_done), end='\r')
        # always show if nfigures is 1
        if nfigures == 1 and clas0['showplot']:
            plt.show()   
    # close figure at end of loop
    plt.close()

# Checkpoint and time
comm.Barrier()
if rank == 0:
    t2 = time.time()
    print('\n' + fill_str('plotting time'), end='')
    print (format_time(t2 - t1))
    print(make_bold(fill_str('total time')), end='')
    print (make_bold(format_time(t2 - t1_glob)))
    print (buff_line)
