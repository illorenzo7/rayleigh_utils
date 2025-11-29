# Author: Loren Matilsky
# Date created: 05/04/2025
# Description: 3d-cutout slice plotting routine

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
from azav_util import plot_azav, kw_plot_azav_default, azav_fig_dimensions
#from get_slice import get_slice, get_label

# Get CLAs
args = sys.argv
clas0, clas = read_clas(args)
dirname = clas0.dirname
dirname_stripped = strip_dirname(dirname)

# SPECIFIC ARGS
kw_default = dotdict({'t0': False, 'prepend': False, 'dpi': 300, 'varnames': 'vr', 'varnames2': 'omzprime', 'movie': False})

kw_make_figure = dotdict(kw_make_figure_default)
kw_make_figure.update(ortho_fig_dimensions)
kw_make_figure.sub_margin_top_inches = 1.25 # more room for labels
if 'movie' in clas.keys():
    movie = kw_default.movie = True
else:
    movie = kw_default.movie = False # movie is now a global parameter for all processes
    kw_make_figure.sub_margin_top_inches = 3/8 # more room for labels


kw_default.update(kw_plot_cutout_3d_default)
kw_default.update(kw_make_figure)
kw_default.update(kw_range_options_default) # this is only to suppress the keyword warning

# now we can update the default kw
if rank == 0:
    print (buff_line)
    find_bad_keys(kw_default, clas, 'plot/slice/cutout_3d', justwarn=True)

# update relevant keyword args
kw = update_dict(kw_default, clas)
kw_plot_cutout_3d = update_dict(kw_plot_cutout_3d_default, clas)
# need to adjust bottom margin depending on number of colorbars
if kw_plot_cutout_3d.numcbar == 0:
    kw_make_figure.sub_margin_bottom_inches = 0.
elif kw_plot_cutout_3d.onecbar:
    kw_make_figure.sub_margin_bottom_inches = 1/4
elif kw_plot_cutout_3d.twocbar:
    kw_make_figure.sub_margin_bottom_inches = 1/2

kw_make_figure = update_dict(kw_make_figure, clas)

# Rayleigh data dirs
ssdir = dirname + '/Shell_Slices/'
merdir = dirname + '/Meridional_Slices/'
if kw.eq:
    eqdir = dirname + '/Equatorial_Slices/'


# figure out all the different plots we need
if rank == 0:
    # get desired file names in datadir and their integer counterparts
    # by default, read in last available file
    # For now just work with Shell_Slices files only and assume
    # that Equatorial_Slices and Meridional_Slices overlap
    file_list, int_file_list, nfiles = get_file_lists(ssdir, clas)

    # get desired varnames
    # again assume that the quantity list between Shell_Slices and
    # everything else is identical
    sliceinfo_ss = get_sliceinfo(dirname, 'Shell_Slices')
    if isall(kw.varnames):
        kw.varnames = array_of_strings(sliceinfo_ss.qv)
    # no matter what, this needs to be an array of strings
    kw.varnames = array_of_strings(make_array(kw.varnames))
    nq = len(kw.varnames)
    if kw.twovars:
        kw.varnames2 = array_of_strings(make_array(kw.varnames2))


    # print file list
    print(buff_line)
    if len(file_list) == 1:
        print ('plotting 1 temporal slice:', arr_to_str(file_list, '%s'))
    else:
        print ('plotting %i %s temporal slices:\n%s through %s'\
            %(nfiles, dataname, file_list[0], file_list[-1]))

    # print varnames
    print (buff_line)
    print (("plotting %i variables:\nvarnames = " %nq) +\
            arr_to_str(kw.varnames, "%s"))
    if kw.twovars:
        print ("with different variables in the inner part:" +\
                arr_to_str(kw.varnames2, "%s"))

    # calculate total number of figures
    nfigures = nq*nfiles
    print (buff_line)
    print ("nfigures = %i x %i = %i" %(nq, nfiles, nfigures))
    print (buff_line)

    # check if we need to store the first time, 
    if movie:
        t0 = translate_times(int_file_list[0], dirname).time

    # prepare the epic loop!
    plotting_instructions = []
    if kw.movie:
        count = 0
    print("flist=", file_list)
    for fname in file_list:
        if kw.movie:
            count += 1
        
        ivar = 0
        for varname in kw.varnames:
            basic = is_basic(varname)
            if basic:
                varlabel = get_label(varname)
                simple_label = varname
            else:
                varlabel, simple_label = get_label(varname)

            if kw.twovars:
                varname2 = kw.varnames2[ivar]
                basic = is_basic(varname2)
                if basic:
                    varlabel2 = get_label(varname2)
                    simple_label2 = varname2
                else:
                    varlabel2, simple_label2 = get_label(varname2)
                varlabel = varlabel + ' and ' + varlabel2
                simple_label = simple_label + '_and_' + simple_label2


            # get metadata labels
            r1, r2 = interpret_rvals(dirname, [kw.r1, kw.r2])
            meta_label = simple_label +\
                ('_clat' + lat_fmt) %kw.clat +\
                ('_clon' + lat_fmt) %kw.clon +\
                '_r1_%.3f' %r1 +\
                '_r2_%.3f' %r2 +\
                ('_dlon1' + lat_fmt) %kw.dlon1 +\
                ('_dlon2' + lat_fmt) %kw.dlon2
            if kw.eq:
                meta_label += '_witheq'
            else:
                meta_label += '_noeq'

            # meta data, for titling the plots
            location_and_perspective =\
                    (r'$\phi_0$' + ' = ' + lat_fmt_tex) %kw.clon +\
                    ('    ' + r'$\chi_0$' + ' = ' + lat_fmt_tex) %kw.clat +\
                    ('\n' + r'$r_1=%0.3f$' %r1) +\
                    ('    ' + r'$r_2=%0.3f$' %r2) +\
                    ('\n' + r'$\Delta\phi_1$' + ' = ' + lat_fmt_tex) %kw.dlon1 +\
                    ('    ' + r'$\Delta\phi_2$' + ' = ' + lat_fmt_tex) %np.abs(kw.dlon2)


            # get plot directory and image name to save the file
            if kw.movie:
                plotdir = 'movie_cut3d/' + meta_label
                savename = 'img%04i.png' %(count+1)
            else:
                plotdir = clas0['plotdir'] + '/cut3d' + clas0['tag']
                savename = 'cut3d_' + fname + '_' + meta_label + '.png'
                if kw.prepend:
                    savename = dirname_stripped + '_' + savename

            savefile = plotdir + '/' + savename

            plotdir = my_mkdir(plotdir, erase=kw.movie)

            varnames_to_plot = [varname]
            if kw.twovars:
                varnames_to_plot += [varname2]

            to_append = [fname, varnames_to_plot, savefile,\
                    varlabel, kw.movie]
            if kw.movie:
                to_append += [t0]

            plotting_instructions.append(to_append)

            ivar += 1

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
            comm.send([my_instructions, my_nfigures, location_and_perspective], dest=k)
else: # recieve my_files, my_nfigures
    my_instructions, my_nfigures, location_and_perspective = comm.recv(source=0)

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
    the_instructions = my_instructions[ifigure]
    if len(the_instructions) == 5:
        fname, varnames_to_plot, savefile, varlabel, movie = the_instructions
    elif len(the_instructions) == 6:
        fname, varnames_to_plot, savefile, varlabel, movie, t0 = the_instructions

    if len(varnames_to_plot) == 1:
        varname = varnames_to_plot[0]
        twovars = False
    elif len(varnames_to_plot) == 2:
        varname, varname2 = varnames_to_plot
        twovars = True

    # make plot
    fig, axs, fpar = make_figure(**kw_make_figure)
    ax = axs[0, 0]

    if twovars:
        kw_plot_cutout_3d.varname2 = varname2
    plot_cutout_3d(dirname, fname, varname, fig, ax, **kw_plot_cutout_3d)

    # label the plot
    the_time = translate_times(int(fname), dirname).time
    SF = 5
    nolabel = False
    if movie:
        the_time -= t0
        SF = 3
        nolabel = True

    time_string = get_time_string(dirname, t1=the_time,SF=SF, nolabel=nolabel)

    if movie:
        title = varlabel + ' '*5 + time_string
    else:
        title = dirname_stripped + '\n' +\
            varlabel + '\n' +\
            time_string + '\n' +\
            location_and_perspective 
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
