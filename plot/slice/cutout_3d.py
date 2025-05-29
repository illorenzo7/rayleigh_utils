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
kw_default = dotdict(dict({'t0': False, 'movie': False, 'prepend': False, 'dpi': 300}))

kw_make_figure = kw_make_figure_default.copy()
kw_make_figure.update(ortho_fig_dimensions)

kw_my_contourf = kw_my_contourf_default.copy()
kw_my_contourf.plotcontours = False

kw_default.update(kw_plot_cutout_3d_default)
kw_default.update(kw_make_figure)
kw_default.update(kw_my_contourf)
kw_default.update(kw_range_options_default)

# now we can update the default kw
if rank == 0:
    print (buff_line)
    find_bad_keys(kw_default, clas, 'plot/slice/cutout_3d', justwarn=True)

# update relevant keyword args
kw = update_dict(kw_default, clas)
kw_plot_cutout_3d = update_dict(kw_plot_cutout_3d_default, clas)
kw_my_contourf = update_dict(kw_my_contourf, clas)
kw_make_figure = update_dict(kw_make_figure_default, clas)
kw_range_options = update_dict(kw_range_options_default, clas)

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
    file_lists, int_file_lists, nfiles = get_file_lists(ssdir, kw_range_options)

    # get desired varnames
    # again assume that the quantity list between Shell_Slices and
    # everything else is identical
    sliceinfo_ss = get_sliceinfo(dirname, 'Shell_Slices')
    if isall(kw.varnames):
        kw.varnames = array_of_strings(sliceinfo_ss.qv)
    # no matter what, this needs to be an array of strings
    kw.varnames = array_of_strings(make_array(kw.varnames))
    nq = len(kw.varnames)

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

    # calculate total number of figures
    nfigures = nq*nfiles
    print (buff_line)
    print ("nfigures = %i x %i = %i" %(nq, nfiles, nfigures))
    print (buff_line)

    # prepare the epic loop!
    plotting_instructions = []
    count = 0
    for fname in file_list:
        # get gridinfo (do this individually for each slice, in case
        # resolution changes)
        sliceinfo_ss = get_sliceinfo(dirname, 'Shell_Slices', fname=fname)
        sliceinfo_mer = get_sliceinfo(dirname, 'Meridional_Slices', fname=fname)
        gi = get_grid_info(dirname)

        # and unpack it
        ntheta = gi.nt
        nphi = gi.nphi
        lon1d = gi.phi - np.pi
        lat1d = gi.tt_lat*np.pi/180.

        # choose spherical cutout parameters

        # get radian versions of stuff
        clon, clat = np.array([kw.clon, kw.clat])*np.pi/180.

        # get indices and recompute values
        # longitudes and latitudes
        iclat = iclosest(lat1d, clat)
        clat = lat1d[iclat]
        iclon = iclosest(lon1d, clon)
        clon = lon1d[iclon]

        # user might specify lon1 and lon2 via dlon1 and dlon2
        # (the distance from clon)
        if not kw.dlon1 is None:
            kw.lon1 = clon + kw.dlon1
        if not kw.dlon2 is None:
            kw.lon2 = clon + kw.dlon2

        # get indices and recompute values for r1, r2, lon1, lon2
        r1, r2 = interpret_rvals(dirname, np.array([kw.r1, kw.r2]))
        rvals_ss = sliceinfo_ss.rvals
        ir1_ss, ir2_ss = inds_from_vals(rvals_ss, [r1, r2])
        r1, r2 = rvals_ss[[ir1_ss, ir2_ss]]
        ir1, ir2 = inds_from_vals(gi.rr, np.array([r1, r2]))
        beta = r1/r2 # this is the aspect ratio of the spherical cutout 

        # need to shift available lons from Meridional Slices
        # to lie in range [-pi, pi)
        lonvals_mer = np.copy(sliceinfo_mer.lonvals)
        for ilon in range(len(lonvals_mer)):
            if lonvals_mer[ilon] >= np.pi:
                lonvals_mer[ilon] -= 2*np.pi
        ilon1_mer, ilon2_mer = inds_from_vals(lonvals_mer, [lon1, lon2])
        lon1, lon2 = lonvals_mer[ilon1_mer, ilon2_mer]
        ilon1, ilon2 = inds_from_vals(lon1d, [lon1, lon2])

        # get degree versions of things
        clon_deg, clat_deg, lon1_deg, lon2_deg = np.array([clon, clat, lon1, lon2])*180./np.pi

        if count == 0:
            # print the parameters of the cutouts
            print(buff_line)
            print("PLOTTING SPHERICAL CUTOUTS!")
            print("clon =", clon_deg)
            print("clat =", clat_deg)
            print("r1 =", r1)
            print("r2 =", r2)
            print("(beta = " + str(beta) + ")")
            print("lon1 =", lon1_deg)
            print("lon2 =", lon2_deg)

        # loop over the variable names and get plotting instructions
        for varname in kw.varnames:
            if is_basic(varname):
                varlabel = get_label(varname)
                simple_label = varname
            else:
                varlabel, simple_label = get_label(varname)

            # get metadata labels
            meta_label = simple_label +\
                ('_ccolat' + lon_fmt) %(90. - clat_deg) +\
                ('_clon' + lon_fmt) %clon_deg +\
                ('_lon1' + lon_fmt) %lon1_deg +\
                ('_lon2' + lon_fmt) %lon2_deg +\
                ('_r1' + flt_fmt) %r1 +\
                ('_r2' + flt_fmt) %r2
            if kw.eq:
                    meta_label += '_witheq'
                else:
                    meta_label += '_noeq'

            # now we're in the big loop
            # get plot directory and image name to save the file
            if kw.movie:
                plotdir = 'movie_cut3d/' + meta_label
                savename = 'img%04i.png' %(count+1)
            else:
                plotdir = clas0['plotdir'] + '/cut3d' + clas0['tag']
                savename = 'cut3d_' + fname + '_' + meta_label + '.png'
                if kw.prepend:
                    savename = dirname_stripped + '_' + savename
                savename += '.png'

            savefile = plotdir + '/' + savename

            plotdir = my_mkdir(plotdir, erase=kw.movie)

            plotting_instructions.append([fname,\
                    varname,\
                    iclon,\
                    clat,\
                    ilon1_mer,\
                    ilon2_mer,\
                    ilon1,\
                    ilon2,\
                    beta,\
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

# start making plots

# spherical cutout plotting routine
kw__default = dict({'clon': 0., 'clat': 20., 'shrinkage': 1., 'plotlonlines': True, 'lonvals': np.arange(0., 360., 30.), 'plotlatlines': True, 'latvals': np.arange(-60., 90., 30.), 'linewidth': 0.75*default_lw, 'plotboundary': True, 'ortho': False})
kw_plot_moll_or_ortho_default.update(kw_my_contourf_default)

dict({'clon': 0., 'clat': 20., 'shrinkage': 1., 'plotlonlines': True, 'lonvals': np.arange(0., 360., 30.), 'plotlatlines': True, 'latvals': np.arange(-60., 90., 30.), 'linewidth': 0.75*default_lw, 'plotboundary': True, 'ortho': False})
kw_plot_cutout_3d_default = dotdict(dict({'r1': 'rmin', 'r2': 'rmax', 'lon1': -30., 'lon2': 60., 'dlon1': None, 'dlon2': None, 'eq': True, 'varnames': 'vr', 'clon': 0., 'clat': 20., 't0': False}))
kw_plot_cutout3d_default.update(kw_my_contourf_default)
kw_plot_cutout3d_default.plotcontours = False

def plot_cutout_3d(dirname, fname, varname, fig, ax, **kw):

    # get the Shell_Slice data
    ss 

    # ORTHO 1
    # get the first spherical slice
    field = np.copy(ss.vals[:, :, ir1_ss, ss.lut[qval], 0])
    # shift the field so that the clon is in the ~center of the array
    difflon = np.pi - clon # basically difflon is the amount the clon
    # must be shifted to arrive at 180, which is near the center of array
    iphi_shift = int(difflon/(2*np.pi)*nphi)
    field = np.roll(field, iphi_shift, axis=0)

    # shift the vals in longitude so iclon is at lon = 0


    # get a "meshgrid" from 1D arrays
    dlon2d, lat2d = np.meshgrid(lon1d, lat1d, indexing='ij')

    # do ortho projection
    xx = np.cos(lat2d)*np.sin(dlon2d)
    yy = np.cos(clat)*np.sin(lat2d) - np.sin(clat)*np.cos(lat2d)*np.cos(dlon2d)

    # mask xx and yy
    cond1 = np.sin(clat)*np.sin(lat2d)+np.cos(clat)*np.cos(lat2d)*np.cos(dlon2d) < 0 
    cond2 = (dlon2d > dlon1 - clon) & (dlon2d < dlon2 - clon) & (lat2d > 0)

    field[cond1] = np.nan
    field[cond2] = np.nan

    kw_my_contourf = dotdict(kw_my_contourf_default.copy())
    kw_my_contourf.plotcontours = False

    my_contourf(xx, yy, field, fig, ax, **kw_my_contourf)

    # ORTHO 2
    field = np.copy(ss.vals[:, :, ir2_ss, ss.lut[qval], 0])
    # shift the field so that the clon is in the ~center of the array
    difflon = np.pi - clon # basically difflon is the amount the clon
    # must be shifted to arrive at 180, which is near the center of array
    iphi_shift = int(difflon/(2*np.pi)*nphi)
    field = np.roll(field, iphi_shift, axis=0)

    # shift the vals in longitude so iclon is at lon = 0

    # do ortho projection
    xx = beta*np.cos(lat2d)*np.sin(dlon2d)
    yy = beta*(np.cos(clat)*np.sin(lat2d) - np.sin(clat)*np.cos(lat2d)*np.cos(dlon2d))

    # mask xx and yy
    cond1 = np.sin(clat)*np.sin(lat2d)+np.cos(clat)*np.cos(lat2d)*np.cos(dlon2d) < 0 
    cond2 = (dlon2d < dlon1 - clon) | (dlon2d > dlon2 - clon) | (lat2d < 0)

    field[cond1] = np.nan
    field[cond2] = np.nan

    kw_my_contourf.plotcbar = False

    my_contourf(xx, yy, field, fig, ax, **kw_my_contourf)


    # MERIDIAN 1
    field = np.copy(mer.vals[idlon1_mer, :, :, mer.lut[qval], 0]).T

    rad1d = np.copy(rr)/r1
    rad2d, lat2d = np.meshgrid(rad1d, lat1d, indexing='ij')

    xx = - rad2d * np.cos(lat2d) * np.sin(clon - dlon1)
    yy = rad2d * (np.cos(clat)*np.sin(lat2d) - np.cos(clon - dlon1)*np.sin(clat)*np.cos(lat2d))
    cond1 = (lat2d < 0) | (rad2d < beta)  | (rad2d > 1.)
    field[cond1] = np.nan

    my_contourf(xx, yy, field, fig, ax, **kw_my_contourf)

    # MERIDIAN 2
    field = np.copy(mer.vals[idlon2_mer, :, ::-1, mer.lut[qval], 0]).T

    rad1d = np.copy(rr[::-1])/r1
    rad2d, lat2d = np.meshgrid(rad1d, lat1d, indexing='ij')

    xx = rad2d * np.cos(lat2d) * np.sin(dlon2 - clon)
    yy = rad2d * (np.cos(clat)*np.sin(lat2d) - np.cos(dlon2 - clon)*np.sin(clat)*np.cos(lat2d))
    cond1 = (lat2d < 0) | (rad2d < beta)  | (rad2d > 1.)
    field[cond1] = np.nan

    my_contourf(xx, yy, field, fig, ax, **kw_my_contourf)

    # EQUATOR
    field = np.copy(eq.vals[:, :, ss.lut[qval], 0])
    field = np.roll(field, iphi_shift, axis=0)
    rad1d = np.copy(rr)/r1
    dlon2d, rad2d = np.meshgrid(lon1d, rad1d, indexing='ij')

    xx = rad2d*np.sin(dlon2d)
    yy = -rad2d*np.sin(clat)*np.cos(dlon2d)
    cond = (dlon2d < dlon1 - clon) | (dlon2d > dlon2 - clon) | (rad2d < beta)  | (rad2d > 1.)
    field[cond] = np.nan

    my_contourf(xx, yy, field, fig, ax, **kw_my_contourf)


    lilbit = 0.01
    ax.set_xlim(-1-lilbit, 1 + lilbit)
    ax.set_ylim(-1-lilbit, 1+lilbit)

    # plot the boundary
    nsvals = 1000
    svals = np.linspace(0, 2*np.pi, nsvals)
    ax.plot(np.cos(svals), np.sin(svals), 'k-', linewidth=linewidth)

    svals = np.linspace(0., np.pi/2.)
    ax.plot(np.cos(svals)*np.sin(dlon1-clon),\
            np.cos(clat)*np.sin(svals) - np.sin(clat)*np.cos(svals)*np.cos(dlon1-clon),\
            'k-', linewidth=linewidth)

    ax.plot(beta*np.cos(svals)*np.sin(dlon1-clon),\
            beta*(np.cos(clat)*np.sin(svals) - np.sin(clat)*np.cos(svals)*np.cos(dlon1-clon)),\
            'k-', linewidth=linewidth)

    ax.plot(np.cos(svals)*np.sin(dlon2-clon),\
            np.cos(clat)*np.sin(svals) - np.sin(clat)*np.cos(svals)*np.cos(dlon2-clon),\
            'k-', linewidth=linewidth)

    ax.plot(beta*np.cos(svals)*np.sin(dlon2-clon),\
            beta*(np.cos(clat)*np.sin(svals) - np.sin(clat)*np.cos(svals)*np.cos(dlon2-clon)),\
            'k-', linewidth=linewidth)

    svals = np.linspace(beta, 1, nsvals)

    ax.plot(svals*np.sin(dlon1-clon), -svals*np.sin(clat)*np.cos(dlon1-clon),\
            'k-', linewidth=linewidth)

    ax.plot(svals*np.sin(dlon2-clon), -svals*np.sin(clat)*np.cos(dlon2-clon),\
            'k-', linewidth=linewidth)

    ax.plot(0*svals, svals*np.cos(clat), 'k-', linewidth=linewidth)

    svals = np.linspace(dlon1 - clon, dlon2 - clon, nsvals)
    ax.plot(np.sin(svals), -np.sin(clat)*np.cos(svals), 'k-', linewidth=linewidth)
    ax.plot(beta*np.sin(svals), -beta*np.sin(clat)*np.cos(svals), 'k-', linewidth=linewidth)


    plt.show()
