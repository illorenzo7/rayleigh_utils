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
kw_default = dotdict(dict({'rvals': ['rmin', 'rmax'], 'irvals': None, 'lonvals': [-30., 60.], 'eq': True, 'varnames': 'vr', 'clon': 0., 'clat': 20., 't0': False, 'movie': False, 'prepend': False, 'dpi': 300}))

kw_make_figure = kw_make_figure_default.copy()
kw_make_figure.update(ortho_fig_dimensions)

kw_my_contourf = kw_my_contourf_default.copy()
kw_my_contourf.plotcontours = False

kw_range_options = dict({})
for key in range_options: # add the range options key
    kw_range_options[key] = None
kw_range_options.iters = dict({'iters': 'last'})

kw_default.update(kw_make_figure)
kw_default.update(kw_my_contourf)
kw_default.update(kw_range_options)

# now we can update the default kw
if rank == 0:
    print (buff_line)
    find_bad_keys(kw_default, clas, 'plot/slice/cutout_3d', justwarn=True)

# Rayleigh data dirs
ssdir = dirname + '/Shell_Slices/'
merdir = dirname + '/Meridional_Slices/'
eqdir = dirname + '/Equatorial_Slices/'

# update relevant keyword args
kw = update_dict(kw_default, clas)
kw_my_contourf = update_dict(kw_my_contourf, clas)
kw_make_figure = update_dict(kw_make_figure_default, clas)
kw_range_options = update_dict(kw_range_options, clas)

# figure out all the different plots we need
if rank == 0:
    # get desired file names in datadir and their integer counterparts
    # by default, read in last available file
    # For now just work with Shell_Slices files only and assume
    # that Equatorial_Slices and Meridional_Slices overlap
    file_lists, int_file_lists, nfiles = get_file_lists(ssdir, kw_range_options)

    # get desired varnames
    # again assume that the quantity list between Shell_Slices and
    # everything is identical
    sliceinfo_ss = get_sliceinfo(dirname, 'Shell_Slices')
    sliceinfo_mer = get_sliceinfo(dirname, 'Meridional_Slices')
    if isall(kw.varnames):
        kw.varnames = array_of_strings(sliceinfo_ss.qv)
    # no matter what, this needs to be an array of strings
    kw.varnames = array_of_strings(make_array(kw.varnames))
    nq = len(kw.varnames)

    # get gridinfo
    gi = get_grid_info('.')

    # and unpack it
    ntheta = gi.nt
    nphi = 2*ntheta
    lon1d = gi.phi - np.pi
    lat1d = gi.tt_lat*np.pi/180.

    # choose spherical cutout parameters

    # get radian versions of stuff
    kw.lonvals = make_array(kw.lonvals)
    lon1_deg, lon2_deg = kw.lonvals
    lon0, lon1, lon2, lat0 = np.array([kw.clon, lon1_deg, lon2_deg, kw.clat])*np.pi/180.
    # lon1 and lon2 are chosen centered at lon0
    ilon0 = np.argmin(np.abs(lon1d - lon0))
    lon0 = lon1d[ilon0]
    lon1 += lon0
    lon2 += lon0 

    # now put lon1, lon2 in the 0 to 2pi range for convenience
    kw.lonvals = np.array([lon1, lon2])*180./np.pi
    for i in range(2):
        if kw.lonvals[i] < 0.:
            kw.lonvals[i] += 360.
        if kw.lonvals[i] >= 360.:
            kw.lonvals[i] -= 360.

    # get indices and update chosen parameters to correspond to actual grid points

    # radii
    kw.rvals = interpret_rvals(dirname, kw.rvals)
    # get irvals from rvals (default)
    # if irvals were already set, then use those
    if kw.irvals is None:
        kw.irvals = inds_from_vals(sliceinfo_ss.rvals, kw.rvals)
    kw.rvals = sliceinfo_ss.rvals[kw.irvals]
    ir1_ss, ir2_ss = kw.irvals
    r1, r2 = kw.rvals
    ir1 = np.argmin(np.abs(gi.rr - r1))
    ir2 = np.argmin(np.abs(gi.rr - r2))

    # longitudes and latitudes
    ilat0 = np.argmin(np.abs(lat1d - lat0))
    lat0 = lat1d[ilat0]

    # need to shift 
    ilon1_mer, ilon2_mer = inds_from_vals(sliceinfo_mer.lonvals, kw.lonvals)
    kw.lonvals = sliceinfo_mer.lonvals[[ilon1_mer, ilon2_mer]]
    lon1_deg, lon2_deg = kw.lonvals
    lon1, lon2 = kw.lonvals*np.pi/180.
    ilon1 = np.argmin(np.abs(lon1d - lon1))
    ilon2 = np.argmin(np.abs(lon1d - lon2))

    # recompute some other parameters
    beta0 = r2/r1
    lon0_deg, lon1_deg, lon2_deg, lat0_deg  = np.array([lon0, lon1, lon2, lat0])*180./np.pi

    print("r1 =", r1)
    print("r2 =", r2)
    print("(beta0 = " + str(beta0) + ")")
    print("lon0 =", lon0_deg)
    print("lat0 =", lat0_deg)
    print("lon1 =", lon1_deg)
    print("lon2 =", lon2_deg)
    # get desired sampling locations
    # we need 
    # r1, r2 = kw.rvals
    # lon1, lon2 = kw.lonvals
    # lon0 from clon

    kw.ilonvals = inds_from_vals(sliceinfo_ss.rvals, kw.rvals)

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
                        plotdir = my_mkdir(plotdir, erase=kw.movie)

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

# start making plots

linewidth = 1.


kw_make_figure = update_dict(kw_make_figure_default, ortho_fig_dimensions)

plt.close()
fig, axs, fpar = make_figure(**kw_make_figure)
ax = axs[0,0]

# ORTHO 1
# get the first spherical slice
field = np.copy(ss.vals[:, :, ir1_ss, ss.lut[qval], 0])
# shift the field so that the clon is in the ~center of the array
difflon = np.pi - lon0 # basically difflon is the amount the clon
# must be shifted to arrive at 180, which is near the center of array
iphi_shift = int(difflon/(2*np.pi)*nphi)
field = np.roll(field, iphi_shift, axis=0)

# shift the vals in longitude so ilon0 is at lon = 0


# get a "meshgrid" from 1D arrays
lon2d, lat2d = np.meshgrid(lon1d, lat1d, indexing='ij')

# do ortho projection
xx = np.cos(lat2d)*np.sin(lon2d)
yy = np.cos(lat0)*np.sin(lat2d) - np.sin(lat0)*np.cos(lat2d)*np.cos(lon2d)

# mask xx and yy
cond1 = np.sin(lat0)*np.sin(lat2d)+np.cos(lat0)*np.cos(lat2d)*np.cos(lon2d) < 0 
cond2 = (lon2d > lon1 - lon0) & (lon2d < lon2 - lon0) & (lat2d > 0)

field[cond1] = np.nan
field[cond2] = np.nan

kw_my_contourf = dotdict(kw_my_contourf_default.copy())
kw_my_contourf.plotcontours = False

my_contourf(xx, yy, field, fig, ax, **kw_my_contourf)

# ORTHO 2
field = np.copy(ss.vals[:, :, ir2_ss, ss.lut[qval], 0])
# shift the field so that the clon is in the ~center of the array
difflon = np.pi - lon0 # basically difflon is the amount the clon
# must be shifted to arrive at 180, which is near the center of array
iphi_shift = int(difflon/(2*np.pi)*nphi)
field = np.roll(field, iphi_shift, axis=0)

# shift the vals in longitude so ilon0 is at lon = 0

# do ortho projection
xx = beta0*np.cos(lat2d)*np.sin(lon2d)
yy = beta0*(np.cos(lat0)*np.sin(lat2d) - np.sin(lat0)*np.cos(lat2d)*np.cos(lon2d))

# mask xx and yy
cond1 = np.sin(lat0)*np.sin(lat2d)+np.cos(lat0)*np.cos(lat2d)*np.cos(lon2d) < 0 
cond2 = (lon2d < lon1 - lon0) | (lon2d > lon2 - lon0) | (lat2d < 0)

field[cond1] = np.nan
field[cond2] = np.nan

kw_my_contourf.plotcbar = False

my_contourf(xx, yy, field, fig, ax, **kw_my_contourf)


# MERIDIAN 1
field = np.copy(mer.vals[ilon1_mer, :, :, mer.lut[qval], 0]).T

rad1d = np.copy(rr)/r1
rad2d, lat2d = np.meshgrid(rad1d, lat1d, indexing='ij')

xx = - rad2d * np.cos(lat2d) * np.sin(lon0 - lon1)
yy = rad2d * (np.cos(lat0)*np.sin(lat2d) - np.cos(lon0 - lon1)*np.sin(lat0)*np.cos(lat2d))
cond1 = (lat2d < 0) | (rad2d < beta0)  | (rad2d > 1.)
field[cond1] = np.nan

my_contourf(xx, yy, field, fig, ax, **kw_my_contourf)

# MERIDIAN 2
field = np.copy(mer.vals[ilon2_mer, :, ::-1, mer.lut[qval], 0]).T

rad1d = np.copy(rr[::-1])/r1
rad2d, lat2d = np.meshgrid(rad1d, lat1d, indexing='ij')

xx = rad2d * np.cos(lat2d) * np.sin(lon2 - lon0)
yy = rad2d * (np.cos(lat0)*np.sin(lat2d) - np.cos(lon2 - lon0)*np.sin(lat0)*np.cos(lat2d))
cond1 = (lat2d < 0) | (rad2d < beta0)  | (rad2d > 1.)
field[cond1] = np.nan

my_contourf(xx, yy, field, fig, ax, **kw_my_contourf)

# EQUATOR
field = np.copy(eq.vals[:, :, ss.lut[qval], 0])
field = np.roll(field, iphi_shift, axis=0)
rad1d = np.copy(rr)/r1
lon2d, rad2d = np.meshgrid(lon1d, rad1d, indexing='ij')

xx = rad2d*np.sin(lon2d)
yy = -rad2d*np.sin(lat0)*np.cos(lon2d)
cond = (lon2d < lon1 - lon0) | (lon2d > lon2 - lon0) | (rad2d < beta0)  | (rad2d > 1.)
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
ax.plot(np.cos(svals)*np.sin(lon1-lon0),\
        np.cos(lat0)*np.sin(svals) - np.sin(lat0)*np.cos(svals)*np.cos(lon1-lon0),\
        'k-', linewidth=linewidth)

ax.plot(beta0*np.cos(svals)*np.sin(lon1-lon0),\
        beta0*(np.cos(lat0)*np.sin(svals) - np.sin(lat0)*np.cos(svals)*np.cos(lon1-lon0)),\
        'k-', linewidth=linewidth)

ax.plot(np.cos(svals)*np.sin(lon2-lon0),\
        np.cos(lat0)*np.sin(svals) - np.sin(lat0)*np.cos(svals)*np.cos(lon2-lon0),\
        'k-', linewidth=linewidth)

ax.plot(beta0*np.cos(svals)*np.sin(lon2-lon0),\
        beta0*(np.cos(lat0)*np.sin(svals) - np.sin(lat0)*np.cos(svals)*np.cos(lon2-lon0)),\
        'k-', linewidth=linewidth)

svals = np.linspace(beta0, 1, nsvals)

ax.plot(svals*np.sin(lon1-lon0), -svals*np.sin(lat0)*np.cos(lon1-lon0),\
        'k-', linewidth=linewidth)

ax.plot(svals*np.sin(lon2-lon0), -svals*np.sin(lat0)*np.cos(lon2-lon0),\
        'k-', linewidth=linewidth)

ax.plot(0*svals, svals*np.cos(lat0), 'k-', linewidth=linewidth)

svals = np.linspace(lon1 - lon0, lon2 - lon0, nsvals)
ax.plot(np.sin(svals), -np.sin(lat0)*np.cos(svals), 'k-', linewidth=linewidth)
ax.plot(beta0*np.sin(svals), -beta0*np.sin(lat0)*np.cos(svals), 'k-', linewidth=linewidth)


plt.show()
