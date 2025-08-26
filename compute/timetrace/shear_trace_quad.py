##################################################################
# Routine to trace Rayleigh data in time, in different quadrants
# of the meridional plane
##################################################################
# This routine computes the trace in time of global averages in various 
# quadrants
#
# By default the "quadrant" is the entire plane and 
# G_Avgs are used
#
# if --nquadr or --rvals or is specified (but not --nquadlat or --latvals)
# the meridional plane is divided into spherical shells and 
# Shell_Avgs are used
#
# if --nquadlat or --latvals is specified, meridional plane is divided 
# into conic (and/or shellular) sections  and
# AZ_Avgs are used
#
# vals is shape
# [ntimes, nq, nquadlat, nquadr] 
# quadrants are ordered in increasing latitude (South to North) 
# and increasing radius (inner shell to outer shell) 
##################################################################

# import common first
import sys, os
sys.path.append(os.environ['raco'])
from common import *

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
    from cla_util import *
    nproc = comm.Get_size()
    t1_glob = time.time()
    t1 = t1_glob + 0.0
    if nproc > 1:
        print ('processing in parallel with %i ranks' %nproc)
        print ('communication initialized')
    else:
        print ('processing in serial with 1 rank')
    print(fill_str('processes importing necessary modules'), end='')

# modules needed by everyone
import numpy as np
# data type and reading function
import sys, os
sys.path.append(os.environ['raco'])
sys.path.append(os.environ['rapp'])
from rayleigh_diagnostics import G_Avgs, Shell_Avgs, AZ_Avgs
from common import inds_from_vals
from grid_util import *

if rank == 0:
    # modules needed only by proc 0 
    import pickle

# Checkpoint and time
comm.Barrier()
if rank == 0:
    t2 = time.time()
    print (format_time(t2 - t1))
    print(fill_str('proc 0 distributing the file lists'), end='')
    t1 = time.time()

# proc 0 reads the file lists and distributes them
if rank == 0:
    # read the arguments
    args = sys.argv
    clas0, clas = read_clas(args)
    dirname = clas0['dirname']
    eq = get_eq(dirname)
    omega0 = eq.omega0
    tag = clas0['tag']

    # get grid info + default kwargs
    kw_default = dict({})
    kw_default['rvals'] = None # can specify radial domain boundaries directly 
    kw_default['latvals'] = None
    kw_default['thinby'] = None # by default don't thin end series

    # update these possibly
    kw = update_dict(kw_default, clas)

    # get desired quadrant boundaries

    # radial boundaries
    if kw.rvals is None: 
        # this is the default, most easily overwritten
        # with --nquadr
        rmin, rmax = get_rminmax(clas0.dirname)
        rvals = np.array([rmin, rmax])
    else:
        rvals = np.sort(kw.rvals)

    # latitudinal boundaries
    if kw.latvals is None:
        # this is the default, most easily overwritten
        # with --nquadlat
        latmin, latmax = get_latminmax(clas0.dirname)
        latvals = np.array([latmin, latmax])
    else:
        latvals = np.sort(kw.latvals)

    # always use the az averages
    dataname = 'AZ_Avgs'

    # update the number of quadrants
    nquadlat = len(latvals) - 1
    nquadr = len(rvals) - 1
    nquad = nquadlat*nquadr

    # compute the volumes of each quadrant
    volumes = np.zeros((nquadlat, nquadr))
    for ilat in range(nquadlat):
        cost1 = np.cos(np.pi/2. - latvals[ilat]*np.pi/180.)
        cost2 = np.cos(np.pi/2. - latvals[ilat+1]*np.pi/180.)
        for ir in range(nquadr):
            r1, r2 = rvals[ir], rvals[ir+1]
            volumes[ilat, ir] = 2.*np.pi/3.*(r2**3. - r1**3.)*(cost2 - cost1)

    # Get the Rayleigh data directory
    radatadir = dirname + '/' + dataname + '/'

    # Get all the file names in datadir and their integer counterparts
    file_list, int_file_list, nfiles = get_file_lists(radatadir, clas)

    # Get the problem size
    nproc_min, nproc_max, n_per_proc_min, n_per_proc_max =\
            opt_workload(nfiles, nproc)

    # Distribute file_list to each process
    for k in range(nproc - 1, -1, -1):
        # distribute the partial file list to other procs 
        if k < nproc_max: # first processes analyzes more files
            my_nfiles = np.copy(n_per_proc_max)
            istart = k*my_nfiles
            iend = istart + my_nfiles
        else: # last processes analyze fewer files
            my_nfiles = np.copy(n_per_proc_min)
            istart = nproc_max*n_per_proc_max + (k - nproc_max)*my_nfiles
            iend = istart + my_nfiles

        # Get the file list portion for rank k
        my_files = np.copy(int_file_list[istart:iend])

        # send  my_files, my_nfiles if nproc > 1
        if k >= 1:
            comm.send([my_files, my_nfiles], dest=k)
else: # recieve my_files, my_nfiles
    my_files, my_nfiles = comm.recv(source=0)

# Broadcast dirname, radatadir, etc.
if rank == 0:
    meta = [\
dirname, dataname, radatadir, latvals, rvals, omega0]
else:
    meta = None
dirname, dataname, radatadir, latvals, rvals, omega0 = comm.bcast(meta, root=0)

# figure out which reading_func to use
reading_func = AZ_Avgs

# Checkpoint and time
comm.Barrier()
if rank == 0:
    t2 = time.time()
    print (format_time(t2 - t1))
    print ("tracing over %i x %i = %i quadrants" %(nquadlat, nquadr, nquad))
    print ("rvals = " + arr_to_str(rvals, "%1.3e"))
    print ("latvals = " + arr_to_str(latvals, "%.1f"))
    print ('Considering %i %s files for the trace: %s through %s'\
        %(nfiles, dataname, file_list[0], file_list[-1]))
    print(fill_str('computing'), end='\r')
    t1 = time.time()

# Now analyze the data (some processes may not have nquadlat, nquadr, etc.)
nquadlat = len(latvals) - 1
nquadr = len(rvals) - 1
my_times = []
my_iters = []
my_vals = []
the_qv = 0 # need this just as a place holder

for i in range(my_nfiles):
    a = reading_func(radatadir + str(my_files[i]).zfill(8), '')

    # want to deal here with the possibility of a changing grid
    # (changing n_theta, ncheby)
    # thus need to compute ilatvals, irvals, tw, rw using current 
    # data slice, instead of above

    # get the grid from the current data slice
    # only needed if not a global average
    nt = a.ntheta
    ncheby, domain_bounds = get_ncheby_from_rr(a.radius)
    nr = a.nr
    rr, rw, tt, tw = compute_grid_info(ncheby, domain_bounds, nt)
    tt_lat = 180./np.pi*(np.pi/2. - tt)
    rr_2d = rr.reshape((1, nr))
    tt_2d = tt.reshape((nt, 1))
    xx = rr_2d*np.sin(tt_2d)

    irvals = inds_from_vals(rr, rvals)
    ilatvals = np.sort(inds_from_vals(tt_lat, latvals))

    for j in range(a.niter):
        vals_loc = a.vals[..., j]


        # Get the values in the separate quadrants
        vals_gav = np.zeros((3, nquadlat, nquadr))
        # dOmega, grad Omega, \Delta
        # note: the default is nquadlat, nquadr = 1, 1 (yes "extra" dimensions)
        for ilat in range(nquadlat): # remember: tt_lat is INCREASING
            it1 = ilatvals[ilat]
            it2 = ilatvals[ilat + 1]

            for ir in range(nquadr): # remember: rvals increase but r-inds decrease
                # compute the shear properties
                ir1 = irvals[ir + 1]
                ir2 = irvals[ir]

                # compute quadrant weights
                tw_quad = (tw[it1:it2+1]/np.sum(tw[it1:it2+1])).\
                        reshape((it2 - it1 + 1, 1))
                rw_quad = (rw[ir1:ir2+1]/np.sum(rw[ir1:ir2+1])).\
                        reshape((1, ir2 - ir1 + 1))
                weights_quad = tw_quad*rw_quad

                # get rotation rate
                uphi = np.copy(vals_loc[it1:it2+1, ir1:ir2+1, a.lut[3]])
                uphi /= 2*omega0
                xx_quad = np.copy(xx[it1:it2+1, ir1:ir2+1])
                rr_quad = np.copy(rr[ir1:ir2+1])
                tt_quad = np.copy(tt[it1:it2+1])
                omega = uphi/xx_quad

                # get standard deviation
                omega_mean = np.sum(weights_quad*omega)
                omega_std = (omega - omega_mean)**2 # get integrand first
                omega_std = np.sqrt(np.sum(weights_quad*omega_std))

                # compute magnitude of gradient of omega
                domdr = drad(omega, rr_quad)
                domdt = dth(omega, tt_quad)/rr_2d[:, ir1:ir2+1]
                gradom = np.sqrt(domdr**2 + domdt**2)
                gradom = np.sum(weights_quad*gradom)

                vals_gav[:, ilat, ir] = omega_std, gradom, omega_std/gradom

        my_vals.append(vals_gav)
        my_times.append(a.time[j])
        my_iters.append(a.iters[j])

    if rank == 0:
        pcnt_done = i/my_nfiles*100.
        print(fill_str('computing') +\
                ('rank 0 %5.1f%% done' %pcnt_done), end='\r')

# Checkpoint and time
comm.Barrier()
if rank == 0:
    t2 = time.time()
    print('\n' + fill_str('computing time'), end='')
    print (format_time(t2 - t1))
    print(fill_str('rank 0 collecting and saving the results'), end='')
    t1 = time.time()

# proc 0 now collects the results from each process
if rank == 0:
    # initialize empty lists for the final arrays
    times = []
    iters = []
    vals = []

    # Gather the results into these "master" arrays
    for j in range(nproc):
        if j >= 1:
            # Get my_my_times, my_iters, my_vals from rank j
            my_times, my_iters, my_vals = comm.recv(source=j)
        times += my_times
        iters += my_iters
        vals += my_vals
else: # other processes send their data
    comm.send([my_times, my_iters, my_vals], dest=0)

# Make sure proc 0 collects all data
comm.Barrier()

# proc 0 saves the data
if rank == 0:
    # create data directory if it doesn't already exist
    datadir = clas0['datadir']
    if not os.path.isdir(datadir):
        os.makedirs(datadir)

    # make things into arrays
    vals = np.array(vals)
    times = np.array(times)
    iters = np.array(iters)

    # Set  savename by what we are saving, quadrant info, 
    # and first and last iteration files for the trace
    basename = 'shear_trace_nquadlat%i_nquadr%i' %(nquadlat, nquadr)
    savename = basename + tag + '-' + file_list[0] + '_' + file_list[-1] + '.pkl'
    savefile = datadir + savename

    # compute the full average over the whole volume (really a sanity check)
    vals_full = np.zeros(np.shape(vals[:, :, 0, 0]))
    for ilat in range(nquadlat):
        for ir in range(nquadr):
            vals_full += volumes[ilat, ir]*vals[:, :, ilat, ir]
    vals_full /= np.sum(volumes)

    # possibly thin data
    if not kw.thinby is None:
        times = times[::kw.thinby]
        iters = iters[::kw.thinby]
        vals = vals[::kw.thinby]
        vals_full = vals_full[::kw.thinby]

    # save the data
    # Get first and last iters of files
    iter1, iter2 = int_file_list[0], int_file_list[-1]
    f = open(savefile, 'wb')
    di_sav = {'vals': vals, 'times': times, 'iters': iters,  'volumes': volumes, 'volume_full': np.sum(volumes), 'rvals': rvals, 'latvals': latvals}
    if 'nquad' in savefile:
        di_sav[ 'vals_full'] = vals_full
    pickle.dump(di_sav, f, protocol=4)
    f.close()
    t2 = time.time()
    print (format_time(t2 - t1))
    print(make_bold(fill_str('total time')), end='')
    print (make_bold(format_time(t2 - t1_glob)))
    print ('data saved at ')
    print (make_bold(savefile))
