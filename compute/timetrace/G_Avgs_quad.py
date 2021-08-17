##################################################################
# Routine to trace Rayleigh G_Avgs data in time
# Author: Loren Matilsky
# Created: 12/18/2018
# Parallelized: 11/26/2020
##################################################################
# This routine computes the trace in time of the values in the G_Avgs data 
# for a particular simulation. 
# This routine computes the trace in time of the values in the G_Avgs data 
# for a particular simulation. 
# traces separtely over different "quadrants" (4-8 in total, depending on
# if rbcz is specified to separate the domains in radius)
##################################################################

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
    # import common here
    from common import *
    from cla_util import *
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
# data type and reading function
import sys, os
sys.path.append(os.environ['rapp'])
from rayleigh_diagnostics import AZ_Avgs
reading_func = AZ_Avgs
dataname = 'AZ_Avgs'

if rank == 0:
    # modules needed only by proc 0 
    import pickle

# Checkpoint and time
comm.Barrier()
if rank == 0:
    t2 = time.time()
    print (format_time(t2 - t1))
    print(fill_str('proc 0 distributing the file lists', lent, char),\
            end='')
    t1 = time.time()

# proc 0 reads the file lists and distributes them
if rank == 0:
    # read the arguments
    args = sys.argv
    clas0, clas = read_clas(args)
    dirname = clas0['dirname']
    tag = clas0['tag']

    # Get the Rayleigh data directory
    radatadir = dirname + '/' + dataname + '/'

    # Get all the file names in datadir and their integer counterparts
    file_list, int_file_list, nfiles = get_file_lists(radatadir, args)

    # Get the problem size
    nproc_min, nproc_max, n_per_proc_min, n_per_proc_max =\
            opt_workload(nfiles, nproc)

    # get grid information
    di_grid = get_grid_info(dirname)
    nt = di_grid['nt']
    nr = di_grid['nr']
    rw = di_grid['rw']
    tw = di_grid['tw']
    rr = di_grid['rr']
    rmin, rmax = np.min(rr)/rsun, np.max(rr)/rsun
    cost = di_grid['cost']
    tt_lat = di_grid['tt_lat']
    latmin, latmax = np.min(tt_lat), np.max(tt_lat)

    ncheby, domain_bounds = get_domain_bounds(dirname)
    domain_bounds = np.array(domain_bounds)/rsun # normalize by rsun
    ndomains = len(ncheby)

    # get grid info + default kwargs
    kwargs_default = dict({})
    kwargs_default['nquadr'] = ndomains # can divide up the radial grid into nquadr equally spaced domains
    kwargs_default['rbounds'] = None # can specify radial domain boundaries directly
    kwargs_default['nquadlat'] = 4 # "high and low" latitudes in both North and South
    kwargs_default['latbounds'] = None

    # update these possibly
    kw = update_dict(kwargs_default, clas)
    if kw.latbounds is None: # this is the default
        kw.latbounds = np.linspace(latmin, latmax, kw.nquadlat + 1)
    if kw.rbounds is None: # this is the default
        kw.rbounds = np.linspace(rmin, rmax, kw.nquadr + 1)

    # convert everything into a sorted array
    latbounds = np.sort(make_array(kw.latbounds))
    rbounds = np.sort(make_array(kw.rbounds))

    # update the number of quadrants
    nquadlat = len(latbounds) - 1
    nquadr = len(rbounds) - 1
    nquad = nquadlat*nquadr

    # get the associated grid indices and update the actual boundary vals
    ilatbounds = inds_from_vals(tt_lat, latbounds)
    irbounds = inds_from_vals(rr/rsun, rbounds)

    latbounds = tt_lat[ilatbounds]
    rbounds = rr[irbounds]/rsun

    # compute the volumes of each quadrant
    volumes = np.zeros((nquadlat, nquadr))
    for ilat in range(nquadlat): # remember: tt_lat is INCREASING
        it1 = ilatbounds[ilat]
        it2 = ilatbounds[ilat + 1]

        for ir in range(nquadr): # remember: rr is DECREASING
            ir2 = irbounds[ir]
            ir1 = irbounds[ir + 1]
            volumes[ilat, ir] = 2.*np.pi/3.*(rr[ir1]**3. - rr[ir2]**3.)*(cost[it2] - cost[it1])

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
dirname, radatadir, ilatbounds, irbounds, tw, rw]
else:
    meta = None
dirname, radatadir, ilatbounds, irbounds, tw, rw = comm.bcast(meta, root=0)

# Checkpoint and time
comm.Barrier()
if rank == 0:
    t2 = time.time()
    print (format_time(t2 - t1))
    print ("tracing over %i x %i = %i quadrants" %(nquadlat, nquadr, nquad))
    print ("rbounds/rsun = " + arr_to_str(rbounds, "%.3f"))
    print ("latbounds = " + arr_to_str(latbounds, "%.1f"))
    print ('Considering %i %s files for the trace: %s through %s'\
        %(nfiles, dataname, file_list[0], file_list[-1]))
    print(fill_str('computing', lent, char), end='\r')
    t1 = time.time()

# Now analyze the data
nquadlat = len(ilatbounds) - 1
nquadr = len(irbounds) - 1
my_times = []
my_iters = []
my_vals = []

for i in range(my_nfiles):
    a = reading_func(radatadir + str(my_files[i]).zfill(8), '')
    for j in range(a.niter):
        vals_loc = a.vals[..., j]

        # Get the values in the separate quadrants
        vals_gav = np.zeros((a.nq, nquadlat, nquadr))
        for ilat in range(nquadlat): # remember: tt_lat is INCREASING
            it1 = ilatbounds[ilat]
            it2 = ilatbounds[ilat + 1]

            for ir in range(nquadr): # remember: rr is DECREASING
                ir2 = irbounds[ir]
                ir1 = irbounds[ir + 1]

                vals_quad = vals_loc[it1:it2+1, ir1:ir2+1, :]
                tw_quad = (tw[it1:it2+1]/np.sum(tw[it1:it2+1])).\
                        reshape((it2 - it1 + 1, 1, 1))
                rw_quad = (rw[ir1:ir2+1]/np.sum(rw[ir1:ir2+1])).\
                        reshape((ir2 - ir1 + 1, 1))

                gav = np.sum(tw_quad*vals_quad, axis=0)
                gav = np.sum(rw_quad*gav, axis=0)
                vals_gav[:, ilat, ir] = gav
        my_vals.append(vals_gav)
        my_times.append(a.time[j])
        my_iters.append(a.iters[j])
    if rank == 0:
        pcnt_done = i/my_nfiles*100.
        print(fill_str('computing', lent, char) +\
                ('rank 0 %5.1f%% done' %pcnt_done), end='\r')
# Checkpoint and time
comm.Barrier()
if rank == 0:
    t2 = time.time()
    print('\n' + fill_str('computing time', lent, char), end='')
    print (format_time(t2 - t1))
    print(fill_str('rank 0 collecting and saving the results',\
            lent, char), end='')
    t1 = time.time()

# proc 0 now collects the results from each process
if rank == 0:
    # initialize empty lists for the final arrays
    times = []
    iters = []
    vals = []

    # Gather the results
    for j in range(nproc):
        if j >= 1:
            # Get my_times, my_iters, my_vals from rank j
            my_times, my_iters, my_vals = comm.recv(source=j)
        times += my_times
        iters += my_iters
        vals += my_vals
    # convert everything to arrays
    vals = np.array(vals)
    times = np.array(times)
    iters = np.array(iters)
else: # other processes send their data
    comm.send([my_times, my_iters, my_vals], dest=0)

# Make sure proc 0 collects all data
comm.Barrier()

# proc 0 saves the data
if rank == 0:
    # create data directory if it doesn't already exist
    datadir = dirname + '/data/'
    if not os.path.isdir(datadir):
        os.makedirs(datadir)

    # Set the timetrace savename by the directory, what we are saving,
    # and first and last iteration files for the trace
    dirname_stripped = strip_dirname(dirname)
    savename = 'G_Avgs_trace_quad' + tag + '-' + file_list[0] + '_' + file_list[-1] + '.pkl'
    savefile = datadir + savename

    # compute the full average over the whole volume (really a sanity check)
    vals_full = np.zeros(np.shape(vals[:, :, 0, 0]))
    for ilat in range(nquadlat):
        for ir in range(nquadr):
            vals_full += volumes[ilat, ir]*vals[:, :, ilat, ir]
    vals_full /= np.sum(volumes)

    # save the data
    # Get first and last iters of files
    iter1, iter2 = int_file_list[0], int_file_list[-1]
    f = open(savefile, 'wb')
    pickle.dump({'vals': vals, 'vals_full': vals_full,'times': times, 'iters': iters, 'lut': a.lut, 'qv': a.qv, 'volumes': volumes, 'rbounds': rbounds, 'latbounds': latbounds}, f, protocol=4)
    f.close()
    t2 = time.time()
    print (format_time(t2 - t1))
    print(make_bold(fill_str('total time', lent, char)), end='')
    print (make_bold(format_time(t2 - t1_glob)))
    print ('data saved at ')
    print (make_bold(savefile))
