# Routine to trace Rayleigh Shell_Slices data in time
# Created by: Loren Matilsky
# On: 08/15/2019
# Parallelized: 04/01/2021
############################################################################
# This routine computes the trace in time/longitude of quantities in the 
# Shell_Slices data for a particular simulation. 
#
# By default, the 8 variables are computed at each time, in a latitude strip
# defined by (clat, dlat) = ([central latitude for average],
# [range of latitudes to average over]) at the depths the shell slices were 
# sampled. 
#
# The strip range can be changed using the options --clat and --dlat, e.g., 
# --clat 60 --dlat 30, for a strip averaged between 45 and 75 degrees (North)
#
# By default, the routine traces over all Shell_Slices in the directory,
# though user can specify an alternate range, by, e.g.,
# -n 10 (last 10 files)
# -range iter1 iter2 (no.s for start/stop data files; iter2 can be "last")
# -centerrange iter0 nfiles (trace about central file iter0 over nfiles)
#
# The final datacube output ('vals') will have shape
# (nphi, niter, nr, nq), where nr and nq are the attributes of the shell slices

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
from rayleigh_diagnostics import Shell_Slices
reading_func = Shell_Slices
dataname = 'Shell_Slices'

if rank == 0:
    # modules needed only by proc 0 
    import pickle
    from rayleigh_diagnostics import GridInfo

# Checkpoint and time
comm.Barrier()
if rank == 0:
    t2 = time.time()
    print (format_time(t2 - t1))
    t1 = time.time()

# proc 0 reads the file lists and distributes them, also the meta data
if rank == 0:
    # Get the name of the run directory
    dirname = sys.argv[1]

    # Get the Rayleigh data directory
    radatadir = dirname + '/' + dataname + '/'

    # Get all the file names in datadir and their integer counterparts
    file_list, int_file_list, nfiles = get_file_lists(radatadir)

    # Read in CLAs
    args = sys.argv[2:]
    nargs = len(args)

    # get desired analysis range
    the_tuple = get_desired_range(int_file_list, args)
    if the_tuple is None:
        index_first, index_last = nfiles - 101, nfiles - 1  
        # By default trace over the last 100 files
    else:
        index_first, index_last = the_tuple

    # Remove parts of file lists we don't need
    file_list = file_list[index_first:index_last + 1]
    int_file_list = int_file_list[index_first:index_last + 1]
    nfiles = index_last - index_first + 1

    # Read in grid info from GridInfo file
    gi = GridInfo(dirname + '/grid_info', '')
    rr = gi.radius
    ri, ro = np.min(rr), np.max(rr)
    d = ro - ri
    rr_depth = (ro - rr)/d
    rr_height = (rr - ri)/d
    sint, cost = gi.sintheta, gi.costheta
    tt = np.arccos(cost)
    tt_lat = (np.pi/2 - tt)*180./np.pi
    nt, nr = len(sint), len(rr)
    nphi = 2*nt
    lons = np.arange(0., 360., 360/nphi)
    tw = gi.tweights

    # By default trace over all quantities and depths
    a0 = reading_func(radatadir + file_list[0], '')
    qvals = a0.qv
    all_depths = (ro - a0.radius)/(ro - ri)
    depths = np.copy(all_depths)
    aspect = 0.

    # Set other defaults
    tag = ''
    clat = 10.
    dlat = 0. # by default do NOT average over latitude

    # reset parameters with CLAs
    for i in range(nargs):
        arg = args[i]
        if arg == '--clat':
            clat = float(args[i+1])
        if arg == '--dlat':
            dlat = float(args[i+1])

    # convert things to arrays
    depths = np.array(depths)
    qvals = np.array(qvals)

    # didn't put this earlier because it got messed up by the message
    # saying which depths we were taking
    print(fill_str('proc 0 distributing the file lists', lent, char),\
            end='')
 
    # Desired latitude range
    ith1 = np.argmin(np.abs(tt_lat - (clat - dlat/2.)))
    ith2 = np.argmin(np.abs(tt_lat - (clat + dlat/2.)))
    lats_strip = tt_lat[ith1:ith2+1]
    tw_strip = tw[ith1:ith2+1]
    tw_strip /= np.sum(tw_strip)
    tw_strip_1d = np.copy(tw_strip)
    tw_strip = tw_strip.reshape((1, ith2 - ith1 + 1, 1, 1))

    # Get the problem size
    nproc_min, nproc_max, n_per_proc_min, n_per_proc_max =\
            opt_workload(nfiles, nproc)

    # Get metadata
    nrec_full = a0.niter
    ndepths = len(depths)
    nq = len(qvals)

    # get r-indices associated with depths
    rinds = []
    rinds_sslice = []
    for depth in depths:
        rinds.append(np.argmin(np.abs(rr_depth - depth)))
        rinds_sslice.append(np.argmin(np.abs(all_depths - depth)))
    rinds = np.array(rinds)
    rinds_sslice = np.array(rinds_sslice)

    # recompute the actual depths we get
    depths = rr_depth[rinds]
    rvalscm = ro - d*depths
    rvals = rvalscm/rsun

    # Will need nrec (last niter) to get proper time axis size
    af = reading_func(radatadir + file_list[-1], '')
    nrec_last = af.niter
    ntimes = (nfiles - 1)*nrec_full + nrec_last

    # Distribute file_list and my_ntimes to each process
    for k in range(nproc - 1, -1, -1):
        # distribute the partial file list to other procs 
        if k >= nproc_min: # last processes analyzes more files
            my_nfiles = np.copy(n_per_proc_max)
            istart = nproc_min*n_per_proc_min + (k - nproc_min)*my_nfiles
            iend = istart + my_nfiles
        else: # first processes analyze fewer files
            my_nfiles = np.copy(n_per_proc_min)
            istart = k*my_nfiles
            iend = istart + my_nfiles

        if k == nproc - 1: # last process may have nrec_last != nrec_full
            my_ntimes = (my_nfiles - 1)*nrec_full + nrec_last
        else:
            my_ntimes = my_nfiles*nrec_full

        # Get the file list portion for rank k
        my_files = np.copy(int_file_list[istart:iend])

        # send  my_files, my_nfiles, my_ntimes if nproc > 1
        if k >= 1:
            comm.send([my_files, my_nfiles, my_ntimes], dest=k)
else: # recieve my_files, my_nfiles, my_ntimes
    my_files, my_nfiles, my_ntimes = comm.recv(source=0)

# Broadcast meta data
if rank == 0:
    meta = [dirname, radatadir, qvals, nq, rinds_sslice, ndepths, nphi,\
            ith1, ith2, tw_strip]
else:
    meta = None
dirname, radatadir, qvals, nq, rinds_sslice, ndepths, nphi,\
        ith1, ith2, tw_strip = comm.bcast(meta, root=0)

# Checkpoint and time
comm.Barrier()
if rank == 0:
    t2 = time.time()
    print (format_time(t2 - t1))
    print ('Considering %i %s files for the trace: %s through %s'\
        %(nfiles, dataname, file_list[0], file_list[-1]))
    print ("ntimes for trace = %i" %ntimes)
    print ("central latitude = %.1f" %clat)
    if ith2 - ith1 + 1 > 1:
        print("averaging over %.2f in latitude" %dlat)
    else:
        print("not averaging in latitude")
    print(fill_str('computing', lent, char), end='\r')
    t1 = time.time()

# Now analyze the data
my_times = np.zeros(my_ntimes)
my_iters = np.zeros(my_ntimes, dtype='int')
my_vals = np.zeros((my_ntimes, nphi, ndepths, nq))

my_count = 0
for i in range(my_nfiles):
    if rank == 0 and i == 0:
        a = a0
    else:   
        a = reading_func(radatadir + str(my_files[i]).zfill(8), '')
    for j in range(a.niter):
        if my_count < my_ntimes: # make sure we don't go over the allotted
            # space in the arrays
            # get desired values in the strip, careful not to change 
            # the dimensionality of the array: (nphi, nlats, ndepths, nq)
            vals_strip = a.vals[:, ith1:ith2+1, :, :, j]
            vals_strip = vals_strip[:, :, rinds_sslice, :]
            vals_strip = vals_strip[:, :, :, a.lut[qvals]]
            my_vals[my_count, :, :, :] = np.sum(vals_strip*tw_strip, axis=1)
            my_times[my_count] = a.time[j] 
            my_iters[my_count] = a.iters[j]
        my_count += 1
    if rank == 0:
        pcnt_done = my_count/my_ntimes*100.
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
    # Initialize zero-filled 'vals/times/iters' arrays to store the data
    vals = np.zeros((ntimes, nphi, ndepths, nq))
    times = np.zeros(ntimes)
    iters = np.zeros(ntimes, dtype='int')

    # Gather the results into these "master" arrays
    istart = 0
    for j in range(nproc):
        if j >= 1:
            # Get my_ntimes, my_times, my_iters, my_vals from rank j
            my_ntimes, my_times, my_iters, my_vals = comm.recv(source=j)
        times[istart:istart+my_ntimes] = my_times
        iters[istart:istart+my_ntimes] = my_iters
        vals[istart:istart+my_ntimes, :, :, :] = my_vals
        istart += my_ntimes
else: # other processes send their data
    comm.send([my_ntimes, my_times, my_iters, my_vals], dest=0)

# Make sure proc 0 collects all data
comm.Barrier()

# proc 0 saves the data
if rank == 0:
    # create data directory if it doesn't already exist
    datadir = dirname + '/data/'
    if not os.path.isdir(datadir):
        os.makedirs(datadir)

    # Set the timetrace savename by the directory, what we are saving, 
    # and first and last iteration files for the trace (and optional tag)
    if clat >= 0.:
        hemisphere = 'N'
    else:
        hemisphere = 'S'
        
    dirname_stripped = strip_dirname(dirname)
    savename = dirname_stripped + '_time-longitude_' + tag +\
            ('clat%s%02.0f_dlat%03.0f' %(hemisphere, np.abs(clat), dlat)) +\
            '_' + file_list[0] + '_' + file_list[-1] + '.pkl'
    savefile = datadir + savename    

    # save the data
    f = open(savefile, 'wb')
    # will also need first and last iteration files
    iter1, iter2 = int_file_list[0], int_file_list[-1]
    pickle.dump({'vals': vals, 'times': times, 'iters': iters, 'qvals': a0.qv,\
            'rinds': rinds, 'rinds_sslice': rinds_sslice, 'qvals': qvals,\
        'rvals': rvals, 'nrvals': len(rvals), 'lut': a0.lut,\
       'niter': len(iters), 'nq': nq, 'iter1': iter1, 'iter2': iter2, 'rr': rr,\
        'nr': nr, 'ri': ri, 'ro': ro, 'tt': tt, 'tt_lat': tt_lat, 'sint': sint,\
        'cost': cost,'nt': nt, 'lats_strip': lats_strip, 'clat': clat,\
        'dlat': dlat, 'lons':lons, 'nphi': nphi, 'tw': tw,\
        'tw_strip': tw_strip_1d}, f, protocol=4)
    f.close()
    t2 = time.time()
    print (format_time(t2 - t1))
    print(make_bold(fill_str('total time', lent, char)), end='')
    print (make_bold(format_time(t2 - t1_glob)))
    print ('data saved at ')
    print (make_bold(savefile))
