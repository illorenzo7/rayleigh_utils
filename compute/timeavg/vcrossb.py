##################################################################
# Routine to average Rayleigh data in time (generic)
# Author: Loren Matilsky
# Created: 04/08/2021
##################################################################
# This routine computes the time average of v X B components from
# Rayleigh data (AZ_Avgs, Meridional Slices), breaking up terms 
# by Reynolds decomposition
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
from rayleigh_diagnostics import AZ_Avgs, Meridional_Slices
reading_func1 = AZ_Avgs
reading_func2 = Meridional_Slices
dataname1 = 'AZ_Avgs'
dataname2 = 'Meridional_Slices'

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
    # Get the name of the run directory
    dirname = sys.argv[1]
    args = sys.argv[2:]

    # Get the Rayleigh data directory
    radatadir1 = dirname + '/' + dataname1 + '/'
    radatadir2 = dirname + '/' + dataname2 + '/'

    # Get all the file names in datadir and their integer counterparts
    # Use mer slices, since there are sometimes fewer of them than azavgs
    file_list, int_file_list, nfiles = get_file_lists(radatadir1, args)

    # Get the problem size
    nproc_min, nproc_max, n_per_proc_min, n_per_proc_max =\
            opt_workload(nfiles, nproc)

    # Will need the first data file for a number of things
    a0 = reading_func1(radatadir1 + file_list[0], '')
    # Make sure you put same number of slices in AZ_Avgs and 
    # Meridional_Slices, or this won't work

    # get grid information
    di_grid = get_grid_info(dirname)
    nt = di_grid['nt']
    nr = di_grid['nr']
    tt = di_grid['tt']
    rr = di_grid['rr']
    rr_3d = di_grid['rr_3d']
    cost_3d = di_grid['cost_3d']
    sint_3d = di_grid['sint_3d']

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

        # Get the file list portion for rank k
        my_files = np.copy(int_file_list[istart:iend])

        # send  my_files, my_nfiles, my_ntimes if nproc > 1
        if k >= 1:
            comm.send([my_files, my_nfiles], dest=k)
else: # recieve my_files, my_nfiles, my_ntimes
    my_files, my_nfiles = comm.recv(source=0)

# Broadcast dirname, radatadirs, rr, etc.
if rank == 0:
    meta = [dirname, radatadir1, radatadir2, tt, rr, nt, nr, rr_3d, cost_3d, sint_3d, nfiles]
else:
    meta = None
dirname, radatadir1, radatadir2, tt, rr, nt, nr, rr_3d, cost_3d, sint_3d, nfiles = comm.bcast(meta, root=0)

# Checkpoint and time
comm.Barrier()
if rank == 0:
    t2 = time.time()
    print (format_time(t2 - t1))
    print ('Considering %i %s and %s files for the average: %s through %s'\
        %(nfiles, dataname1, dataname2, file_list[0], file_list[-1]))
    print ("no. slices = %i" %nfiles)
    print(fill_str('computing', lent, char), end='\r')
    t1 = time.time()

# Now analyze the data
nq = 9*3
my_vals = np.zeros((nt, nr, nq))

for i in range(my_nfiles):
    a = reading_func1(radatadir1 + str(my_files[i]).zfill(8), '')
    mer = reading_func2(radatadir2 + str(my_files[i]).zfill(8), '')

    my_weight = 1.0/(nfiles*a.niter)
    for j in range(a.niter):
        # full v
        vr = mer.vals[:, :, :, mer.lut[1], j]
        vt = mer.vals[:, :, :, mer.lut[2], j]
        vp = mer.vals[:, :, :, mer.lut[3], j]

        # full B
        br = mer.vals[:, :, :, mer.lut[801], j]
        bt = mer.vals[:, :, :, mer.lut[802], j]
        bp = mer.vals[:, :, :, mer.lut[803], j]

        # mean v
        vr_m = a.vals[:, :, a.lut[1], j].reshape((1, nt, nr))
        vt_m = a.vals[:, :, a.lut[2], j].reshape((1, nt, nr))
        vp_m = a.vals[:, :, a.lut[3], j].reshape((1, nt, nr))

        # mean B
        br_m = a.vals[:, :, a.lut[801], j].reshape((1, nt, nr))
        bt_m = a.vals[:, :, a.lut[802], j].reshape((1, nt, nr))
        bp_m = a.vals[:, :, a.lut[803], j].reshape((1, nt, nr))

        # compute v X B

        # full (r, theta, phi)
        ind_off_full = 0
        my_vals[:, :, ind_off_full + 0] += np.mean(vt*bp - vp*bt, axis=0)*my_weight
        my_vals[:, :, ind_off_full + 1] += np.mean(vp*br - vr*bp, axis=0)*my_weight
        my_vals[:, :, ind_off_full + 2] += np.mean(vr*bt - vt*br, axis=0)*my_weight

        my_vals[:, :, ind_off_full + 3] += np.mean(vt*bp, axis=0)*my_weight
        my_vals[:, :, ind_off_full + 4] += np.mean(-vp*bt, axis=0)*my_weight
        my_vals[:, :, ind_off_full + 5] += np.mean(vp*br, axis=0)*my_weight
        my_vals[:, :, ind_off_full + 6] += np.mean(-vr*bp, axis=0)*my_weight
        my_vals[:, :, ind_off_full + 7] += np.mean(vr*bt, axis=0)*my_weight
        my_vals[:, :, ind_off_full + 8] += np.mean(-vt*br, axis=0)*my_weight

        # mean (r, theta, phi)
        ind_off_mean = 9
        my_vals[:, :, ind_off_mean + 0] += np.mean(vt_m*bp_m - vp_m*bt_m, axis=0)*my_weight
        my_vals[:, :, ind_off_mean + 1] += np.mean(vp_m*br_m - vr_m*bp_m, axis=0)*my_weight
        my_vals[:, :, ind_off_mean + 2] += np.mean(vr_m*bt_m - vt_m*br_m, axis=0)*my_weight

        my_vals[:, :, ind_off_mean + 3] += np.mean(vt_m*bp_m, axis=0)*my_weight
        my_vals[:, :, ind_off_mean + 4] += np.mean(-vp_m*bt_m, axis=0)*my_weight
        my_vals[:, :, ind_off_mean + 5] += np.mean(vp_m*br_m, axis=0)*my_weight
        my_vals[:, :, ind_off_mean + 6] += np.mean(-vr_m*bp_m, axis=0)*my_weight
        my_vals[:, :, ind_off_mean + 7] += np.mean(vr_m*bt_m, axis=0)*my_weight
        my_vals[:, :, ind_off_mean + 8] += np.mean(-vt_m*br_m, axis=0)*my_weight

        # fluc terms
        ind_off_fluc = 18
        for k in range(9):
            my_vals[:, :, ind_off_fluc + k] = my_vals[:, :, ind_off_full + k] -\
                my_vals[:, :, ind_off_mean + k]

        if rank == 0:
            pcnt_done = (i + 1)/my_nfiles*100.
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
    vals = np.zeros((nt, nr, 18))

    # Gather the results into this "master" array
    for j in range(nproc):
        if j >= 1:
            # Get my_ntimes, my_times, my_iters, my_vals from rank j
            my_vals = comm.recv(source=j)
    
else: # other processes send their data
    comm.send(my_vals, dest=0)

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
    savename = 'vcrossb-' +\
            file_list[0] + '_' + file_list[-1] + '.pkl'
    savefile = datadir + savename

    # save the data
    # Get first and last iters of files
    iter1, iter2 = int_file_list[0], int_file_list[-1]
    f = open(savefile, 'wb')
    pickle.dump({'vals': vals}, f, protocol=4)
    f.close()
    t2 = time.time()
    print (format_time(t2 - t1))
    print(make_bold(fill_str('total time', lent, char)), end='')
    print (make_bold(format_time(t2 - t1_glob)))
    print ('data saved at ')
    print (make_bold(savefile))
