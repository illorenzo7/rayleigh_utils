# Routine to trace Rayleigh G_Avgs data in time
# Created by: Loren Matilsky
# On: 12/18/2018
# Parallelized: 11/26/2020
##################################################################
# This routine computes the trace in time of the values in the G_Avgs data 
# for a particular simulation. 

# By default, the routine traces over the last 100 files of datadir, though
# the user can specify a different range in sevaral ways:
# -n 10 (last 10 files)
# -range iter1 iter2 (names of start and stop data files; iter2 can be 'last')
# -centerrange iter0 nfiles (trace about central file iter0 over nfiles)

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
from rayleigh_diagnostics import Shell_Avgs
reading_func = Shell_Avgs
dataname = 'Shell_Avgs'

if rank == 0:
    # modules needed only by proc 0 
    import pickle
    from get_eq import get_eq

# Checkpoint and time
comm.Barrier()
if rank == 0:
    t2 = time.time()
    print ('%8.2e s' %(t2 - t1))
    print(fill_str('proc 0 distributing the file lists', lent, char),\
            end='')
    t1 = time.time()

# proc 0 reads the file lists and distributes them
if rank == 0:
    # Get the name of the run directory
    dirname = sys.argv[1]

    # Get the Rayleigh data directory
    radatadir = dirname + '/' + dataname + '/'

    # Get all the file names in datadir and their integer counterparts
    file_list, int_file_list, nfiles = get_file_lists(radatadir)

    # Get desired file list from command-line arguments
    args = sys.argv[2:]
    nargs = len(args)
    if nargs == 0:
        index_first, index_last = nfiles - 101, nfiles - 1  
        # By default trace over the last 100 files
    else:
        index_first, index_last = get_desired_range(int_file_list, args)

    # Remove parts of file lists we don't need
    file_list = file_list[index_first:index_last + 1]
    int_file_list = int_file_list[index_first:index_last + 1]
    nfiles = index_last - index_first + 1

    # Get the problem size
    nproc_min, nproc_max, n_per_proc_min, n_per_proc_max =\
            opt_workload(nfiles, nproc)

    # Will need the first data file for a number of things
    a0 = reading_func(radatadir + file_list[0], '')
    nrec_full = a0.niter
    nq = a0.nq + 3 # will add the three internal energies after

    # Will need nrec (last niter) to get proper time axis size
    af = reading_func(radatadir + file_list[-1], '')
    nrec_last = af.niter
    ntimes = (nfiles - 1)*nrec_full + nrec_last

    # get grid information
    gi = GridInfo(dirname + '/grid_info', '')
    rw = gi.rweights
    nr = gi.nr
    ir_bcz = get_parameter(dirname, 'ncheby')[1] - 1
    nr_cz = ir_bcz + 1
    nr_rz = nr - nr_cz

    # get averaging weights for CZ and RZ separately
    rw_cz = np.copy(rw[:nr_cz])
    rw_rz = np.copy(rw[nr_cz:])
    rw_cz /= np.sum(rw_cz)
    rw_rz /= np.sum(rw_rz)
    rw = rw.reshape((nr, 1))
    rw_cz = rw_cz.reshape((nr_cz, 1))
    rw_rz = rw_rz.reshape((nr_rz, 1))

    # Get rho*T
    eq = get_eq(dirname)
    rhot = (eq.density*eq.temperature).reshape((1, nr))

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

# Broadcast dirname, radatadir, nq
if rank == 0:
    meta = [dirname, radatadir, nq, nr, nr_cz, nr_rz, rw, rw_cz, rw_rz,\
            ir_bcz, rhot]
else:
    meta = None
dirname, radatadir, nq, nr, nr_cz, nr_rz, rw, rw_cz, rw_rz, ir_bcz, rhot =\
    comm.bcast(meta, root=0)

# Checkpoint and time
comm.Barrier()
if rank == 0:
    t2 = time.time()
    print ('%8.2e s' %(t2 - t1))
    print ('Considering %i %s files for the trace: %s through %s'\
        %(nfiles, dataname, file_list[0], file_list[-1]))
    print ("ntimes for trace = %i" %ntimes)
    print(fill_str('computing', lent, char), end='\r')
    t1 = time.time()

# Now analyze the data
my_times = np.zeros(my_ntimes)
my_iters = np.zeros(my_ntimes, dtype='int')
my_vals = np.zeros((my_ntimes, nq))
my_vals_cz = np.zeros((my_ntimes, nq))
my_vals_rz = np.zeros((my_ntimes, nq))

my_count = 0
for i in range(my_nfiles):
    if rank == 0 and i == 0:
        a = a0
    else:   
        a = reading_func(radatadir + str(my_files[i]).zfill(8), '')
    for j in range(a.niter):
        if my_count < my_ntimes: # make sure we don't go over the allotted
            # space in the arrays
            vals_loc = np.copy(a.vals[:, 0, :, j])
            # add in internal energy
            inte_loc = rhot*vals_loc[:, a.lut[501]]
            # top S subtracted
            inte_loc_subt = rhot*(vals_loc[:, a.lut[501]] -\
                    vals_loc[0, a.lut[501]])
            # bottom S subtracted
            inte_loc_subb = rhot*(vals_loc[:, a.lut[501]] -\
                    vals_loc[-1, a.lut[501]])

            # add in the three energies
            vals_loc = np.hstack((vals_loc, inte_loc.T, inte_loc_subt.T,\
                    inte_loc_subb.T))

            # Get the values in the CZ/RZ separately
            vals_cz_loc = vals_loc[:ir_bcz + 1]
            vals_rz_loc = vals_loc[ir_bcz + 1:]

            gav = np.sum(rw*vals_loc, axis=0)
            gav_cz = np.sum(rw_cz*vals_cz_loc, axis=0)
            gav_rz = np.sum(rw_rz*vals_rz_loc, axis=0)

            my_vals[my_count, :] = gav
            my_vals_cz[my_count, :] = gav_cz
            my_vals_rz[my_count, :] = gav_rz
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
    print(fill_str('\ncomputing time', lent, char), end='')
    print ('%8.2e s                                 ' %(t2 - t1))
    print(fill_str('rank 0 collecting and saving the results',\
            lent, char), end='')
    t1 = time.time()

# proc 0 now collects the results from each process
if rank == 0:
    # Initialize zero-filled 'vals/times/iters' arrays to store the data
    times = np.zeros(ntimes)
    iters = np.zeros(ntimes, dtype='int')
    vals = np.zeros((ntimes, nq))
    vals_cz = np.zeros((ntimes, nq))
    vals_rz = np.zeros((ntimes, nq))

    # Gather the results into these "master" arrays
    istart = 0
    for j in range(nproc):
        if j >= 1:
            # Get my_ntimes, my_times, my_iters, my_vals from rank j
            my_ntimes, my_times, my_iters,\
                    my_vals, my_vals_cz, my_vals_rz = comm.recv(source=j)
        times[istart:istart+my_ntimes] = my_times
        iters[istart:istart+my_ntimes] = my_iters
        vals[istart:istart+my_ntimes, :] = my_vals
        vals_cz[istart:istart+my_ntimes, :] = my_vals_cz
        vals_rz[istart:istart+my_ntimes, :] = my_vals_rz
        istart += my_ntimes
else: # other processes send their data
    comm.send([my_ntimes, my_times, my_iters,\
            my_vals, my_vals_cz, my_vals_rz], dest=0)

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
    savename = dirname_stripped + '_trace_2dom_G_Avgs_' +\
            file_list[0] + '_' + file_list[-1] + '.pkl'
    savefile = datadir + savename

    # save the data
    # Get first and last iters of files
    iter1, iter2 = int_file_list[0], int_file_list[-1]
    # append the lut (making the inte, inte_subt, and inte_subb quantities
    # (4000, 4001, 4002)
    lut_app = np.array([a0.nq, a0.nq + 1, a0.nq + 2])
    lut = np.hstack((a0.lut, lut_app))
    qv_app = np.array([4000, 4001, 4002])
    qv = np.hstack((a0.qv, qv_app))
    f = open(savefile, 'wb')
    pickle.dump({'vals': vals, 'vals_cz': vals_cz, 'vals_rz': vals_rz, 'times': times, 'iters': iters, 'lut': lut, 'ntimes': ntimes, 'iter1': iter1, 'iter2': iter2, 'rr': a0.radius, 'nr': a0.nr, 'qv': qv, 'nq': nq}, f, protocol=4)
    f.close()
    t2 = time.time()
    print ('%8.2e s' %(t2 - t1))
    print(fill_str('total time', lent, char), end='')
    print ('%8.2e s' %(t2 - t1_glob))
