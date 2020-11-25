# Routine to trace Rayleigh G_Avgs data in time
# Created by: Loren Matilsky
# Before: 11/25/2020
##################################################################
# This routine computes the trace in time of the values in the Shell_Avgs
# data integrated separately in two domains (2dom; CZ and RZ) 
# for a particular simulation. 

# By default, the routine traces over the last 100 files of datadir, though
# the user can specify a different range in sevaral ways:
# -n 10 (last 10 files)
# -range iter1 iter2 (names of start and stop data files; iter2 can be 'last')
# -centerrange iter0 nfiles (trace about central file iter0 over nfiles)

# Import relevant modules
from mpi4py import MPI
import numpy as np
import pickle
import sys, os
sys.path.append(os.environ['rapp'])
sys.path.append(os.environ['raco'])
from mpi_util import opt_workload
from common import get_file_lists, get_desired_range, fill_str,\
    strip_dirname, get_parameter
from get_eq import get_eq
import time

# Get data type and reading function
from rayleigh_diagnostics import Shell_Avgs, GridInfo

dataname = 'Shell_Avgs'
reading_func = Shell_Avgs

# initialize communication
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
nproc = comm.Get_size()

if rank == 0:
    if nproc > 1:
        print ('processing in parallel with %i ranks' %nproc)
    else:
        print ('processing in serial with 1 rank')

# info for future print messages
lent = 50
char = '.'

# Get the name of the run directory (all procs)
dirname = sys.argv[1]

# Get the Rayleigh data directory (all procs)
radatadir = dirname + '/' + dataname + '/'

# Each proc gets grid information
gi = GridInfo(dirname + '/grid_info', '')
rw = gi.rweights
nr = gi.nr
ir_bcz = get_parameter(dirname, 'ncheby')[1] - 1
nr_cz = ir_bcz + 1
nr_rz = nr - nr_cz

# Each proc gets averaging weights for CZ and RZ separately
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

# Checkpoint and time
if rank == 0:
    for k in range(1, nproc): # check each process has reached this block
        pi_int = np.array([0])
        comm.Recv(pi_int, source=k)
        pi_int = pi_int[0]
    # they're all at this block, start the clock
    t1_glob = time.time()
    t1 = np.copy(t1_glob)
    print(fill_str('proc 0 distributing the file lists', lent, char),\
            end='')
    for k in range(1, nproc): # tell each process to continue
        comm.Send(np.array([pi_int]), dest=k)
else:
    pi_int = 314
    comm.Send(np.array([pi_int]), dest=0)
    pi_int = np.array([0])
    comm.Recv(pi_int, source=0)
    pi_int = pi_int[0]

# proc 0 reads the file lists and distributes them
if rank == 0:
    # Get all the file names in datadir and their integer counterparts
    file_list, int_file_list, nfiles = get_file_lists(radatadir)

    # Read in CLAs
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

    for k in range(1, nproc):
        # distribute the partial file list to other procs 
        if k < nproc_min: # first processes analyze fewer files
            my_nfiles = np.copy(n_per_proc_min)
            istart = k*my_nfiles
            iend = istart + my_nfiles
        else: # last processes analyze more files
            my_nfiles = np.copy(n_per_proc_max)
            istart = nproc_min*n_per_proc_min + (k - nproc_min)*my_nfiles
            iend = istart + my_nfiles

        # Send some meta parameters
        arr = np.array([nproc_min, nproc_max, n_per_proc_min,\
                n_per_proc_max, nq, nrec_full, my_nfiles])
        comm.Send(arr, dest=k)

        # send the file_list
        my_files = np.copy(int_file_list[istart:iend])
        comm.Send(my_files, dest=k)

    # Will need nrec (last niter) to get proper time axis size
    nrec_last = np.array([0])
    if nproc > 1: # can't receive from myself!
        comm.Recv(nrec_last, source=nproc-1)
        nrec_last = nrec_last[0]
    else:
        af = reading_func(radatadir + file_list[-1], '')
        nrec_last = af.niter
        print ('nrec_last = ', nrec_last)

    # Make sure proc 0 also analyzes n_per_proc_min files
    if nproc > 1:
        my_ntimes = n_per_proc_min*nrec_full
    else:
        my_ntimes = (n_per_proc_min - 1)*nrec_full + nrec_last
    my_nfiles = n_per_proc_min
    my_files = int_file_list[:my_nfiles]

else: # other processes receive the data
    # all processes need parameters from proc 0
    arr = np.zeros(7, dtype='int') 
    comm.Recv(arr, source=0)
    nproc_min, nproc_max, n_per_proc_min, n_per_proc_max,\
            nq, nrec_full, my_nfiles = arr

    # Get the partial file list from proc 0 
    my_files = np.zeros(my_nfiles, dtype='int')
    comm.Recv(my_files, source=0)

    if rank == nproc - 1:
        # last process gets nrec_last
        af = reading_func(radatadir + str(my_files[-1]).zfill(8), '')
        nrec_last = af.niter
        comm.Send(np.array([nrec_last]), dest=0)
        my_ntimes = (my_nfiles - 1)*nrec_full + nrec_last
    else:
        my_ntimes = my_nfiles*nrec_full

# Checkpoint and time
if rank == 0:
    for k in range(1, nproc): # check each process has reached this block
        pi_int = np.array([0])
        comm.Recv(pi_int, source=k)
        pi_int = pi_int[0]
    # they're all at this block, restart the clock
    t2 = time.time()
    print ('%8.2e s' %(t2 - t1))
    print ('Considering %i %s files for the trace: %s through %s'\
        %(nfiles, dataname, file_list[0], file_list[-1]))
    print(fill_str('computing', lent, char), end='\r')
    t1 = time.time()
    for k in range(1, nproc): # tell each process to continue
        comm.Send(np.array([pi_int]), dest=k)
else:
    pi_int = 314
    comm.Send(np.array([pi_int]), dest=0)
    pi_int = np.array([0])
    comm.Recv(pi_int, source=0)
    pi_int = pi_int[0]

# Now analyze the data
my_times = np.zeros(my_ntimes)
my_iters = np.zeros(my_ntimes, dtype='int')
my_vals = np.zeros((nq, my_ntimes))
my_vals_cz = np.zeros((nq, my_ntimes))
my_vals_rz = np.zeros((nq, my_ntimes))

my_count = 0
for i in range(my_nfiles):
    if rank == 0:
        print(fill_str('computing', lent, char) + 'rank 0 on ' + dataname +\
                '/' + str(my_files[i]).zfill(8), end='\r')
    if rank == 0 and i == 0:
        a = a0
    elif rank == nproc - 1 and i == my_nfiles - 1:
        a = af
    else:   
        a = reading_func(radatadir + str(my_files[i]).zfill(8), '')
    for j in range(a.niter):
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

        my_vals[:, my_count] = gav
        my_vals_cz[:, my_count] = gav_cz
        my_vals_rz[:, my_count] = gav_rz
        my_times[my_count] = a.time[j] 
        my_iters[my_count] = a.iters[j]
        my_count += 1

if rank == 0:
    print(fill_str('computing', lent, char) + 'rank 0 done      ', end='\r')

# Checkpoint and time
if rank == 0:
    for k in range(1, nproc): # check each process has reached this block
        pi_int = np.array([0])
        comm.Recv(pi_int, source=k)
        pi_int = pi_int[0]
    # they're all at this block, restart the clock
    t2 = time.time()
    print(fill_str('computing', lent, char), end='')
    print ('%8.2e s                                 ' %(t2 - t1))
    print(fill_str('rank 0 collecting and saving the results',\
            lent, char), end='')
    t1 = time.time()
    for k in range(1, nproc): # tell each process to continue
        comm.Send(np.array([pi_int]), dest=k)
else:
    pi_int = 314
    comm.Send(np.array([pi_int]), dest=0)
    pi_int = np.array([0])
    comm.Recv(pi_int, source=0)
    pi_int = pi_int[0]

# proc 0 now collects the results from each process
if rank == 0:
    # Get dimensions of output
    ntimes = (nfiles - 1)*nrec_full + nrec_last

    # Initialize zero-filled 'vals/times/iters' arrays to store the data
    times = np.zeros(ntimes)
    iters = np.zeros(ntimes, dtype='int')
    vals = np.zeros((nq, ntimes))
    vals_cz = np.zeros((nq, ntimes))
    vals_rz = np.zeros((nq, ntimes))

    # First, proc 0 already has some results
    times[:my_ntimes] = my_times
    iters[:my_ntimes] = my_iters
    vals[:, :my_ntimes] = my_vals
    vals_cz[:, :my_ntimes] = my_vals_cz
    vals_rz[:, :my_ntimes] = my_vals_rz
    # now get the results from the other processes
    istart = np.copy(my_ntimes)
    for j in range(1, nproc):
        if j == nproc - 1:
            my_nfiles = n_per_proc_max
            my_ntimes = (my_nfiles - 1)*nrec_full + nrec_last
        elif j < nproc_min:
            my_nfiles = n_per_proc_min
            my_ntimes = my_nfiles*nrec_full
        else:
            my_nfiles = n_per_proc_max
            my_ntimes = my_nfiles*nrec_full
            
        my_times = np.zeros(my_ntimes)
        my_iters = np.zeros(my_ntimes, dtype='int')
        my_vals = np.zeros((nq, my_ntimes))
        my_vals_cz = np.zeros((nq, my_ntimes))
        my_vals_rz = np.zeros((nq, my_ntimes))

        comm.Recv(my_times, source=j)
        comm.Recv(my_iters, source=j)
        comm.Recv(my_vals, source=j)
        comm.Recv(my_vals_cz, source=j)
        comm.Recv(my_vals_rz, source=j)

        times[istart:istart+my_ntimes] = my_times
        iters[istart:istart+my_ntimes] = my_iters
        vals[:, istart:istart+my_ntimes] = my_vals
        vals_cz[:, istart:istart+my_ntimes] = my_vals_cz
        vals_rz[:, istart:istart+my_ntimes] = my_vals_rz

        istart += my_ntimes
else: # other processes send their data
    comm.Send(my_times, dest=0)
    comm.Send(my_iters, dest=0)
    comm.Send(my_vals, dest=0)
    comm.Send(my_vals_cz, dest=0)
    comm.Send(my_vals_rz, dest=0)

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
    pickle.dump({'vals': vals.T, 'vals_cz': vals_cz.T, 'vals_rz': vals_rz.T, 'times': times, 'iters': iters, 'lut': lut, 'ntimes': ntimes, 'iter1': iter1, 'iter2': iter2, 'rr': a0.radius, 'nr': a0.nr, 'qv': qv, 'nq': nq}, f, protocol=4)
    f.close()
    t2 = time.time()
    print ('%8.2e s' %(t2 - t1))
    print(fill_str('total time', lent, char), end='')
    print ('%8.2e s' %(t2 - t1_glob))
