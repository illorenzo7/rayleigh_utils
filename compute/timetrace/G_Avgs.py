# Routine to trace Rayleigh G_Avgs data in time
# Created by: Loren Matilsky
# On: 12/18/2018
##################################################################
# This routine computes the trace in time of the values in the G_Avgs data 
# for a particular simulation. 

# By default, the routine traces over the last 100 files of datadir, though
# the user can specify a different range in sevaral ways:
# -n 10 (last 10 files)
# -range iter1 iter2 (names of start and stop data files; iter2 can be "last")
# -centerrange iter0 nfiles (trace about central file iter0 over nfiles)

# Import relevant modules
from mpi4py import MPI
import numpy as np
import pickle
import sys, os
sys.path.append(os.environ['rapp'])
sys.path.append(os.environ['raco'])
from mpi_util import opt_workload
from rayleigh_diagnostics import G_Avgs
from common import get_file_lists, get_desired_range, strip_dirname,\
        fill_str
import time

# initialize communication
comm = MPI.COMM_WORLD
rank = comm.Get_rank()

# Get the name of the run directory
dirname = sys.argv[1]
# Get the stripped name to use in file naming
dirname_stripped = strip_dirname(dirname)

# Get the Rayleigh data directory
radatadir = dirname + '/G_Avgs/'

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

# Get Parallel problem size
nproc = comm.Get_size()
nproc_min, nproc_max, n_per_proc_min, n_per_proc_max =\
        opt_workload(nfiles, nproc)

# Find the relevant place to store the data, and create the directory if it
# doesn't already exist
if rank == 0:
    # start timing
    t1 = time.time()
    # parameters for timing messages
    lent = 40
    char = "."

    datadir = dirname + '/data/'
    if not os.path.isdir(datadir):
        os.makedirs(datadir)

    # Set the timetrace savename by the directory, what we are saving,
    # and first and last iteration files for the trace
    savename = dirname_stripped + '_trace_G_Avgs_' +\
            file_list[0] + '_' + file_list[-1] + '.pkl'
    savefile = datadir + savename

    # Will need the first G_Avgs file for a number of things
    g0 = G_Avgs(radatadir + file_list[0], '')

    # Will need the last G_Avgs file to get the size of the arrays
    gf = G_Avgs(radatadir + file_list[-1], '')
    nrec_full = g0.niter
    nrec_last = gf.niter
    ntimes = (nfiles - 1)*nrec_full + nrec_last
    nq = g0.nq

    # State what we will analyze
    print ('Considering %i G_Avgs files for the trace: %s through %s'\
        %(ntasks, file_list[0], file_list[-1]))

    # will also need first and last iteration files
    iter1, iter2 = int_file_list[0], int_file_list[-1]

    # Initialize empty "vals/times/iters" arrays to store the data
    vals = np.zeros((nq, ntimes))
    times = np.zeros(ntimes)
    iters = np.zeros(ntimes, dtype='int')

    # Distribute the file list (tasks) to other processes
    for k in range(1, nproc):
        # Send smaller portion of file list to other processes
        if k < nproc_min: # first processes analyze fewer files
            istart = k*n_per_proc_min
            iend = istart + n_per_proc_min
        else: # last processes analyze more files
            istart = nproc_min*n_per_proc_min +\
                    (k - nproc_min)*n_per_proc_max
            iend = istart + n_per_proc_max
        my_files = int_file_list[istart:iend]
        comm.Send(my_files, dest=k)
        comm.Send(np.array([nq]), dest=k)
        comm.Send(np.array([nrec_full]), dest=k)
        comm.Send(np.array([nrec_last]), dest=k)

    # Make sure proc 0 also analyzes n_per_proc_min files
    my_ntimes = n_per_proc_min*nrec_full
    my_nfiles = n_per_proc_min
    my_files = int_file_list[:my_nfiles]
    my_times = np.zeros(my_ntimes)
    my_iters = np.zeros(my_ntimes, dtype='int')
    my_vals = np.zeros((nq, my_ntimes))

    for i in range(my_nfiles):
        if i == 0:
            g = g0
        else:   
            g = G_Avgs(radatadir + str(my_files[i]).zfill(8), '')
        my_count = 0
        for j in range(nrec_full):
            my_vals[:, my_count] = g.vals[j, :]
            my_times[my_count] = g.time[j] 
            my_iters[my_count] = g.iters[j]
            my_count += 1

else: # other processes receive the data
    # all processes need nrec_full and nq from proc 0
    nrec_full = np.array([0])
    comm.Recv(nrec_full, source=0)
    nrec_full = nrec_full[0]

    nq = np.array([0])
    comm.Recv(nq, source=0)
    nq = nq[0]

    if rank == nproc - 1:
        # last process gets nrec_last
        gf = G_Avgs(radatadir + file_list[-1], '')
        nrec_last = np.array([gf.niter])
        comm.Send(nrec_last, dest=0)
        nrec_last = nrec_last[0]

        # get length of file list and total number of slices from
        # proc 0 
        my_nfiles = n_per_proc_max
        my_ntimes = (my_nfiles - 1)*nrec_full + nrec_last
    elif rank < nproc_min:
        my_nfiles = n_per_proc_min
        my_ntimes = my_nfiles*nrec_full
    else:
        my_nfiles = n_per_proc_max
        my_ntimes = my_nfiles*nrec_full

    # Get the partial file list from proc 0 
    my_files = np.empty(my_nfiles, dtype='int')
    comm.Recv(my_files, source=0)
   
    # inialize local vals array with proper size
    my_times = np.zeros(my_ntimes)
    my_iters = np.zeros(my_ntimes, dtype='int')
    my_vals = np.zeros((nq, my_ntimes))

    for i in range(my_nfiles):
        if rank == nproc - 1 and i == my_nfiles - 1:
            g = gf
        else:   
            g = G_Avgs(radatadir + str(my_files[i]).zfill(8), '')
        my_count = 0
        for j in range(g.niter):
            my_vals[:, my_count] = g.vals[j, :]
            my_times[my_count] = g.time[j] 
            my_iters[my_count] = g.iters[j]
            my_count += 1

# Now collect the results from each process
if rank == 0:
    # First, proc 0 already has some results
    vals[nq, :my_ntimes] = my_vals
    # now get the results from the other processes
    istart = np.copy(my_ntimes)
    for j in range(1, nproc):
        if j == nproc - 1:
            my_nfiles = n_per_proc_max
            my_ntimes = (my_nfiles - 1)*nrec_full + nrec_last
        elif rank < nproc_min:
            my_nfiles = n_per_proc_min
            my_ntimes = my_nfiles*nrec_full
        else:
            my_nfiles = n_per_proc_max
            my_ntimes = my_nfiles*nrec_full
        my_vals = np.zeros((nq, my_ntimes))
        comm.Recv(my_vals, source=j)
        vals[nq, istart:istart + my_ntimes] = my_vals
        istart += my_ntimes
    else: # other processes receive the data
        comm.Send(my_vals, dest=0)
        comm.Send(my_times, dest=0)
        comm.Send(my_iters, dest=0)



print ('Traced over %i G_Avgs slice(s)' %count)

# Save the avarage
print ('Saving file at ' + savefile)
f = open(savefile, 'wb')
pickle.dump({'vals': vals, 'times': times, 'iters': iters, 'lut': g0.lut, 'count': count, 'iter1': iter1, 'iter2': iter2, 'qv': g0.qv, 'nq': g0.nq},\
        f, protocol=4)
f.close()
t2 = time.time()
print(fill_str("total time", lent, char), end="")
print ("%8.2e s" %(t2 - t1))
