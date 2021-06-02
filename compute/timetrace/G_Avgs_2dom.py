##################################################################
# Routine to trace Rayleigh G_Avgs data in time
# Author: Loren Matilsky
# Created: 12/18/2018
# Parallelized: 11/26/2020
##################################################################
# This routine computes the trace in time of the values in the G_Avgs data 
# for a particular simulation. Computes global average for:
# the whole shell (vals)
# the CZ (vals_cz)
# the RZ (vals_rz)
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
from rayleigh_diagnostics import Shell_Avgs
reading_func = Shell_Avgs
dataname = 'Shell_Avgs'

if rank == 0:
    # modules needed only by proc 0 
    import pickle
    from rayleigh_diagnostics import GridInfo

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
    dirname = sys.argv[1]
    args = sys.argv[2:]

    # Get the Rayleigh data directory
    radatadir = dirname + '/' + dataname + '/'

    # Get all the file names in datadir and their integer counterparts
    file_list, int_file_list, nfiles = get_file_lists(radatadir, args)

    # Get the problem size
    nproc_min, nproc_max, n_per_proc_min, n_per_proc_max =\
            opt_workload(nfiles, nproc)

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

        # Get the file list portion for rank k
        my_files = np.copy(int_file_list[istart:iend])

        # send appropriate file info if nproc > 1
        if k >= 1:
            comm.send([my_files, my_nfiles], dest=k)
else: # recieve appropriate file info if rank > 1
    my_files, my_nfiles = comm.recv(source=0)

# Broadcast meta data
if rank == 0:
    meta = [dirname, radatadir, nr, nr_cz, nr_rz, rw, rw_cz, rw_rz,\
            ir_bcz, rhot]
else:
    meta = None
dirname, radatadir, nr, nr_cz, nr_rz, rw, rw_cz, rw_rz, ir_bcz, rhot =\
    comm.bcast(meta, root=0)

# Checkpoint and time
comm.Barrier()
if rank == 0:
    t2 = time.time()
    print (format_time(t2 - t1))
    print ('Considering %i %s files for the trace: %s through %s'\
        %(nfiles, dataname, file_list[0], file_list[-1]))
    print(fill_str('computing', lent, char), end='\r')
    t1 = time.time()

# Now analyze the data
my_times = []
my_iters = []
my_vals = []
my_vals_cz = []
my_vals_rz = []

for i in range(my_nfiles):
    a = reading_func(radatadir + str(my_files[i]).zfill(8), '')
    for j in range(a.niter):
        vals_loc = np.copy(a.vals[:, 0, :, j])

        # Get the values in the CZ/RZ separately
        vals_cz_loc = vals_loc[:ir_bcz + 1]
        vals_rz_loc = vals_loc[ir_bcz + 1:]

        gav = np.sum(rw*vals_loc, axis=0)
        gav_cz = np.sum(rw_cz*vals_cz_loc, axis=0)
        gav_rz = np.sum(rw_rz*vals_rz_loc, axis=0)

        my_vals.append(gav)
        my_vals_cz.append(gav_cz)
        my_vals_rz.append(gav_rz)

        my_times.append(a.time[j])
        my_iters.append(a.iters[j])
    if rank == 0:
        pcnt_done = i/my_nfiles*100.0
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
    vals_cz = []
    vals_rz = []

    # Gather the results into these "master" arrays
    for j in range(nproc):
        if j >= 1:
            # Get data from rank j
            my_times, my_iters, my_vals, my_vals_cz, my_vals_rz = comm.recv(source=j)
        times += my_times
        iters += my_iters
        vals += my_vals
        vals_cz += my_vals_cz
        vals_rz += my_vals_rz
else: # other processes send their data
    comm.send([my_times, my_iters, my_vals, my_vals_cz, my_vals_rz], dest=0)

# Make sure proc 0 collects all data
comm.Barrier()

# proc 0 saves the data
if rank == 0:
    # create data directory if it doesn't already exist
    datadir = dirname + '/data/'
    if not os.path.isdir(datadir):
        os.makedirs(datadir)

    # Set the timetrace savename
    savename = 'G_Avgs_trace_2dom-' + file_list[0] + '_' +\
            file_list[-1] + '.pkl'
    savefile = datadir + savename

    # save the data
    f = open(savefile, 'wb')
    # convert everything to arrays
    vals = np.array(vals)
    vals_cz = np.array(vals_cz)
    vals_rz = np.array(vals_rz)
    times = np.array(times)
    iters = np.array(iters)
    pickle.dump({'vals': vals, 'vals_cz': vals_cz, 'vals_rz': vals_rz, 'times': times, 'iters': iters, 'lut': a.lut, 'qv': a.qv}, f, protocol=4)
    f.close()
    t2 = time.time()
    print (format_time(t2 - t1))
    print(make_bold(fill_str('total time', lent, char)), end='')
    print (make_bold(format_time(t2 - t1_glob)))
    print ('data saved at ')
    print (make_bold(savefile))
