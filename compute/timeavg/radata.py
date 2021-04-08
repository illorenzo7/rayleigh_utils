# Routine to average Rayleigh data in time (generic)
# Author: Loren Matilsky
# 04/08/2021
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

# import the rading routines
from rayleigh_diagnostics import AZ_Avgs, Shell_Avgs, G_Avgs,\
        Shell_Spectra, Shell_Slices, Meridional_Slices, Equatorial_Slices

# Broadcast the desired datatype
if rank == 0:
    args = sys.argv[2:]
    clas = read_clas(args)
    dtype = clas['dtype'].val
else:
    dtype = None
dtype = comm.bcast(dtype, root=0)

if dtype == 'azav':
    reading_func = AZ_Avgs
    dataname = 'AZ_Avgs'
if dtype == 'shav':
    reading_func = Shell_Avgs
    dataname = 'Shell_Avgs'
if dtype == 'gav':
    reading_func = G_Avgs
    dataname = 'G_Avgs'
if dtype == 'specav':
    reading_func = Shell_Spectra
    dataname = 'Shell_Spectra'
if dtype == 'ssav':
    reading_func = Shell_Slices
    dataname = 'Shell_Slices'
if dtype == 'merav':
    reading_func = Meridional_Slices
    dataname = 'Meridional_Slices'
if dtype == 'eqav':
    reading_func = Equatorial_Slices
    dataname = 'Equatorial_Slices'

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

    # Get the Rayleigh data directory
    radatadir = dirname + '/' + dataname + '/'

    # Get all the file names in datadir and their integer counterparts
    file_list, int_file_list, nfiles = get_file_lists(radatadir)

    # get CLAs
    args = sys.argv[2:]

    # get desired analysis range
    the_tuple = get_desired_range(int_file_list, args)
    if the_tuple is None:
        index_first, index_last = nfiles - 101, nfiles - 1  
        # By default trace over the last 100 files
    else:
        index_first, index_last = the_tuple

    datadir = clas['datadir'].val
    if datadir is None:
        datadir = dirname + '/data/'

    # Remove parts of file lists we don't need
    file_list = file_list[index_first:index_last + 1]
    int_file_list = int_file_list[index_first:index_last + 1]
    nfiles = index_last - index_first + 1
    weight = 1.0/nfiles

    # Get the problem size
    nproc_min, nproc_max, n_per_proc_min, n_per_proc_max =\
            opt_workload(nfiles, nproc)

    # Will need the first data file for array shape
    # (exclude the last axis---time axis---from the shape
    a0 = reading_func(radatadir + file_list[0], '')
    shape = np.shape(a0.vals[..., 0])

    # Distribute file_list to each process
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

        # send  my_files nproc > 1
        if k >= 1:
            comm.send(my_files, dest=k)
else: # recieve my_files
    my_files = comm.recv(source=0)

# Broadcast dirname, radatadir, nq, etc.
if rank == 0:
    meta = [dirname, radatadir, shape, weight]
else:
    meta = None
dirname, radatadir, shape, weight = comm.bcast(meta, root=0)

# Checkpoint and time
comm.Barrier()
if rank == 0:
    t2 = time.time()
    print (format_time(t2 - t1))
    print ('Considering %i %s files for the average: %s through %s'\
        %(nfiles, dataname, file_list[0], file_list[-1]))
    print(fill_str('computing', lent, char), end='\r')
    t1 = time.time()

# Now analyze the data
my_vals = np.zeros(shape)
# "my_vals will be a weighted sum"
my_nfiles = len(my_files)

for i in range(my_nfiles):
    if rank == 0 and i == 0:
        a = a0
    else:   
        a = reading_func(radatadir + str(my_files[i]).zfill(8), '')
    # take mean along the time axis, which is always the last axis
    my_vals += np.mean(a.vals, axis=len(shape))*weight
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
    # Initialize zero-filled 'vals' array to store the data
    vals = np.zeros(shape)

    # Gather the results into this "master" array
    for j in range(nproc):
        if j >= 1:
            # Get my_vals from rank j
            my_vals = comm.recv(source=j)
        # "my_vals" are all weighted: their sum equals the overall average
        vals += my_vals 

else: # other processes send their data
    comm.send(my_vals, dest=0)

# Make sure proc 0 collects all data
comm.Barrier()

# proc 0 saves the data
if rank == 0:
    # create data directory if it doesn't already exist
    if datadir is None:
        datadir = dirname + '/data/'
        if not os.path.isdir(datadir):
            os.makedirs(datadir)

    # set the save name by what we are saving (dataname)
    # and first and last iteration files for the trace
    iter1, iter2 = int_file_list[0], int_file_list[-1]
    savename = dataname + '_' + str(iter1).zfill(8) + '_' +\
            str(iter2).zfill(8) + '.pkl'
    savefile = datadir + savename

    # save the data
    f = open(savefile, 'wb')
    pickle.dump({'vals': vals, 'lut': a0.lut, 'count': nfiles, 'iter1': iter1, 'iter2': iter2, 'qv': a0.qv, 'shape': shape}, f, protocol=4)
    f.close()
    t2 = time.time()
    print (format_time(t2 - t1))
    print(make_bold(fill_str('total time', lent, char)), end='')
    print (make_bold(format_time(t2 - t1_glob)))
    print ('data saved at ')
    print (make_bold(savefile))
