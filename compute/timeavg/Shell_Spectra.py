# Routine to average Rayleigh Shell_Spectra data in time
# For the full power (both l, m) the "power" starts complex, so the 
# square is taken before averaging, and the square root may be taken 
# at the end to get rms power in (l, m) space
#
# For the lpower, the square has already been taken in the Rayleigh
# diagnostics.
# Created by: Loren Matilsky
# On: 12/07/2018
# Parallelized: 12/04/2020
##################################################################
# This routine computes the time average of the Shell_Spectra data 
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
    from common import fill_str
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
from rayleigh_diagnostics import Shell_Spectra, GridInfo
reading_func = Shell_Spectra
dataname = 'Shell_Spectra'

if rank == 0:
    # modules needed only by proc 0 
    import pickle
from common import *

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
        # By default average over the last 100 files
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
    nq = a0.nq  # will add the three internal energies after

    # Will need nrec (last niter) to get proper time axis size
    af = reading_func(radatadir + file_list[-1], '')
    nrec_last = af.niter
    ntimes = (nfiles - 1)*nrec_full + nrec_last

    # get grid information
    # note: full radius info not contained in Shell_Spectra objects
    # Get from the file grid_info
    gi = GridInfo(dirname + '/grid_info')
    rr = gi.radius
    ri, ro = np.min(rr), np.max(rr)
    shell_depth = ro - ri

    # spectral l/m values: useful arrays/paremeters
    lmax = a0.lmax
    mmax = a0.mmax

    nell = a0.nell
    nm = a0.nm
    nr = a0.nr
    nq = a0.nq

    lvals = np.arange(lmax + 1, dtype='float')
    mvals = np.arange(mmax + 1, dtype='float')

    lvals_2d, mvals_2d = np.meshgrid(lvals, mvals, indexing='ij')

    # Compute some quantities having to do with the radii sampled
    rvals = a0.radius
    rvals_depth = (ro - rvals)/shell_depth
    rvals_height = (rvals - ri)/shell_depth

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
            comm.send([my_files, my_nfiles, my_ntimes, ntimes], dest=k)
else: # recieve my_files, my_nfiles, my_ntimes
    my_files, my_nfiles, my_ntimes, ntimes = comm.recv(source=0)

# Broadcast dirname, radatadir, nq, etc.
if rank == 0:
    meta = [dirname, radatadir, nq, nell, nm, nr, ntimes]
else:
    meta = None
dirname, radatadir, nq, nell, nm, nr, ntimes = comm.bcast(meta, root=0)

# Checkpoint and time
comm.Barrier()
if rank == 0:
    t2 = time.time()
    print ('%8.2e s' %(t2 - t1))
    print ('Considering %i %s files for the average: %s through %s'\
        %(nfiles, dataname, file_list[0], file_list[-1]))
    print ("no. slices = %i" %ntimes)
    print(fill_str('computing', lent, char), end='\r')
    t1 = time.time()

# Now analyze the data
# empty arrays to hold partial results
my_fullpower = np.zeros((nell, nm, nr, nq)) 
my_lpower = np.zeros((nell, nr, nq, 3))
# values will be a portion of a weighted average
my_weight = my_ntimes/ntimes

for i in range(my_nfiles):
    if rank == 0 and i == 0:
        a = a0
    else:   
        a = reading_func(radatadir + str(my_files[i]).zfill(8), '')
    # take means along the time axis, and weight them
    my_fullpower = np.mean(np.abs(a.vals)**2.0, axis=4)*my_weight
    my_lpower = np.mean(a.lpower, axis=3)*my_weight
    if rank == 0:
        pcnt_done = (i + 1)/my_nfiles*100.
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
    # empty arrays to hold full results
    fullpower = np.zeros((nell, nm, nr, nq)) 
    lpower = np.zeros((nell, nr, nq, 3))

    # Gather the results into this "master" array
    for j in range(nproc):
        if j >= 1:
            # Get my_ntimes, my_times, my_iters, my_vals from rank j
            my_fullpower = comm.recv(source=j)
            my_lpower = comm.recv(source=j)
        # "my" values are all weighted: their sum equals the overall average
        fullpower += my_fullpower
        lpower += my_lpower
else: # other processes send their data
    comm.send(my_fullpower, dest=0)
    comm.send(my_lpower, dest=0)

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
    savename = dirname_stripped + '_Shell_Spectra_' +\
            file_list[0] + '_' + file_list[-1] + '.pkl'
    savefile = datadir + savename

    # save the data
    # Get first and last iters of files
    iter1, iter2 = int_file_list[0], int_file_list[-1]
    f = open(savefile, 'wb')
    pickle.dump({'fullpower': fullpower, 'lpower': lpower, 'lut': a0.lut, 'count': ntimes, 'iter1': iter1, 'iter2': iter2, 'qv': a0.qv, 'nq': nq, 'rinds': a0.inds, 'rvals': rvals, 'rvals_depth': rvals_depth, 'rvals_height': rvals_height, 'rr': rr, 'nr': nr, 'ri': ri, 'ro': ro, 'shell_depth': shell_depth, 'lvals': lvals, 'lvals_2d': lvals_2d, 'nell': nell, 'lmax': lmax, 'mvals': mvals, 'mvals_2d': mvals_2d, 'nm': nm, 'mmax': mmax}, f, protocol=4)
    f.close()
    t2 = time.time()
    print ('%8.2e s' %(t2 - t1))
    print(fill_str('total time', lent, char), end='')
    print ('%8.2e s' %(t2 - t1_glob))
