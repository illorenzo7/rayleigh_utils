# Routine to trace AZ_Avgs in time/latitude space (pick different radii)
# Created by: Loren Matilsky
# On: 02/27/2019
# Parallelized: 12/12/2020
##################################################################
# This routine computes the trace in time/latitude of quantities in the 
# AZ_Avgs data for a particular simulation. 
# 
# By default, the quantities traced are the nq (=5, or 8) fluid variables 
# (vr, vt, vp, s, p, [bp, bt, and bp] (if magnetism present).
# This may be changed via the '-qvals' CLA, e.g., -qvals '1 2 3 1425 1427'.
#
# By default, the 8 variables are computed at each time (for all latitudes)
# at ndepths = 9 depths equally spaced throughout the shell (like the shell
# slice levels).
# This may be changed via the '-depths' CLA, e.g., -depths '0 0.4 0.75 0.95'
# or -rzquarter, -rzhalf, -rz75, -rz1 (depth RZ = 1/4 depth CZ, 1/2, 3/4, 1)
# will generate 9 depths in each zone
# (assuming rm = 5.0 x 10^10 cm, ro = 6.5860209 x 10^10 cm)
# also can be changed with -rvals '.719 .863' [radii units of rsun]
# or -rvalscm '5e10 5.8e10' [radii units of cm]
#
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

# proc 0 reads the file lists and distributes them, also the meta data
if rank == 0:
    # get the name of the run directory + CLAs
    dirname = sys.argv[1]
    args = sys.argv[2:]
    clas = read_clas(dirname, args)

    # get the Rayleigh data directory
    radatadir = dirname + '/' + dataname + '/'

    # get desired analysis range
    file_list, int_file_list, nfiles = get_file_lists(radatadir)
    index_first, index_last = get_desired_range(int_file_list, args)
    file_list = file_list[index_first:index_last + 1]
    int_file_list = int_file_list[index_first:index_last + 1]
    nfiles = index_last - index_first + 1

    # get the problem size
    nproc_min, nproc_max, n_per_proc_min, n_per_proc_max =\
            opt_workload(nfiles, nproc)

    # Get metadata
    a0 = reading_func(radatadir + file_list[0], '')
    nrec_full = a0.niter
    qvals = clas['qvals']
    nq = len(qvals)

    # Get bunch of grid info
    rr = a0.radius
    ri, ro = np.min(rr), np.max(rr)
    d = ro - ri
    rr_depth = (ro - rr)/d
    rr_height = (rr - ri)/d
    sint = a0.sintheta
    cost = a0.costheta
    tt = np.arccos(cost)
    tt_lat = (np.pi/2 - tt)*180./np.pi
    nr = a0.nr
    nt = a0.ntheta

    # compute some derivative quantities for the grid
    tt_2d, rr_2d = np.meshgrid(tt, rr, indexing='ij')
    sint_2d = np.sin(tt_2d); cost_2d = np.cos(tt_2d)
    xx = rr_2d*sint_2d
    zz = rr_2d*cost_2d

    # get r-indices associated with desired rvals
    irvals = []
    for rval in clas['rvals']:
        irvals.append(np.argmin(np.abs(rr/rsun - rval)))
    irvals = np.array(irvals)
    # recompute the actual rvals we get
    rvals = rr[irvals]/rsun
    nrvals = len(rvals)
    print ("sampling radii at")
    print ("rvals = " + arr_to_str(rvals, '%1.3f'))

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

# broadcast meta data
if rank == 0:
    meta = [dirname, radatadir, qvals, nq, irvals, nrvals, nt]
else:
    meta = None
dirname, radatadir, qvals, nq, irvals, nrvals, nt = comm.bcast(meta, root=0)

# Checkpoint and time
comm.Barrier()
if rank == 0:
    t2 = time.time()
    print (format_time(t2 - t1))
    print ('Considering %i %s files for the trace: %s through %s'\
        %(nfiles, dataname, file_list[0], file_list[-1]))
    print ("ntimes for trace = %i" %ntimes)
    print(fill_str('computing', lent, char), end='\r')
    t1 = time.time()

# Now analyze the data
my_times = np.zeros(my_ntimes)
my_iters = np.zeros(my_ntimes, dtype='int')
my_vals = np.zeros((my_ntimes, nt, nrvals, nq))

my_count = 0
for i in range(my_nfiles):
    if rank == 0 and i == 0:
        a = a0
    else:   
        a = reading_func(radatadir + str(my_files[i]).zfill(8), '')
    for j in range(a.niter):
        if my_count < my_ntimes: # make sure we don't go over the allotted
            # space in the arrays
            my_vals[my_count, :, :, :] =\
                    a.vals[:, :, :, j][:, irvals, :][:, :, a.lut[qvals]]
            my_times[my_count] = a.time[j] 
            my_iters[my_count] = a.iters[j]
        my_count += 1
    if rank == 0:
        pcnt_done = my_count/my_ntimes*100.0
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
    vals = np.zeros((ntimes, nt, nrvals, nq))
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
    # and first and last iteration files for the trace
    savename = 'timelat' + clas['tag'] + '-' +\
            file_list[0] + '_' + file_list[-1] + '.pkl'
    savefile = datadir + savename

    # save the data
    f = open(savefile, 'wb')
    di_sav = {'vals': vals, 'times': times, 'iters': iters,\
    'rvals': rvals, 'qvals': qvals}
    pickle.dump(di_sav, f, protocol=4)
    f.close()
    t2 = time.time()
    print (format_time(t2 - t1))
    print(make_bold(fill_str('total time', lent, char)), end='')
    print (make_bold(format_time(t2 - t1_glob)))
    print ('data saved at ')
    print (make_bold(savefile))
