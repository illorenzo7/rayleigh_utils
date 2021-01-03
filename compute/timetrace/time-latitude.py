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
# This may be changed via the '-vars' CLA, e.g., -vars '1 2 3 1425 1427'.
#
# By default, the 8 variables are computed at each time (for all latitudes)
# at ndepths = 9 depths equally spaced throughout the shell (like the shell
# slice levels).
# This may be changed via the '-depths' CLA, e.g., -depths '0 0.4 0.75 0.95'
# or -rzquarter, -rzhalf, -rz75, -rz1 (depth RZ = 1/4 depth CZ, 1/2, 3/4, 1)
# will generate 9 depths in each zone
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
    from get_parameter import get_parameter
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
from rayleigh_diagnostics import AZ_Avgs
reading_func = AZ_Avgs
dataname = 'AZ_Avgs'

if rank == 0:
    # modules needed only by proc 0 
    import pickle
    from common import get_file_lists, get_desired_range, strip_dirname
    from mpi_util import opt_workload

# Checkpoint and time
comm.Barrier()
if rank == 0:
    t2 = time.time()
    print ('%8.2e s' %(t2 - t1))
    print(fill_str('proc 0 distributing the file lists', lent, char),\
            end='')
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

    the_tuple = get_desired_range(int_file_list, args)
    if the_tuple is None:
        index_first, index_last = nfiles - 101, nfiles - 1  
        # By default trace over the last 100 files
    else:
        index_first, index_last = the_tuple

    # Set other defaults
    qvals = [1, 2, 3, 301, 302]
    magnetism = get_parameter(dirname, 'magnetism')
    if magnetism:
        qvals.append(801)
        qvals.append(802)
        qvals.append(803)
    depths = [0.05, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 0.95]
    tag = ''
    for i in range(nargs):
        arg = args[i]
        if arg == '-vars':
            qvals = []
            qvals_str = args[i+1].split() 
            for qval_str in qvals_str:
                qvals.append(int(qval_str))
        elif arg == '-depths':
            depths = []
            depths_str = args[i+1].split()
            for depth_str in depths_str:
                depths.append(float(depth_str))
        elif arg == '-tag':
            tag = args[i+1] + '_'
        elif arg == '-rzquarter': # 9 depths in RZ and CZ, with RZ depth
            # 0.25 of CZ depth
            print("Taking 9 depths in CZ and RZ each")
            print("assuming depth RZ = (1/4) depth CZ")
            depths = [0.04, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.76,\
                   0.81, 0.825, 0.85, 0.875, 0.9, 0.925, 0.95, 0.975, 0.990]
        elif arg == '-rzhalf': # 9 depths in RZ and CZ, with RZ depth
            # 0.5 of CZ depth
            print("Taking 9 depths in CZ and RZ each")
            print("assuming depth RZ = (1/2) depth CZ")
            depths =[0.03333333, 0.08333333, 0.16666667, 0.25, 0.33333333,\
                0.41666667, 0.5, 0.58333333, 0.63333333, 0.68333333,\
                0.70833333, 0.75, 0.79166667, 0.83333333, 0.875,\
                0.91666667, 0.95833333, 0.9833333]
        elif arg == '-rz75': # 9 depths in RZ and CZ, with RZ depth
            # 0.75 of CZ depth
            print("Taking 9 depths in CZ and RZ each")
            print("assuming depth RZ = (3/4) depth CZ")
            depths = 1.0 - np.array([0.02142857, 0.05357143, 0.10714286,\
                    0.16071429, 0.21428571, 0.26785714, 0.32142857, 0.375,\
                    0.40714286, 0.45714286, 0.5, 0.57142857, 0.64285714,\
                    0.71428571, 0.78571429, 0.85714286, 0.92857143,\
                    0.97142857])
            depths = depths.tolist()
        elif arg == '-rz1': # 9 depths in RZ and CZ, with RZ depth
            # equal to CZ depth
            print("Taking 9 depths in CZ and RZ each")
            print("assuming depth RZ = depth CZ")
            depths = 1.0 - np.array([0.025, 0.0625, 0.125, 0.1875, 0.25,\
                    0.3125, 0.375, 0.4375, 0.475, 0.525, 0.5625, 0.625,\
                    0.6875, 0.75, 0.8125, 0.875, 0.9375, 0.975])
            depths = depths.tolist()
        elif arg == '-torques':
            print("tracing over TORQUES")
            qvals = [3, 1801, 1802, 1803, 1804, 1819]
            if magnetism:
                qvals.append(1805)
                qvals.append(1806)
            tag = 'torques' + '_'
        elif arg == '-induction':
            print("tracing over INDUCTION QUANTITIES")
            qvals = [1604, 1605, 1609, 1610, 1614, 1615, 1619, 1620, 1623,\
                    1624, 1629, 1630, 1601, 1602, 1603, 1606, 1607, 1608,\
                    1611, 1612, 1613, 1616, 1617, 1618, 1621, 1622, 1623,\
                    1626, 1627, 1628]
            tag = 'induction' + '_'
        elif arg == '-tag':
            tag = args[i+1] + '_'

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

    # Get metadata
    a0 = reading_func(radatadir + file_list[0], '')
    nrec_full = a0.niter
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
    ndepths = len(depths)
    nt = a0.ntheta

    # compute some derivative quantities for the grid
    tt_2d, rr_2d = np.meshgrid(tt, rr, indexing='ij')
    sint_2d = np.sin(tt_2d); cost_2d = np.cos(tt_2d)
    xx = rr_2d*sint_2d
    zz = rr_2d*cost_2d

    # get r-indices associated with depths
    rinds = []
    for depth in depths:
        rinds.append(np.argmin(np.abs(rr_depth - depth)))

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

# Broadcast dirname, radatadir, qvals, nq, rinds, ndepths, nt
if rank == 0:
    meta = [dirname, radatadir, qvals, nq, rinds, ndepths, nt]
else:
    meta = None
dirname, radatadir, qvals, nq, rinds, ndepths, nt = comm.bcast(meta, root=0)

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
my_vals = np.zeros((my_ntimes, nt, ndepths, nq))

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
                    a.vals[:, :, :, j][:, rinds, :][:, :, a.lut[qvals]]
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
    vals = np.zeros((ntimes, nt, ndepths, nq))
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
    dirname_stripped = strip_dirname(dirname)
    savename = dirname_stripped + '_time-latitude_' + tag +\
            file_list[0] + '_' + file_list[-1] + '.pkl'
    savefile = datadir + savename

    # save the data
    f = open(savefile, 'wb')
    # will also need first and last iteration files
    iter1, iter2 = int_file_list[0], int_file_list[-1]
    pickle.dump({'vals': vals, 'times': times, 'iters': iters,\
    'depths': depths,'qvals': qvals, 'rinds': rinds, 'niter': len(iters),\
    'ndepths': ndepths, 'nq': nq, 'iter1': iter1, 'iter2': iter2, 'rr': rr,\
    'rr_depth': rr_depth, 'rr_height': rr_height, 'nr': nr, 'ri': ri,\
    'ro': ro, 'd': d, 'tt': tt, 'tt_lat': tt_lat,\
    'sint': sint, 'cost': cost, 'ntheta': nt,\
    'rr_2d': rr_2d, 'tt_2d': tt_2d, 'sint_2d': sint_2d,\
    'cost_2d': cost_2d, 'xx': xx, 'zz': zz}, f, protocol=4)
    f.close()
    t2 = time.time()
    print ('%8.2e s' %(t2 - t1))
    print(fill_str('total time', lent, char), end='')
    print ('%8.2e s' %(t2 - t1_glob))
