# Routine to trace Rayleigh G_Avgs data in time
# Created by: Loren Matilsky
# On: 12/18/2018
# Parallelized: 11/30/2020
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
sys.path.append(os.environ['raco'])
from rayleigh_diagnostics import AZ_Avgs, Meridional_Slices
from common import drad, dth
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
    print ('%8.2e s' %(t2 - t1))
    print(fill_str('proc 0 distributing the file lists', lent, char),\
            end='')
    t1 = time.time()

# proc 0 reads the file lists and distributes them
if rank == 0:
    # Get the name of the run directory
    dirname = sys.argv[1]

    # Get the Rayleigh data directory
    radatadir1 = dirname + '/' + dataname1 + '/'
    radatadir2 = dirname + '/' + dataname2 + '/'

    # Get all the file names in datadir and their integer counterparts
    file_list, int_file_list, nfiles = get_file_lists(radatadir1)

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
    a0 = reading_func1(radatadir1 + file_list[0], '')
    nrec_full = a0.niter

    # Will need nrec (last niter) to get proper time axis size
    af = reading_func1(radatadir1 + file_list[-1], '')
    nrec_last = af.niter
    ntimes = (nfiles - 1)*nrec_full + nrec_last

    # get grid information
    rr = a0.radius
    ri, ro = np.min(rr), np.max(rr)
    d = ro - ri
    rr_depth = (ro - rr)/d
    rr_height = (rr - ri)/d
    sint = a0.sintheta
    cost = a0.costheta
    tt = np.arccos(cost)
    tt_lat = (np.pi/2 - tt)*180/np.pi
    nr = a0.nr
    nt = a0.ntheta

    # compute some derivative quantities for the grid
    tt_2d, rr_2d = np.meshgrid(tt, rr, indexing='ij')
    sint_2d = np.sin(tt_2d); cost_2d = np.cos(tt_2d)
    xx = rr_2d*sint_2d
    zz = rr_2d*cost_2d

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
    meta = [dirname, radatadir1, radatadir2, nt, nr, ntimes, rr, rr_2d,\
            tt, tt_2d, sint_2d]
else:
    meta = None
dirname, radatadir1, radatadir2, nt, nr, ntimes, rr, rr_2d,\
    tt, tt_2d, sint_2d = comm.bcast(meta, root=0)

# Checkpoint and time
comm.Barrier()
if rank == 0:
    t2 = time.time()
    print ('%8.2e s' %(t2 - t1))
    print ('Considering %i %s/%s files for the trace: %s through %s'\
        %(nfiles, dataname1, dataname2, file_list[0], file_list[-1]))
    print ("no. slices = %i" %ntimes)
    print(fill_str('computing', lent, char), end='\r')
    t1 = time.time()

# Now analyze the data
nq = 16
my_vals = np.zeros((nt, nr, nq))
# "my_vals will be a weighted sum"
my_weight = 1./ntimes

for i in range(my_nfiles):
    a = reading_func1(radatadir1 + str(my_files[i]).zfill(8), '')
    mer = reading_func2(radatadir2 + str(my_files[i]).zfill(8), '')
    # take mean along the time axis;

    for j in range(a.niter):
        # mean B
        br_m = a.vals[:, :, a.lut[801], j].reshape((1, nt, nr))
        bt_m = a.vals[:, :, a.lut[802], j].reshape((1, nt, nr))
        bp_m = a.vals[:, :, a.lut[803], j].reshape((1, nt, nr))

        # mean v
        vr_m = a.vals[:, :, a.lut[1], j].reshape((1, nt, nr))
        vt_m = a.vals[:, :, a.lut[2], j].reshape((1, nt, nr))
        vp_m = a.vals[:, :, a.lut[3], j].reshape((1, nt, nr))

        # full B
        br = mer.vals[:, :, :, mer.lut[801], j]
        bt = mer.vals[:, :, :, mer.lut[802], j]
        bp = mer.vals[:, :, :, mer.lut[803], j]

        # full v
        vr = mer.vals[:, :, :, mer.lut[1], j]
        vt = mer.vals[:, :, :, mer.lut[2], j]
        vp = mer.vals[:, :, :, mer.lut[3], j]

        # fluc B
        br_p = br - br_m
        bt_p = bt - bt_m
        bp_p = bp - bp_m

        # fluc v
        vr_p = vr - vr_m
        vt_p = vt - vt_m
        vp_p = vp - vp_m

        # correlations
        vrbt_mm = np.mean(vr_m*bt_m, axis=0)
        vrbt_pp = np.mean(vr_p*bt_p, axis=0)
        vrbp_mm = np.mean(vr_m*bp_m, axis=0)
        vrbp_pp = np.mean(vr_p*bp_p, axis=0)
    
        vtbr_mm = np.mean(vt_m*br_m, axis=0)
        vtbr_pp = np.mean(vt_p*br_p, axis=0)
        vtbp_mm = np.mean(vt_m*bp_m, axis=0)
        vtbp_pp = np.mean(vt_p*bp_p, axis=0)

        vpbr_mm = np.mean(vp_m*br_m, axis=0)
        vpbr_pp = np.mean(vp_p*br_p, axis=0)
        vpbt_mm = np.mean(vp_m*bt_m, axis=0)
        vpbt_pp = np.mean(vp_p*bt_p, axis=0)
        
        # induction terms
        # radial
        my_vals[:, :, 0] += 1./rr_2d/sint_2d*\
            dth(sint_2d*vrbt_mm, tt)*my_weight
        my_vals[:, :, 1] += -1./rr_2d/sint_2d*\
            dth(sint_2d*vtbr_mm, tt)*my_weight
        my_vals[:, :, 2] += 1./rr_2d/sint_2d*\
            dth(sint_2d*vrbt_pp, tt)*my_weight
        my_vals[:, :, 3] += -1./rr_2d/sint_2d*\
            dth(sint_2d*vtbr_pp, tt)*my_weight
        # theta
        ind_off = 4
        my_vals[:, :, ind_off + 0] +=\
            1./rr_2d*drad(rr_2d*vtbr_mm, rr)*my_weight
        my_vals[:, :, ind_off + 1] +=\
            -1./rr_2d*drad(rr_2d*vrbt_mm, rr)*my_weight
        my_vals[:, :, ind_off + 2] +=\
            1./rr_2d*drad(rr_2d*vtbr_pp, rr)*my_weight
        my_vals[:, :, ind_off + 3] +=\
            -1./rr_2d*drad(rr_2d*vrbt_pp, rr)*my_weight
        # phi
        ind_off += 4
        my_vals[:, :, ind_off + 0] +=\
            1./rr_2d*drad(rr_2d*vpbr_mm, rr)*my_weight
        my_vals[:, :, ind_off + 1] +=\
            -1./rr_2d*drad(rr_2d*vrbp_mm, rr)*my_weight
        my_vals[:, :, ind_off + 2] += 1./rr_2d*dth(vpbt_mm, tt)*my_weight
        my_vals[:, :, ind_off + 3] += -1./rr_2d*dth(vtbp_mm, tt)*my_weight
        my_vals[:, :, ind_off + 4] +=\
            1./rr_2d*drad(rr_2d*vpbr_pp, rr)*my_weight
        my_vals[:, :, ind_off + 5] +=\
            -1./rr_2d*drad(rr_2d*vrbp_pp, rr)*my_weight
        my_vals[:, :, ind_off + 6] += 1./rr_2d*dth(vpbt_pp, tt)*my_weight
        my_vals[:, :, ind_off + 7] += -1./rr_2d*dth(vtbp_pp, tt)*my_weight

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
    # Initialize zero-filled 'vals' array to store the data
    vals = np.zeros((nt, nr, nq))

    # Gather the results into this "master" array
    for j in range(nproc):
        if j >= 1:
            # Get my_ntimes, my_times, my_iters, my_vals from rank j
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
    datadir = dirname + '/data/'
    if not os.path.isdir(datadir):
        os.makedirs(datadir)

    # Set the timetrace savename by the directory, what we are saving,
    # and first and last iteration files for the trace
    dirname_stripped = strip_dirname(dirname)
    savename = dirname_stripped + '_induction_terms_' +\
            file_list[0] + '_' + file_list[-1] + '.pkl'
    savefile = datadir + savename

    # save the data
    # Get first and last iters of files
    iter1, iter2 = int_file_list[0], int_file_list[-1]
    f = open(savefile, 'wb')
    pickle.dump({'vals': vals, 'lut': a0.lut, 'count': ntimes, 'iter1': iter1, 'iter2': iter2, 'qv': a0.qv, 'nq': a0.nq,  'rr': rr, 'rr_depth': rr_depth, 'rr_height': rr_height, 'nr': nr, 'ri': ri, 'ro': ro, 'd': d, 'tt': tt, 'tt_lat': tt_lat, 'sint': sint, 'cost': cost,'nt': nt, 'rr_2d': rr_2d, 'tt_2d': tt_2d, 'sint_2d': sint_2d, 'cost_2d': cost_2d, 'xx': xx, 'zz': zz}, f, protocol=4)
    f.close()
    t2 = time.time()
    print ('%8.2e s' %(t2 - t1))
    print(fill_str('total time', lent, char), end='')
    print ('%8.2e s' %(t2 - t1_glob))
    print ('data saved at ')
    print (savefile)
