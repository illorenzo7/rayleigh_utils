##################################################################
# Routine to average Rayleigh data in time (generic)
# Author: Loren Matilsky
# Created: 04/08/2021
##################################################################
# This routine computes the time average of magnetic torques from
# Rayleigh data (AZ_Avgs, Meridional Slices), breaking up terms 
# into separate contributions from radial/latitudinal fluxes
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
my_vals = np.zeros((nt, nr, 18)) # mm and pp: full torques, r and theta
# contributions to full torques, 3-term contributions to r and theta:
# 2 + 2 x 2 + 2 x 2 x 3 = 2 x (1 + 2 + 6) = 2 x 9 = 18

for i in range(my_nfiles):
    a = reading_func1(radatadir1 + str(my_files[i]).zfill(8), '')
    mer = reading_func2(radatadir2 + str(my_files[i]).zfill(8), '')

    my_weight = 1.0/(nfiles*a.niter)
    for j in range(a.niter):
        # mean fields
        br_m = a.vals[:, :, a.lut[801], 0].reshape((1, nt, nr))
        bt_m = a.vals[:, :, a.lut[802], 0].reshape((1, nt, nr))
        bp_m = a.vals[:, :, a.lut[803], 0].reshape((1, nt, nr))

        # mean derivs
        dbrdr_m = drad(br_m, rr)
        dbtdt_m = dth(bt_m, tt)
        dbpdr_m = drad(bp_m, rr)
        dbpdt_m = dth(bp_m, tt)

        # full fields
        br = mer.vals[:, :, :, mer.lut[801], 0]
        bt = mer.vals[:, :, :, mer.lut[802], 0]
        bp = mer.vals[:, :, :, mer.lut[803], 0]

        # full derivs
        dbrdr = drad(br, rr)
        dbtdt = dth(bt, tt)
        dbpdr = drad(bp, rr)
        dbpdt = dth(bp, tt)

        # fluc fields + derivs
        br_p = br - br_m
        bt_p = bt - bt_m
        bp_p = bp - bp_m

        dbrdr_p = dbrdr - dbrdr_m
        dbtdt_p = dbtdt - dbtdt_m
        dbpdr_p = dbpdr - dbpdr_m
        dbpdt_p = dbpdt - dbpdt_m

        # various contributions to torque
        c4 = 1./(4.*np.pi)

        tau_mm_r_1 = np.mean(c4*sint_3d*rr_3d*br_m*dbpdr_m, axis=0)*my_weight
        tau_mm_r_2 = np.mean(c4*sint_3d*rr_3d*bp_m*dbrdr_m, axis=0)*my_weight
        tau_mm_r_3 = np.mean(3.*c4*sint_3d*br_m*bp_m, axis=0)*my_weight
        tau_mm_r = tau_mm_r_1 + tau_mm_r_2 + tau_mm_r_3

        tau_pp_r_1 = np.mean(c4*sint_3d*rr_3d*br_p*dbpdr_p, axis=0)*my_weight
        tau_pp_r_2 = np.mean(c4*sint_3d*rr_3d*bp_p*dbrdr_p, axis=0)*my_weight
        tau_pp_r_3 = np.mean(3.*c4*sint_3d*br_p*bp_p, axis=0)*my_weight
        tau_pp_r = tau_pp_r_1 + tau_pp_r_2 + tau_pp_r_3

        tau_mm_t_1 = np.mean(c4*sint_3d*bt_m*dbpdt_m, axis=0)*my_weight
        tau_mm_t_2 = np.mean(c4*sint_3d*bp_m*dbtdt_m, axis=0)*my_weight
        tau_mm_t_3 = np.mean(2.*c4*cost_3d*bt_m*bp_m, axis=0)*my_weight
        tau_mm_t = tau_mm_t_1 + tau_mm_t_2 + tau_mm_t_3

        tau_pp_t_1 = np.mean(c4*sint_3d*bt_p*dbpdt_p, axis=0)*my_weight
        tau_pp_t_2 = np.mean(c4*sint_3d*bp_p*dbtdt_p, axis=0)*my_weight
        tau_pp_t_3 = np.mean(2.*c4*cost_3d*bt_p*bp_p, axis=0)*my_weight
        tau_pp_t = tau_pp_t_1 + tau_pp_t_2 + tau_pp_t_3

        tau_mm = tau_mm_r + tau_mm_t
        tau_pp = tau_pp_r + tau_pp_t

        # mean + fluc decomp
        my_vals[:, :, 0] += tau_mm
        my_vals[:, :, 1] += tau_pp

        # r + theta decomp.
        my_vals[:, :, 2] += tau_mm_r
        my_vals[:, :, 3] += tau_mm_t

        my_vals[:, :, 4] += tau_pp_r
        my_vals[:, :, 5] += tau_pp_t

        # 3-term decomp
        my_vals[:, :, 6] += tau_mm_r_1
        my_vals[:, :, 7] += tau_mm_r_2
        my_vals[:, :, 8] += tau_mm_r_3

        my_vals[:, :, 9] += tau_mm_t_1
        my_vals[:, :, 10] += tau_mm_t_2
        my_vals[:, :, 11] += tau_mm_t_3

        my_vals[:, :, 12] += tau_pp_r_1
        my_vals[:, :, 13] += tau_pp_r_2
        my_vals[:, :, 14] += tau_pp_r_3

        my_vals[:, :, 15] += tau_pp_t_1
        my_vals[:, :, 16] += tau_pp_t_2
        my_vals[:, :, 17] += tau_pp_t_3

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
    savename = 'mag_torque-' +\
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
