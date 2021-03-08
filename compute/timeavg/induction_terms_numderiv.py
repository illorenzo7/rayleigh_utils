# more breakdown of poloidal terms (separate derivatives) than 
# induction_terms
# Created by: Loren Matilsky
# On: 01/16/2021
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
    cott = cost_2d/sint_2d
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
    meta = [dirname, radatadir1, radatadir2, nt, nr, ntimes, rr_2d,  cott]
else:
    meta = None
dirname, radatadir1, radatadir2, nt, nr, ntimes, rr_2d, cott =\
        comm.bcast(meta, root=0)

# Checkpoint and time
comm.Barrier()
if rank == 0:
    t2 = time.time()
    print (format_time(t2 - t1))
    print ('Considering %i %s/%s files for the trace: %s through %s'\
        %(nfiles, dataname1, dataname2, file_list[0], file_list[-1]))
    print ("no. slices = %i" %ntimes)
    print(fill_str('computing', lent, char), end='\r')
    t1 = time.time()

# Now analyze the data
nq = 45
my_vals = np.zeros((nt, nr, nq))
# "my_vals will be a weighted sum"
my_weight = 1./ntimes

for i in range(my_nfiles):
    a = reading_func1(radatadir1 + str(my_files[i]).zfill(8), '')
    mer = reading_func2(radatadir2 + str(my_files[i]).zfill(8), '')
    # take mean along the time axis;

    for j in range(a.niter):
        # full v
        vr = mer.vals[:, :, :, mer.lut[1], j]
        vt = mer.vals[:, :, :, mer.lut[2], j]
        vp = mer.vals[:, :, :, mer.lut[3], j]

        # full B
        br = mer.vals[:, :, :, mer.lut[801], j]
        bt = mer.vals[:, :, :, mer.lut[802], j]
        bp = mer.vals[:, :, :, mer.lut[803], j]

        # mean v
        vr_m = a.vals[:, :, a.lut[1], j].reshape((1, nt, nr))
        vt_m = a.vals[:, :, a.lut[2], j].reshape((1, nt, nr))
        vp_m = a.vals[:, :, a.lut[3], j].reshape((1, nt, nr))

        # mean B
        br_m = a.vals[:, :, a.lut[801], j].reshape((1, nt, nr))
        bt_m = a.vals[:, :, a.lut[802], j].reshape((1, nt, nr))
        bp_m = a.vals[:, :, a.lut[803], j].reshape((1, nt, nr))

        # full derivatives v
        dvrdr = drad_3d(vr, rr)
        dvtdr = drad_3d(vt, rr)
        dvpdr = drad_3d(vp, rr)

        # full derivatives v
        dvrdr = mer.vals[:, :, :, mer.lut[10], j]
        dvtdr = mer.vals[:, :, :, mer.lut[11], j]
        dvpdr = mer.vals[:, :, :, mer.lut[12], j]
            1./rr_2d*drad(rr_2d*vtbr_mm, rr)*my_weight

        dvrdt = mer.vals[:, :, :, mer.lut[37], j]
        dvtdt = mer.vals[:, :, :, mer.lut[38], j]
        dvpdt = mer.vals[:, :, :, mer.lut[39], j]

        dvrdp = mer.vals[:, :, :, mer.lut[46], j]
        dvtdp = mer.vals[:, :, :, mer.lut[47], j]
        dvpdp = mer.vals[:, :, :, mer.lut[48], j]

        # full derivatives B
        dbrdr = mer.vals[:, :, :, mer.lut[810], j]
        dbtdr = mer.vals[:, :, :, mer.lut[811], j]
        dbpdr = mer.vals[:, :, :, mer.lut[812], j]

        dbrdt = mer.vals[:, :, :, mer.lut[837], j]
        dbtdt = mer.vals[:, :, :, mer.lut[838], j]
        dbpdt = mer.vals[:, :, :, mer.lut[839], j]

        dbrdp = mer.vals[:, :, :, mer.lut[846], j]
        dbtdp = mer.vals[:, :, :, mer.lut[847], j]
        dbpdp = mer.vals[:, :, :, mer.lut[848], j]

        # mean derivatives v
        dvrdr_m = a.vals[:, :, a.lut[10], j].reshape((1, nt, nr))
        dvtdr_m = a.vals[:, :, a.lut[11], j].reshape((1, nt, nr))
        dvpdr_m = a.vals[:, :, a.lut[12], j].reshape((1, nt, nr))

        dvrdt_m = a.vals[:, :, a.lut[37], j].reshape((1, nt, nr))
        dvtdt_m = a.vals[:, :, a.lut[38], j].reshape((1, nt, nr))
        dvpdt_m = a.vals[:, :, a.lut[39], j].reshape((1, nt, nr))

        dvrdp_m = a.vals[:, :, a.lut[46], j].reshape((1, nt, nr))
        dvtdp_m = a.vals[:, :, a.lut[47], j].reshape((1, nt, nr))
        dvpdp_m = a.vals[:, :, a.lut[48], j].reshape((1, nt, nr))

        # mean derivatives B
        dbrdr_m = a.vals[:, :, a.lut[810], j].reshape((1, nt, nr))
        dbtdr_m = a.vals[:, :, a.lut[811], j].reshape((1, nt, nr))
        dbpdr_m = a.vals[:, :, a.lut[812], j].reshape((1, nt, nr))

        dbrdt_m = a.vals[:, :, a.lut[837], j].reshape((1, nt, nr))
        dbtdt_m = a.vals[:, :, a.lut[838], j].reshape((1, nt, nr))
        dbpdt_m = a.vals[:, :, a.lut[839], j].reshape((1, nt, nr))

        dbrdp_m = a.vals[:, :, a.lut[846], j].reshape((1, nt, nr))
        dbtdp_m = a.vals[:, :, a.lut[847], j].reshape((1, nt, nr))
        dbpdp_m = a.vals[:, :, a.lut[848], j].reshape((1, nt, nr))

        # compute induction terms

        # full radial
        ind_off = 0
        my_vals[:, :, ind_off + 0] += np.mean(dvrdt*bt + vr*dbtdt, axis=0)*my_weight
        my_vals[:, :, ind_off + 1] += -np.mean(dvtdt*br + vt*dbrdt, axis=0)*my_weight
        my_vals[:, :, ind_off + 2] += -np.mean(dvpdp*br + vp*dbrdp, axis=0)*my_weight
        my_vals[:, :, ind_off + 3] += np.mean(dvrdp*bp + vr*dbpdp, axis=0)*my_weight
        my_vals[:, :, ind_off + 4] += cott/rr_2d*np.mean(vr*bt - vt*br, axis=0)*my_weight
        ind_off += 5

        # full theta
        my_vals[:, :, ind_off + 0] += np.mean(dvtdp*bp + vt*dbpdp, axis=0)*my_weight
        my_vals[:, :, ind_off + 1] += -np.mean(dvpdp*bt + vp*dbtdp, axis=0)*my_weight
        my_vals[:, :, ind_off + 2] += -np.mean(dvrdr*bt + vr*dbtdr, axis=0)*my_weight
        my_vals[:, :, ind_off + 3] += np.mean(dvtdr*br + vt*dbrdr, axis=0)*my_weight
        my_vals[:, :, ind_off + 4] += 1/rr_2d*np.mean(vt*br - vr*bt, axis=0)*my_weight
        ind_off += 5

        # full phi
        my_vals[:, :, ind_off + 0] += np.mean(dvpdr*br + vp*dbrdr, axis=0)*my_weight
        my_vals[:, :, ind_off + 1] += -np.mean(dvrdr*bp + vr*dbpdr, axis=0)*my_weight
        my_vals[:, :, ind_off + 2] += -np.mean(dvtdt*bp + vt*dbpdt, axis=0)*my_weight
        my_vals[:, :, ind_off + 3] += np.mean(dvpdt*bt + vp*dbtdt, axis=0)*my_weight
        my_vals[:, :, ind_off + 4] += 1/rr_2d*np.mean(vp*br - vr*bp, axis=0)*my_weight
        ind_off += 5

        # mean radial
        my_vals[:, :, ind_off + 0] += np.mean(dvrdt_m*bt_m + vr*dbtdt_m, axis=0)*my_weight
        my_vals[:, :, ind_off + 1] += -np.mean(dvtdt_m*br_m + vt_m*dbrdt_m, axis=0)*my_weight
        my_vals[:, :, ind_off + 2] += -np.mean(dvpdp_m*br_m + vp_m*dbrdp_m, axis=0)*my_weight
        my_vals[:, :, ind_off + 3] += np.mean(dvrdp_m*bp_m + vr_m*dbpdp_m, axis=0)*my_weight
        my_vals[:, :, ind_off + 4] += cott/rr_2d*np.mean(vr_m*bt_m - vt_m*br_m, axis=0)*my_weight
        ind_off += 5

        # mean theta
        my_vals[:, :, ind_off + 0] += np.mean(dvtdp_m*bp_m + vt_m*dbpdp_m, axis=0)*my_weight
        my_vals[:, :, ind_off + 1] += -np.mean(dvpdp_m*bt_m + vp_m*dbtdp_m, axis=0)*my_weight
        my_vals[:, :, ind_off + 2] += -np.mean(dvrdr_m*bt_m + vr_m*dbtdr_m, axis=0)*my_weight
        my_vals[:, :, ind_off + 3] += np.mean(dvtdr_m*br_m + vt_m*dbrdr_m, axis=0)*my_weight
        my_vals[:, :, ind_off + 4] += 1/rr_2d*np.mean(vt_m*br_m - vr_m*bt_m, axis=0)*my_weight
        ind_off += 5

        # mean phi
        my_vals[:, :, ind_off + 0] += np.mean(dvpdr_m*br_m + vp_m*dbrdr_m, axis=0)*my_weight
        my_vals[:, :, ind_off + 1] += -np.mean(dvrdr_m*bp_m + vr_m*dbpdr_m, axis=0)*my_weight
        my_vals[:, :, ind_off + 2] += -np.mean(dvtdt_m*bp_m + vt_m*dbpdt_m, axis=0)*my_weight
        my_vals[:, :, ind_off + 3] += np.mean(dvpdt_m*bt_m + vp_m*dbtdt_m, axis=0)*my_weight
        my_vals[:, :, ind_off + 4] += 1/rr_2d*np.mean(vp_m*br_m - vr_m*bp_m, axis=0)*my_weight
        ind_off += 5

        ind_off_full = 0
        ind_off_mean = 15

        # fluc terms
        for k in range(15):
            my_vals[:, :, ind_off + k] = my_vals[:, :, ind_off_full + k] -\
                my_vals[:, :, ind_off_mean + k]

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
    print (format_time(t2 - t1))
    print(make_bold(fill_str('total time', lent, char)), end='')
    print (make_bold(format_time(t2 - t1_glob)))
    print ('data saved at ')
    print (make_bold(savefile))
    print ('data saved at ')
    print (savefile)
