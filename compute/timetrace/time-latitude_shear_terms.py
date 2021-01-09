# Routine to trace AZ_Avgs in time/latitude space (pick different radii)
# Created by: Loren Matilsky
# On: 01/05/2020
##################################################################
# This routine computes the trace in time/latitude of derived nonlinear
# quantities (shear terms appearing in the induction equation) for 
# Meridional_Slices data for a particular simulation. 
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

# Derivative functions needed by everyone
def drad(arr, radius):
    nphi, nt, nr = np.shape(arr)
    two_dr = np.zeros((1, 1, nr-2))
    two_dr[0, 0, :] = radius[:nr-2] - radius[2:nr]     
    deriv = np.zeros_like(arr)
    deriv[:, :, 1:nr-1] = (arr[:, :, :nr-2] - arr[:, :, 2:nr])/two_dr
    deriv[:, :, 0] = deriv[:, :, 1]
    deriv[:, :, nr-1] = deriv[:, :, nr-2]
    return deriv

def dth(arr, theta):
    nphi, nt, nr = np.shape(arr)
    two_dt = np.zeros((1, nt-2, 1))
    two_dt[0, :, 0] = theta[:nt-2] - theta[2:nt]     
    deriv = np.zeros_like(arr)
    deriv[:, 1:nt-1, :] = (arr[:, :nt-2, :] - arr[:, 2:nt, :])/two_dt
    deriv[:, 0, :] = deriv[:, 1, :]
    deriv[:, nt-1, :] = deriv[:, nt-2, :]
    return deriv

def dph(arr):
    nphi, nt, nr = np.shape(arr)
    dphi = 2.*np.pi/nphi
    arr2 = np.roll(arr, -1, axis=0)
    arr1 = np.roll(arr, 1, axis=0)
    deriv = (arr2 - arr1)/2./dphi
    return deriv

# modules needed by everyone
import numpy as np
# data type and reading function
import sys, os
sys.path.append(os.environ['rapp'])
from rayleigh_diagnostics import AZ_Avgs, Meridional_Slices
reading_func1 = AZ_Avgs
reading_func2 = Meridional_Slices
dataname1 = 'AZ_Avgs'
dataname2 = 'Meridional_Slices'
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
    radatadir1 = dirname + '/' + dataname1 + '/'
    radatadir2 = dirname + '/' + dataname2 + '/'

    # Get all the file names in datadir and their integer counterparts
    file_list, int_file_list, nfiles = get_file_lists(radatadir2)

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
    depths = [0.05, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 0.95]
    phi_deriv = False
    tag = ''
    for i in range(nargs):
        arg = args[i]
        if arg == '-depths':
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
        elif arg == '-dp': # include terms related to phi-derivatives
            phi_deriv = True

    # Get desired file list from command-line arguments
    args = sys.argv[2:]
    nargs = len(args)
    the_tuple = get_desired_range(int_file_list, args)
    if the_tuple is None:
        index_first, index_last = nfiles - 101, nfiles - 1  
        # By default trace over the last 100 files
    else:
        index_first, index_last = the_tuple

    # Remove parts of file lists we don't need
    file_list = file_list[index_first:index_last + 1]
    int_file_list = int_file_list[index_first:index_last + 1]
    nfiles = index_last - index_first + 1

    # Get the problem size
    nproc_min, nproc_max, n_per_proc_min, n_per_proc_max =\
            opt_workload(nfiles, nproc)

    # Get metadata
    a0 = reading_func1(radatadir1 + file_list[0], '')
    nrec_full = a0.niter

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

    # get r-indices associated with depths
    rinds = []
    for depth in depths:
        rinds.append(np.argmin(np.abs(rr_depth - depth)))

    # compute some derivative quantities for the grid
    cost_2d = cost.reshape((nt, 1))
    sint_2d = sint.reshape((nt, 1))
    rr_2d = rr.reshape((1, nr))
    rvals_2d = rr_2d[:, rinds]
    xx = rr_2d*sint_2d
    zz = rr_2d*cost_2d

    # Will need nrec (last niter) to get proper time axis size
    af = reading_func1(radatadir1 + file_list[-1], '')
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

# Broadcast meta data
if rank == 0:
    meta = [dirname, radatadir1, radatadir2, tt, rr, nt, nr, rvals_2d,\
        cost_2d, sint_2d, rinds, ndepths, phi_deriv]
else:
    meta = None
dirname, radatadir1, radatadir2, tt, rr, nt, nr, rvals_2d,\
    cost_2d, sint_2d, rinds, ndepths, phi_deriv = comm.bcast(meta, root=0)

# Checkpoint and time
comm.Barrier()
if rank == 0:
    t2 = time.time()
    print ('%8.2e s' %(t2 - t1))
    print ('Considering %i %s/%s files for the trace: %s through %s'\
        %(nfiles, dataname1, dataname2, file_list[0], file_list[-1]))
    print ("ntimes for trace = %i" %ntimes)
    print(fill_str('computing', lent, char), end='\r')
    t1 = time.time()

# Now analyze the data
my_times = np.zeros(my_ntimes)
my_iters = np.zeros(my_ntimes, dtype='int')
nq = 33
if phi_deriv:
    nq += 3
my_vals = np.zeros((my_ntimes, nt, ndepths, nq))

# need some further grid quantities (that won't be saved in the final dict)
rr_3d = rr.reshape((1, 1, nr))
sint = np.sin(tt)
sint_3d = sint.reshape((1, nt, 1))

my_count = 0
for i in range(my_nfiles):
    a = reading_func1(radatadir1 + str(my_files[i]).zfill(8), '')
    mer = reading_func2(radatadir2 + str(my_files[i]).zfill(8), '')

    for j in range(a.niter):
        # mean B
        br_m = a.vals[:, rinds, a.lut[801], j].reshape((1, nt, ndepths))
        bt_m = a.vals[:, rinds, a.lut[802], j].reshape((1, nt, ndepths))
        bp_m = a.vals[:, rinds, a.lut[803], j].reshape((1, nt, ndepths))

        # mean v
        vr_m = a.vals[:, :, a.lut[1], j].reshape((1, nt, nr))
        vt_m = a.vals[:, :, a.lut[2], j].reshape((1, nt, nr))
        vp_m = a.vals[:, :, a.lut[3], j].reshape((1, nt, nr))

        # mean v derivs
        dvrdr_m = drad(vr_m, rr)
        dvrdt_m = dth(vr_m, tt)/rr_3d
        dvtdr_m = drad(vt_m, rr)
        dvtdt_m = dth(vt_m, tt)/rr_3d
        dvpdr_m = drad(vp_m, rr)
        dvpdt_m = dth(vp_m, tt)/rr_3d

        # now get only the radial indices we need
        vr_m = vr_m[:, :, rinds]
        vt_m = vt_m[:, :, rinds]
        vp_m = vp_m[:, :, rinds]

        dvrdr_m = dvrdr_m[:, :, rinds]
        dvrdt_m = dvrdt_m[:, :, rinds]
        dvtdr_m = dvtdr_m[:, :, rinds]
        dvtdt_m = dvtdt_m[:, :, rinds]
        dvpdr_m = dvpdr_m[:, :, rinds]
        dvpdt_m = dvpdt_m[:, :, rinds]

        # total shear
        shear_r = a.vals[:, rinds, a.lut[1601], j]
        shear_t = a.vals[:, rinds, a.lut[1606], j]
        shear_p = a.vals[:, rinds, a.lut[1611], j]

        # mean shear
        shear_mm_r = a.vals[:, rinds, a.lut[1616], j]
        shear_mm_t = a.vals[:, rinds, a.lut[1621], j]
        shear_mm_p = a.vals[:, rinds, a.lut[1626], j]

        # fluc shear
        shear_pp_r = shear_r - shear_mm_r
        shear_pp_t = shear_t - shear_mm_t
        shear_pp_p = shear_p - shear_mm_p

        # full B
        br = mer.vals[:, :, rinds, mer.lut[801], j]
        bt = mer.vals[:, :, rinds, mer.lut[802], j]
        bp = mer.vals[:, :, rinds, mer.lut[803], j]

        # full v
        vr = mer.vals[:, :, :, mer.lut[1], j]
        vt = mer.vals[:, :, :, mer.lut[2], j]
        vp = mer.vals[:, :, :, mer.lut[3], j]

        # full v derivs
        dvrdr = drad(vr, rr)
        dvrdt = dth(vr, tt)/rr_3d
        dvtdr = drad(vt, rr)
        dvtdt = dth(vt, tt)/rr_3d
        dvpdr = drad(vp, rr)
        dvpdt = dth(vp, tt)/rr_3d
        if phi_deriv:
            dvrdp = mer.vals[:, :, :, mer.lut[28], j]/rr_3d/sint_3d
            dvtdp = mer.vals[:, :, :, mer.lut[29], j]/rr_3d/sint_3d
            dvpdp = mer.vals[:, :, :, mer.lut[30], j]/rr_3d/sint_3d

        # now get only the radial indices we need
        vr = vr[:, :, rinds]
        vt = vt[:, :, rinds]
        vp = vp[:, :, rinds]

        dvrdr = dvrdr[:, :, rinds]
        dvrdt = dvrdt[:, :, rinds]
        dvtdr = dvtdr[:, :, rinds]
        dvtdt = dvtdt[:, :, rinds]
        dvpdr = dvpdr[:, :, rinds]
        dvpdt = dvpdt[:, :, rinds]
        if phi_deriv:
            dvrdp = dvrdp[:, :, rinds]
            dvtdp = dvtdp[:, :, rinds]
            dvpdp = dvpdp[:, :, rinds]

        # fluc B
        br_p = br - br_m
        bt_p = bt - bt_m
        bp_p = bp - bp_m

        # fluc v
        vr_p = vr - vr_m
        vt_p = vt - vt_m
        vp_p = vp - vp_m

        # fluc deriv v
        dvrdr_p = dvrdr - dvrdr_m
        dvrdt_p = dvrdt - dvrdt_m
        dvtdr_p = dvtdr - dvtdr_m
        dvtdt_p = dvtdt - dvtdt_m
        dvpdr_p = dvpdr - dvpdr_m
        dvpdt_p = dvpdt - dvpdt_m

        # compute 33 terms!
        shear_mm_r_1 = np.mean(br_m*dvrdr_m, axis=0)
        shear_mm_r_2 = np.mean(bt_m*dvrdt_m, axis=0)
        shear_mm_r_3 = -np.mean(bt_m*vt_m, axis=0)/rvals_2d
        shear_mm_r_4 = -np.mean(bp_m*vp_m, axis=0)/rvals_2d

        shear_pp_r_1 = np.mean(br_p*dvrdr_p, axis=0)
        shear_pp_r_2 = np.mean(bt_p*dvrdt_p, axis=0)
        shear_pp_r_4 = -np.mean(bt_p*vt_p, axis=0)/rvals_2d
        shear_pp_r_5 = -np.mean(bp_p*vp_p, axis=0)/rvals_2d
        shear_pp_r_3 = shear_pp_r -\
            (shear_pp_r_1 + shear_pp_r_2 + shear_pp_r_4 + shear_pp_r_5)
#        shear_pp_r_3 = shear_r -\
#            (shear_pp_r_1 + shear_pp_r_2 + shear_pp_r_4 + shear_pp_r_5) -\
#            (shear_mm_r_1 + shear_mm_r_2 + shear_mm_r_3 + shear_mm_r_4)

        shear_mm_t_1 = np.mean(br_m*dvtdr_m, axis=0)
        shear_mm_t_2 = np.mean(bt_m*dvtdt_m, axis=0)
        shear_mm_t_3 = np.mean(bt_m*vr_m, axis=0)/rvals_2d
        shear_mm_t_4 = -np.mean(bp_m*vp_m, axis=0)*cost_2d/rvals_2d/sint_2d

        shear_pp_t_1 = np.mean(br_p*dvtdr_p, axis=0)
        shear_pp_t_2 = np.mean(bt_p*dvtdt_p, axis=0)
        shear_pp_t_4 = np.mean(bt_p*vr_p, axis=0)/rvals_2d
        shear_pp_t_5 = -np.mean(bp_p*vp_p, axis=0)*cost_2d/rvals_2d/sint_2d
        shear_pp_t_3 = shear_pp_t -\
            (shear_pp_t_1 + shear_pp_t_2 + shear_pp_t_4 + shear_pp_t_5)
#        shear_pp_t_3 = shear_t -\
#            (shear_pp_t_1 + shear_pp_t_2 + shear_pp_t_4 + shear_pp_t_5) -\
#            (shear_mm_t_1 + shear_mm_t_2 + shear_mm_t_3 + shear_mm_t_4)

        shear_mm_p_1 = np.mean(br_m*dvpdr_m, axis=0)
        shear_mm_p_2 = np.mean(bt_m*dvpdt_m, axis=0)
        shear_mm_p_3 = np.mean(bp_m*vr_m, axis=0)/rvals_2d
        shear_mm_p_4 = np.mean(bp_m*vt_m, axis=0)*cost_2d/rvals_2d/sint_2d

        shear_pp_p_1 = np.mean(br_p*dvpdr_p, axis=0)
        shear_pp_p_2 = np.mean(bt_p*dvpdt_p, axis=0)
        shear_pp_p_4 = np.mean(bp_p*vr_p, axis=0)/rvals_2d
        shear_pp_p_5 = np.mean(bp_p*vt_p, axis=0)*cost_2d/rvals_2d/sint_2d
        shear_pp_p_3 = shear_pp_p -\
            (shear_pp_p_1 + shear_pp_p_2 + shear_pp_p_4 + shear_pp_p_5)
#        shear_pp_p_3 = shear_p -\
#            (shear_pp_p_1 + shear_pp_p_2 + shear_pp_p_4 + shear_pp_p_5) -\
#            (shear_mm_p_1 + shear_mm_p_2 + shear_mm_p_3 + shear_mm_p_4)

        if phi_deriv:
            shear_pp_r_6 = np.mean(bp_p*dvrdp, axis=0)
            shear_pp_t_6 = np.mean(bp_p*dvtdp, axis=0)
            shear_pp_p_6 = np.mean(bp_p*dvpdp, axis=0)

        if my_count < my_ntimes: # make sure we don't go over the allotted
            # space in the arrays
            my_times[my_count] = a.time[j] 
            my_iters[my_count] = a.iters[j]

            ind_off = 0

            my_vals[my_count, :, :, ind_off + 0] = shear_mm_r
            my_vals[my_count, :, :, ind_off + 1] = shear_mm_r_1
            my_vals[my_count, :, :, ind_off + 2] = shear_mm_r_2
            my_vals[my_count, :, :, ind_off + 3] = shear_mm_r_3
            my_vals[my_count, :, :, ind_off + 4] = shear_mm_r_4
            ind_off += 5

            my_vals[my_count, :, :, ind_off + 0] = shear_pp_r
            my_vals[my_count, :, :, ind_off + 1] = shear_pp_r_1
            my_vals[my_count, :, :, ind_off + 2] = shear_pp_r_2
            my_vals[my_count, :, :, ind_off + 3] = shear_pp_r_3
            my_vals[my_count, :, :, ind_off + 4] = shear_pp_r_4
            my_vals[my_count, :, :, ind_off + 5] = shear_pp_r_5
            if phi_deriv:
                my_vals[my_count, :, :, ind_off + 6] = shear_pp_r_6
                ind_off += 1
            ind_off += 6

            my_vals[my_count, :, :, ind_off + 0] = shear_mm_t
            my_vals[my_count, :, :, ind_off + 1] = shear_mm_t_1
            my_vals[my_count, :, :, ind_off + 2] = shear_mm_t_2
            my_vals[my_count, :, :, ind_off + 3] = shear_mm_t_3
            my_vals[my_count, :, :, ind_off + 4] = shear_mm_t_4
            ind_off += 5

            my_vals[my_count, :, :, ind_off + 0] = shear_pp_t
            my_vals[my_count, :, :, ind_off + 1] = shear_pp_t_1
            my_vals[my_count, :, :, ind_off + 2] = shear_pp_t_2
            my_vals[my_count, :, :, ind_off + 3] = shear_pp_t_3
            my_vals[my_count, :, :, ind_off + 4] = shear_pp_t_4
            my_vals[my_count, :, :, ind_off + 5] = shear_pp_t_5
            if phi_deriv:
                my_vals[my_count, :, :, ind_off + 6] = shear_pp_t_6
                ind_off += 1
            ind_off += 6

            my_vals[my_count, :, :, ind_off + 0] = shear_mm_p
            my_vals[my_count, :, :, ind_off + 1] = shear_mm_p_1
            my_vals[my_count, :, :, ind_off + 2] = shear_mm_p_2
            my_vals[my_count, :, :, ind_off + 3] = shear_mm_p_3
            my_vals[my_count, :, :, ind_off + 4] = shear_mm_p_4
            ind_off += 5

            my_vals[my_count, :, :, ind_off + 0] = shear_pp_p
            my_vals[my_count, :, :, ind_off + 1] = shear_pp_p_1
            my_vals[my_count, :, :, ind_off + 2] = shear_pp_p_2
            my_vals[my_count, :, :, ind_off + 3] = shear_pp_p_3
            my_vals[my_count, :, :, ind_off + 4] = shear_pp_p_4
            my_vals[my_count, :, :, ind_off + 5] = shear_pp_p_5
            if phi_deriv:
                my_vals[my_count, :, :, ind_off + 6] = shear_pp_p_6
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
    savename = dirname_stripped + '_time-latitude_shear_terms_' + tag +\
            file_list[0] + '_' + file_list[-1] + '.pkl'
    savefile = datadir + savename

    # save the data
    f = open(savefile, 'wb')
    # will also need first and last iteration files
    iter1, iter2 = int_file_list[0], int_file_list[-1]
    pickle.dump({'vals': vals, 'times': times, 'iters': iters,\
    'depths': depths,'rinds': rinds, 'niter': len(iters),\
    'ndepths': ndepths, 'nq': nq, 'iter1': iter1, 'iter2': iter2, 'rr': rr,\
    'rr_depth': rr_depth, 'rr_height': rr_height, 'nr': nr, 'ri': ri,\
    'ro': ro, 'd': d, 'tt': tt, 'tt_lat': tt_lat,\
    'sint': sint, 'cost': cost, 'ntheta': nt,\
    'rr_2d': rr_2d, 'sint_2d': sint_2d,\
    'cost_2d': cost_2d, 'xx': xx, 'zz': zz}, f, protocol=4)
    f.close()
    t2 = time.time()
    print ('%8.2e s' %(t2 - t1))
    print(fill_str('total time', lent, char), end='')
    print ('%8.2e s' %(t2 - t1_glob))
    print ('saving data at')
    print (savefile)
