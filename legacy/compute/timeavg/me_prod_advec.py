##################################################################
# Author: Loren Matilsky
# Created: 01/16/2021
##################################################################
# This routine computes the time average of mag. energy production terms
# breaking down by:
# nonzero terms in triple product: tot, pmp, ppm, mmm, mpp, ppp
# production of B_[r, theta, phi]
# adv components: (B_r, B_theta, B_phi) * (derivs v), 2 curvature terms
# 6 x 3 x 5 = 90 terms
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
    from cla_util import *
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
    # get the name of the run directory + CLAs
    args = sys.argv
    clas0, clas = read_clas(args)
    dirname = clas0['dirname']

    # Get the Rayleigh data directory
    radatadir1 = dirname + '/' + dataname1 + '/'
    radatadir2 = dirname + '/' + dataname2 + '/'

    # Get all the file names in datadir and their integer counterparts
    file_list, int_file_list, nfiles = get_file_lists(radatadir1, args)

    # Get the problem size
    nproc_min, nproc_max, n_per_proc_min, n_per_proc_max =\
            opt_workload(nfiles, nproc)

    # get grid information
    di_grid = get_grid_info(dirname)
    nt = di_grid['nt']
    nr = di_grid['nr']
    rr = di_grid['rr']
    tt = di_grid['tt']
    rr_2d = di_grid['rr_2d']
    cott_2d = di_grid['cott_2d']
    rr_3d = di_grid['rr_3d']
    cott_3d = di_grid['cott_3d']

    # need density derivative (for dvp/dP)
    eq = get_eq(dirname)
    dlnrho = eq.dlnrho.reshape((1, 1, nr))

    # Distribute file lists to each process
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
    meta = [dirname, radatadir1, radatadir2, nt, nr, rr, tt,\
            rr_2d,  cott_2d, rr_3d, cott_3d, dlnrho, nfiles]
else:
    meta = None
dirname, radatadir1, radatadir2, nt, nr, rr, tt, rr_2d,\
    cott_2d, rr_3d, cott_3d, dlnrho, nfiles = comm.bcast(meta, root=0)

# Checkpoint and time
comm.Barrier()
if rank == 0:
    t2 = time.time()
    print (format_time(t2 - t1))
    print ('Considering %i %s/%s files for the trace: %s through %s'\
        %(nfiles, dataname1, dataname2, file_list[0], file_list[-1]))
    print(fill_str('computing', lent, char), end='\r')
    t1 = time.time()

# Now analyze the data
nq = 90
my_vals = np.zeros((nt, nr, nq))

my_nfiles = len(my_files)
bad_jt = False # first assume we have J_theta = 1004 (not 1002)
for i in range(my_nfiles):
    a = reading_func1(radatadir1 + str(my_files[i]).zfill(8), '')
    mer = reading_func2(radatadir2 + str(my_files[i]).zfill(8), '')

    # take mean along the time axis;
    niter = min(a.niter, mer.niter)
    my_weight = 1.0/(nfiles*niter)
    for j in range(niter -1, -1, -1): # go last to first in case "niters" 
        # don't agree
        # full v
        vr = mer.vals[:, :, :, mer.lut[1], j]
        vt = mer.vals[:, :, :, mer.lut[2], j]
        vp = mer.vals[:, :, :, mer.lut[3], j]

        # full B
        br = mer.vals[:, :, :, mer.lut[801], j]
        bt = mer.vals[:, :, :, mer.lut[802], j]
        bp = mer.vals[:, :, :, mer.lut[803], j]

        # full poloidal vort. and current (for phi-derivatives)
        omr = mer.vals[:, :, :, mer.lut[301], j]
        omt = mer.vals[:, :, :, mer.lut[302], j]
        jr = mer.vals[:, :, :, mer.lut[1001], j]
        try: # had 1002 for a while...which was bad :-(
            jt = mer.vals[:, :, :, mer.lut[1004], j]
        except:
            jt = np.zeros_like(jr)
            bad_jt = True

        # mean v
        vr_m = a.vals[:, :, a.lut[1], j].reshape((1, nt, nr))
        vt_m = a.vals[:, :, a.lut[2], j].reshape((1, nt, nr))
        vp_m = a.vals[:, :, a.lut[3], j].reshape((1, nt, nr))

        # mean B
        br_m = a.vals[:, :, a.lut[801], j].reshape((1, nt, nr))
        bt_m = a.vals[:, :, a.lut[802], j].reshape((1, nt, nr))
        bp_m = a.vals[:, :, a.lut[803], j].reshape((1, nt, nr))

        # mean poloidal vort. and current (for phi-derivatives)
        omr_m = a.vals[:, :, a.lut[301], j].reshape((1, nt, nr))
        omt_m = a.vals[:, :, a.lut[302], j].reshape((1, nt, nr))
        jr_m = a.vals[:, :, a.lut[1001], j].reshape((1, nt, nr))
        if not bad_jt:
            jt_m = a.vals[:, :, a.lut[1004], j].reshape((1, nt, nr))
        else:
            jt_m = np.zeros_like(jr_m)

        # full (poloidal) derivatives v
        dvrdr = drad(vr, rr)
        dvtdr = drad(vt, rr)
        dvpdr = drad(vp, rr)

        dvrdt = dth(vr, tt)/rr_3d
        dvtdt = dth(vt, tt)/rr_3d
        dvpdt = dth(vp, tt)/rr_3d

        # full (toroidal) derivatives v
        dvrdp = dvpdr + (1./rr_3d)*vp + omt
        dvtdp = dvpdt + (cott_3d/rr_3d)*vp - omr
        dvpdp = -dlnrho*vr - dvrdr - (2./rr_3d)*vr - dvtdt -\
                (cott_3d/rr_3d)*vt

        # full (poloidal) derivatives B
        dbrdr = drad(br, rr)
        dbtdr = drad(bt, rr)
        dbpdr = drad(bp, rr)

        dbrdt = dth(br, tt)/rr_3d
        dbtdt = dth(bt, tt)/rr_3d
        dbpdt = dth(bp, tt)/rr_3d

        # full (toroidal) derivatives B
        dbrdp = dbpdr + (1./rr_3d)*bp + jt
        dbtdp = dbpdt + (cott_3d/rr_3d)*bp - jr
        dbpdp = -dbrdr - (2./rr_3d)*br - dbtdt - (cott_3d/rr_3d)*bt

        # mean (poloidal) derivatives v
        dvrdr_m = drad(vr_m, rr)
        dvtdr_m = drad(vt_m, rr)
        dvpdr_m = drad(vp_m, rr)

        dvrdt_m = dth(vr_m, tt)/rr_3d
        dvtdt_m = dth(vt_m, tt)/rr_3d
        dvpdt_m = dth(vp_m, tt)/rr_3d

        # mean (toroidal) derivatives v
        dvrdp_m = dvpdr_m + (1./rr_3d)*vp_m + omt_m
        dvtdp_m = dvpdt_m + (cott_3d/rr_3d)*vp_m - omr_m
        dvpdp_m = -dlnrho*vr_m - dvrdr_m - (2./rr_3d)*vr_m - dvtdt_m -\
                (cott_3d/rr_3d)*vt_m

        # mean (poloidal) derivatives B
        dbrdr_m = drad(br_m, rr)
        dbtdr_m = drad(bt_m, rr)
        dbpdr_m = drad(bp_m, rr)

        dbrdt_m = dth(br_m, tt)/rr_3d
        dbtdt_m = dth(bt_m, tt)/rr_3d
        dbpdt_m = dth(bp_m, tt)/rr_3d

        # mean (toroidal) derivatives B
        dbrdp_m = dbpdr_m + (1./rr_3d)*bp_m + jt_m
        dbtdp_m = dbpdt_m + (cott_3d/rr_3d)*bp_m - jr_m
        dbpdp_m = -dbrdr_m - (2./rr_3d)*br_m - dbtdt_m -\
                (cott_3d/rr_3d)*bt_m

        # make vectors/tensors: 
        # row indexes vector components, col indexes derivative direction
        vv = [vr, vt, vp]
        dvdx = [    [dvrdr, dvrdt, dvrdp],\
                    [dvtdr, dvtdt, dvtdp],\
                    [dvpdr, dvpdt, dvpdp] ]
        vv_m = [vr_m, vt_m, vp_m]
        dvdx_m = [  [dvrdr_m, dvrdt_m, dvrdp_m],\
                    [dvtdr_m, dvtdt_m, dvtdp_m],\
                    [dvpdr_m, dvpdt_m, dvpdp_m] ]

        bb = [br, bt, bp]
        dbdx = [    [dbrdr, dbrdt, dbrdp],\
                    [dbtdr, dbtdt, dbtdp],\
                    [dbpdr, dbpdt, dbpdp] ]
        bb_m = [br_m, bt_m, bp_m]
        dbdx_m = [  [dbrdr_m, dbrdt_m, dbrdp_m],\
                    [dbtdr_m, dbtdt_m, dbtdp_m],\
                    [dbpdr_m, dbpdt_m, dbpdp_m] ]

        # also get fluc vectors/tensors
        vv_p = [0, 0, 0]
        dvdx_p = [[0 for k in range(3)] for l in range(3)]
        bb_p = [0, 0, 0]
        dbdx_p = [[0 for k in range(3)] for l in range(3)]

        for k in range(3):
            vv_p[k] = vv[k] - vv_m[k]
            bb_p[k] = bb[k] - bb_m[k]
            for l in range(3):
                dvdx_p[k][l] = dvdx[k][l] - dvdx_m[k][l]
                dbdx_p[k][l] = dbdx[k][l] - dbdx_m[k][l]

        # get tot, pmp, ppm, mmm, mpp, ppp, energy terms
        # adv break up: 3 b dot v terms, + 2 curvature terms
        count = 0 
        for k in range(6): 
            # get appropriate fields for tot, pmp, ppm, mmm, mpp, ppp
            if k == 0:
                bb1_loc = bb.copy() # the guy getting dotted into 
                                    # inductive terms
                vv_loc = vv.copy()
                dvdx_loc = dvdx.copy()
                bb_loc = bb.copy()
                dbdx_loc = dbdx.copy()
            if k == 1:
                bb1_loc = bb_p.copy()
                vv_loc = vv_m.copy()
                dvdx_loc = dvdx_m.copy()
                bb_loc = bb_p.copy()
                dbdx_loc = dbdx_p.copy()
            if k == 2:
                bb1_loc = bb_p.copy()
                vv_loc = vv_p.copy()
                dvdx_loc = dvdx_p.copy()
                bb_loc = bb_m.copy()
                dbdx_loc = dbdx_m.copy()
            if k == 3:
                bb1_loc = bb_m.copy()
                vv_loc = vv_m.copy()
                dvdx_loc = dvdx_m.copy()
                bb_loc = bb_m.copy()
                dbdx_loc = dbdx_m.copy()
            if k == 4 or k == 5:
                vv_loc = vv_p.copy()
                dvdx_loc = dvdx_p.copy()
                bb_loc = bb_p.copy()
                dbdx_loc = dbdx_p.copy()
            if k == 4:
                bb1_loc = bb_m.copy()
            if k == 5:
                bb1_loc = bb_p.copy()

            # loop over r, theta, phi
            for l in range(3):
                adv1 = -vv_loc[0]*dbdx_loc[l][0]
                adv2 = -vv_loc[1]*dbdx_loc[l][1]
                adv3 = -vv_loc[2]*dbdx_loc[l][2]
                # need curvature terms for adv + adv
                if l == 0:
                    curv1 = (1./rr_3d)*vv_loc[1]*bb_loc[1]
                    curv2 = (1./rr_3d)*vv_loc[2]*bb_loc[2]
                if l == 1:
                    curv1 = -(1./rr_3d)*vv_loc[1]*bb_loc[0]
                    curv2 = (1./rr_3d)*cott_3d*vv_loc[2]*bb_loc[2]
                if l == 2:
                    curv1 = -(1./rr_3d)*vv_loc[2]*bb_loc[0]
                    curv2 = -(1./rr_3d)*cott_3d*vv_loc[2]*bb_loc[1]

                my_vals[:, :, count] += np.mean(bb1_loc[l]*adv1, axis=0)*my_weight
                count += 1
                my_vals[:, :, count] += np.mean(bb1_loc[l]*adv2, axis=0)*my_weight
                count += 1
                my_vals[:, :, count] += np.mean(bb1_loc[l]*adv3, axis=0)*my_weight
                count += 1
                my_vals[:, :, count] += np.mean(bb1_loc[l]*curv1, axis=0)*my_weight
                count += 1
                my_vals[:, :, count] += np.mean(bb1_loc[l]*curv2, axis=0)*my_weight
                count += 1

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
            # Get vals from rank j
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

    # Set the savename what we are saving,
    # and first and last iteration files
    savename = 'me_prod_advec-' +\
            file_list[0] + '_' + file_list[-1] + '.pkl'
    savefile = datadir + savename

    # save the data
    # Get first and last iters of files
    f = open(savefile, 'wb')
    pickle.dump({'vals': vals}, f, protocol=4)
    f.close()
    t2 = time.time()
    print (format_time(t2 - t1))
    print(make_bold(fill_str('total time', lent, char)), end='')
    print (make_bold(format_time(t2 - t1_glob)))
    if bad_jt:
        print (buff_line)
        print ("J_theta = 1004 was not given and so was set to zero! ")
        print("dBr/dphi is WRONG!!!")
        print ("so radial advection and total induction terms are WRONG!!!")
        print (buff_line)
    print ('data saved at ')
    print (make_bold(savefile))
