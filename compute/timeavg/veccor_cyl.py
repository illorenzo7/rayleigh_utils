##################################################################
# Routine to average Rayleigh data in time (generic)
# Author: Loren Matilsky
# Created: 05/28/2021
##################################################################
# This routine computes the time average of cylindrical vector components
# (lambda, phi, z) of either u, om, b, or j (specify --type)
# Rayleigh data (AZ_Avgs, Meridional Slices), breaking up terms 
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
    # get the name of the run directory + CLAs
    args = sys.argv
    clas0, clas = read_clas(args)
    dirname = clas0['dirname']

    # overwrite defaults
    kw_default = dict({'type': 'u', 'ntheta': None})
    kw = update_dict(kw_default, clas)
    vectype = kw.type

    # Get the Rayleigh data directory
    radatadir1 = dirname + '/' + dataname1 + '/'
    radatadir2 = dirname + '/' + dataname2 + '/'

    # Get all the file names in datadir and their integer counterparts
    # Use mer slices, since there are sometimes fewer of them than azavgs
    file_list, int_file_list, nfiles = get_file_lists(radatadir1, clas)

    # Get the problem size
    nproc_min, nproc_max, n_per_proc_min, n_per_proc_max =\
            opt_workload(nfiles, nproc)

    # Will need the first data file for a number of things
    a0 = reading_func1(radatadir1 + file_list[0], '')
    # Make sure you put same number of slices in AZ_Avgs and 
    # Meridional_Slices, or this won't work

    # update ntheta unless told otherwise
    if kw.ntheta is None:
        kw.ntheta, dummy, dummy, dummy = a0.vals.shape

    # get grid information
    di_grid = get_grid_info(dirname, ntheta=kw.ntheta)
    nt = di_grid['nt']
    nr = di_grid['nr']
    tt = di_grid['tt']
    rr = di_grid['rr']
    rr_3d = di_grid['rr_3d']
    cost_2d = di_grid['cost_2d']
    sint_2d = di_grid['sint_2d']
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
    meta = [dirname, radatadir1, radatadir2, tt, rr, nt, nr, cost_2d, sint_2d, rr_3d, cost_3d, sint_3d, nfiles, vectype]
else:
    meta = None
dirname, radatadir1, radatadir2, tt, rr, nt, nr, cost_2d, sint_2d, rr_3d, cost_3d, sint_3d, nfiles, vectype = comm.bcast(meta, root=0)

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
nq = 6*3
my_vals = np.zeros((nt, nr, nq))

for i in range(my_nfiles):
    a = reading_func1(radatadir1 + str(my_files[i]).zfill(8), '')
    mer = reading_func2(radatadir2 + str(my_files[i]).zfill(8), '')

    if vectype == 'u':
        iqr, iqt, iqp = 1, 2, 3
    elif vectype == 'om':
        iqr, iqt, iqp = 301, 302, 303
    elif vectype == 'b':
        iqr, iqt, iqp = 801, 802, 803
    elif vectype == 'j':
        iqr, iqt, iqp = 1001, 1004, 1007
    else:
        if rank==0:
            print("you specified --type " + vectype)
            print("but type must be 'u', 'om', 'b', or 'j'. Exiting...")
        sys.exit(1)

    # take mean along the time axis;
    niter = min(a.niter, mer.niter)
    my_weight = 1.0/(nfiles*niter)
    for j in range(niter -1, -1, -1): # go last to first in case "niters" 
        # full vector
        vecr = mer.vals[:, :, :, mer.lut[iqr], j]
        vect = mer.vals[:, :, :, mer.lut[iqt], j]
        vecp = mer.vals[:, :, :, mer.lut[iqp], j]

        # mean vector
        vecr_m = a.vals[:, :, a.lut[iqr], j].reshape((1, nt, nr))
        vect_m = a.vals[:, :, a.lut[iqt], j].reshape((1, nt, nr))
        vecp_m = a.vals[:, :, a.lut[iqp], j].reshape((1, nt, nr))

        # compute cylindrical components
        vecl = vecr*sint_3d + vect*cost_3d
        vecz = vecr*cost_3d - vect*sint_3d

        vecl_m = vecr_m*sint_2d + vect_m*cost_2d
        vecz_m = vecr_m*cost_2d - vect_m*sint_2d


        # compute vector amplitudes

        # full 
        ind_off = 0
        # diagonal pieces
        my_vals[:, :, ind_off + 0] += np.mean(vecl**2, axis=0)*my_weight
        my_vals[:, :, ind_off + 1] += np.mean(vecp**2, axis=0)*my_weight
        my_vals[:, :, ind_off + 2] += np.mean(vecz**2, axis=0)*my_weight
        # non-diag
        my_vals[:, :, ind_off + 3] += np.mean(vecl*vecp, axis=0)*my_weight
        my_vals[:, :, ind_off + 4] += np.mean(vecl*vecz, axis=0)*my_weight 
        my_vals[:, :, ind_off + 5] += np.mean(vecp*vecz, axis=0)*my_weight

        # mean 
        ind_off += 6
        # diagonal pieces
        my_vals[:, :, ind_off + 0] += np.mean(vecl_m**2, axis=0)*my_weight
        my_vals[:, :, ind_off + 1] += np.mean(vecp_m**2, axis=0)*my_weight
        my_vals[:, :, ind_off + 2] += np.mean(vecz_m**2, axis=0)*my_weight
        # non-diag
        my_vals[:, :, ind_off + 3] += np.mean(vecl_m*vecp_m, axis=0)*my_weight
        my_vals[:, :, ind_off + 4] += np.mean(vecl_m*vecz_m, axis=0)*my_weight 
        my_vals[:, :, ind_off + 5] += np.mean(vecp_m*vecz_m, axis=0)*my_weight

        if rank == 0:
            pcnt_done = (i + 1)/my_nfiles*100.
            print(fill_str('computing', lent, char) +\
                    ('rank 0 %5.1f%% done' %pcnt_done), end='\r')

# easily compute the flucuating terms
# fluc terms
for k in range(6):
    my_vals[:, :, 12 + k] = my_vals[:, :, k] - my_vals[:, :, 6 + k]


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
    savename =  vectype + 'corcyl-' +\
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
