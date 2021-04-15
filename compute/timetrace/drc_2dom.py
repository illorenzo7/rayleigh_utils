##################################################################
# Routine to trace DR contrast in time (for Tacho paper)
# Author: Loren Matilsky
# Created: 12/18/2018
# Parallelized: 11/26/2020
##################################################################
# This routine computes the trace in time of sigma_Omega (spherical stddev
# of Omega from equator to latcut = 60.)
# separates into full, CZ (r/rsun = 0.9 to outer surface), and RZ (full RZ)
# shape(vals) = (ntimes, 3) 
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
from rayleigh_diagnostics import AZ_Avgs
reading_func = AZ_Avgs
dataname = 'AZ_Avgs'

if rank == 0:
    # modules needed only by proc 0 
    import pickle
    from rayleigh_diagnostics import GridInfo

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
    # read the arguments
    dirname = sys.argv[1]
    args = sys.argv[2:]

    # Get the Rayleigh data directory
    radatadir = dirname + '/' + dataname + '/'

    # Get all the file names in datadir and their integer counterparts
    file_list, int_file_list, nfiles = get_file_lists(radatadir, args)

    # Get the problem size
    nproc_min, nproc_max, n_per_proc_min, n_per_proc_max =\
            opt_workload(nfiles, nproc)

    # get grid information
    gi = GridInfo(dirname + '/grid_info', '')
    rw = gi.rweights
    nr = gi.nr
    ir_bcz = get_parameter(dirname, 'ncheby')[1] - 1
    nr_cz = ir_bcz + 1
    nr_rz = nr - nr_cz

    # compute moment arm
    rr = gi.radius
    sint = gi.sintheta
    nt = gi.ntheta
    xx = rr.reshape((1, nr))*sint.reshape((nt, 1))
    xx_cz = xx[:, :nr_cz]
    xx_rz = xx[:, nr_cz:]

    # get latitudes, and upper bound (lat 60)
    cost = gi.costheta
    tt = np.arccos(cost)
    tt_lat = (np.pi/2.0 - tt)*180.0/np.pi
    latcut = 60.0
    icut1 = np.argmin(np.abs(tt_lat + latcut))
    icut2 = np.argmin(np.abs(tt_lat - latcut))
    tw = gi.tweights[icut1:icut2+1]
    tw /= np.sum(tw)
    ntcut = len(tw)
    tw = tw.reshape((ntcut, 1))

    # get averaging weights for CZ and RZ separately
    r0 = 0.9
    ir0 = np.argmin(np.abs(rr/rsun - r0))
    rw_cz = np.copy(rw[:ir0 + 1])
    rw_rz = np.copy(rw[nr_cz:])
    rw_cz /= np.sum(rw_cz)
    rw_rz /= np.sum(rw_rz)

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

        # send  my_files, my_nfiles if nproc > 1
        if k >= 1:
            comm.send([my_files, my_nfiles], dest=k)
else: # recieve my_files, my_nfiles
    my_files, my_nfiles = comm.recv(source=0)

# Broadcast dirname, radatadir, nq
if rank == 0:
    meta = [dirname, radatadir, nr, nr_cz, nr_rz, rw, rw_cz, rw_rz,\
            ir_bcz, tw, nt, xx, ir0, icut1, icut2]
else:
    meta = None
dirname, radatadir, nr, nr_cz, nr_rz, rw, rw_cz, rw_rz,\
        ir_bcz, tw, nt, xx, ir0, icut1, icut2 = comm.bcast(meta, root=0)

# Checkpoint and time
comm.Barrier()
if rank == 0:
    t2 = time.time()
    print (format_time(t2 - t1))
    print ('Considering %i %s files for the trace: %s through %s'\
        %(nfiles, dataname, file_list[0], file_list[-1]))
    print(fill_str('computing', lent, char), end='\r')
    t1 = time.time()

# Now analyze the data
my_times = []
my_iters = []
my_vals = []

for i in range(my_nfiles):
    a = reading_func(radatadir + str(my_files[i]).zfill(8), '')
    for j in range(a.niter):
        om_merplane = (np.copy(a.vals[:, :, a.lut[3], j])/xx)[icut1:icut2+1,:]
        om_mean = np.sum(tw*om_merplane, axis=0).reshape((1, nr))
        om2 = np.sum((om_merplane - om_mean)**2.*tw, axis=0)
        Dom_r = np.sqrt(om2)

        # Get the values in the CZ/RZ separately
        gav = np.sqrt(np.sum(rw*Dom_r))
        gav_cz = np.sqrt(np.sum(rw_cz*Dom_r[:ir0+1]))
        gav_rz = np.sqrt(np.sum(rw_rz*Dom_r[nr_cz:]))

        my_vals.append(np.array([gav, gav_cz, gav_rz]))

        my_times.append(a.time[j])
        my_iters.append(a.iters[j])
    if rank == 0:
        pcnt_done = i/my_nfiles*100.0
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
    # initialize empty lists for the final arrays
    times = []
    iters = []
    vals = []

    # Gather the results into these "master" arrays
    for j in range(nproc):
        if j >= 1:
            # Get data from rank j
            my_times, my_iters, my_vals = comm.recv(source=j)
        times += my_times
        iters += my_iters
        vals += my_vals
else: # other processes send their data
    comm.send([my_times, my_iters, my_vals], dest=0)

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
    savename = 'drc_2dom-' + file_list[0] + '_' + file_list[-1] + '.pkl'
    savefile = datadir + savename

    # save the data
    f = open(savefile, 'wb')
    # convert everything to arrays
    vals = np.array(vals)
    times = np.array(times)
    iters = np.array(iters)
    di_sav = {'vals': vals, 'times': times, 'iters': iters}
    pickle.dump(di_sav, f, protocol=4)
    f.close()
    t2 = time.time()
    print (format_time(t2 - t1))
    print(make_bold(fill_str('total time', lent, char)), end='')
    print (make_bold(format_time(t2 - t1_glob)))
    print ('data saved at ')
    print (make_bold(savefile))
