##################################################################
# Routine to trace Rayleighdata in time
# Author: Loren Matilsky
# Created: 12/18/2018
# Parallelized: 11/26/2020
##################################################################
# This routine computes the trace in time of the values in various 
# Rayleigh data types
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
# import the reading routines
from rayleigh_diagnostics import AZ_Avgs, Shell_Avgs, G_Avgs, Point_Probes,\
        Shell_Spectra, Shell_Slices, Meridional_Slices, Equatorial_Slices

# Broadcast the desired datatype (azav by default)
if rank == 0:
    # get the CLAs
    args = sys.argv
    clas0, clas = read_clas(args)
    if 'radtype' in clas:
        radtype = clas['radtype']
    else:
        radtype = 'gav' # default
else:
    radtype = None
radtype = comm.bcast(radtype, root=0)

if radtype == 'azav':
    reading_func = AZ_Avgs
    dataname = 'AZ_Avgs'
if radtype == 'shav':
    reading_func = Shell_Avgs
    dataname = 'Shell_Avgs'
if radtype == 'gav':
    reading_func = G_Avgs
    dataname = 'G_Avgs'
if radtype == 'specav':
    reading_func = Shell_Spectra
    dataname = 'Shell_Spectra'
if radtype == 'ssav':
    reading_func = Shell_Slices
    dataname = 'Shell_Slices'
if radtype == 'merav':
    reading_func = Meridional_Slices
    dataname = 'Meridional_Slices'
if radtype == 'eqav':
    reading_func = Equatorial_Slices
    dataname = 'Equatorial_Slices'
if radtype == 'pp':
    reading_func = Point_Probes
    dataname = 'Point_Probes'

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

    # Distribute file_list and my_ntimes to each process
    for k in range(nproc - 1, -1, -1):
        # distribute the partial file list to other procs 
        if k < nproc_max: # first processes analyzes more files
            my_nfiles = np.copy(n_per_proc_max)
            istart = k*my_nfiles
            iend = istart + my_nfiles
        else: # last processes analyze fewer files
            my_nfiles = np.copy(n_per_proc_min)
            istart = nproc_max*n_per_proc_max + (k - nproc_max)*my_nfiles
            iend = istart + my_nfiles

        # Get the file list portion for rank k
        my_files = np.copy(int_file_list[istart:iend])

        # send appropriate file info if nproc > 1
        if k >= 1:
            comm.send([my_files, my_nfiles], dest=k)
else: # recieve appropriate file info if rank > 1
    my_files, my_nfiles = comm.recv(source=0)

# Broadcast metadata
if rank == 0:
    meta = [dirname, radatadir]
else:
    meta = None
dirname, radatadir = comm.bcast(meta, root=0)

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
        if radtype == 'gav':
            my_vals.append(a.vals[j, :])
        else:
            my_vals.append(a.vals[..., j])
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

    # Gather the results
    for j in range(nproc):
        if j >= 1:
            # Get my_ntimes, my_times, my_iters, my_vals from rank j
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

    # Set the timetrace savename 
    savename = dataname + '_trace-' + file_list[0] + '_' + file_list[-1] + '.pkl'
    savefile = datadir + savename

    # save the data
    f = open(savefile, 'wb')
    # convert everything to arrays
    vals = np.array(vals)
    times = np.array(times)
    iters = np.array(iters)
    di_sav = {'vals': vals, 'times': times, 'iters': iters, 'lut': a.lut, 'qv': a.qv}
    if radtype == 'pp':
        di_sav['rvals'] = a.radius
        di_sav['thetavals'] = np.arccos(a.costheta)
        di_sav['phivals'] = a.phi

    pickle.dump(di_sav, f, protocol=4)
    f.close()
    t2 = time.time()
    print (format_time(t2 - t1))
    print(make_bold(fill_str('total time', lent, char)), end='')
    print (make_bold(format_time(t2 - t1_glob)))
    print ('data saved at ')
    print (make_bold(savefile))
