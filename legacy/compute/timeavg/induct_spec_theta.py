##################################################################
# Routine to average Rayleigh spectra, multiplying induction 
# by B* (this gets time rate of change of various spectral components)
# Author: Loren Matilsky
# Created: 12/12/2021
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

# import the reading routines
from rayleigh_diagnostics import AZ_Avgs, Shell_Avgs, G_Avgs,\
        Shell_Spectra, Shell_Slices, Meridional_Slices, Equatorial_Slices

# Get clas
if rank == 0:
    # get the CLAs
    args = sys.argv
    clas0, clas = read_clas(args)
    # modules needed only by proc 0 
    import pickle

# Get the desired datatype (specav)
radtype = 'specav'
reading_func = Shell_Spectra
dataname = 'Shell_Spectra'

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
    dirname = clas0['dirname']

    # Get the Rayleigh data directory
    radatadir = dirname + '/' + dataname + '/'

    # Get all the file names in datadir and their integer counterparts
    file_list, int_file_list, nfiles = get_file_lists(radatadir, args)
    weight = 1.0/nfiles # this is the averaging weight

    # Get the problem size
    nproc_min, nproc_max, n_per_proc_min, n_per_proc_max =\
            opt_workload(nfiles, nproc)

    # Will need the first data file for array shape
    # (exclude the last axis---time axis---from the shape
    a0 = reading_func(radatadir + file_list[0], '')
    shape = np.shape(a0.vals[..., 0]) #

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

        # send  my_files nproc > 1
        if k >= 1:
            comm.send(my_files, dest=k)
else: # recieve my_files
    my_files = comm.recv(source=0)

# Broadcast dirname, radatadir, nq, etc.
if rank == 0:
    meta = [dirname, radatadir, shape, weight]
else:
    meta = None
dirname, radatadir, shape, weight = comm.bcast(meta, root=0)

# Checkpoint and time
comm.Barrier()
if rank == 0:
    t2 = time.time()
    print (format_time(t2 - t1))
    print ('Considering %i %s files for the average: %s through %s'\
        %(nfiles, dataname, file_list[0], file_list[-1]))
    print(fill_str('computing', lent, char), end='\r')
    t1 = time.time()

# Now analyze the data
my_vals = np.zeros(shape)
my_nfiles = len(my_files)
for i in range(my_nfiles):
    if rank == 0 and i == 0:
        a = a0
    else:   
        a = reading_func(radatadir + str(my_files[i]).zfill(8), '')

    # loop over the induction terms, multiplying by B*
    nell, nm, nr, nq, niter = np.shape(a.vals)
    bfield = np.conjugate(a.vals[:,:,:,a.lut[802],:])
    induct_set = [3002,3005,3008,3011,3102,3105,3108,3111,3114,3117,3120, 2912]
    for qval in induct_set:
        # take mean along the time axis, which is usually the last axis
        my_vals[..., a.lut[qval]] += np.mean(np.real(a.vals[..., a.lut[qval], :]*bfield), axis=3)*weight

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
    vals = np.zeros(shape)

    # Gather the results into this "master" array
    for j in range(nproc):
        if j >= 1:
            # Get my_vals from rank j
            my_vals = comm.recv(source=j, tag=0)
        # "my_vals" are all weighted: their sum equals the overall average
        vals += my_vals 

else: # other processes send their data
    comm.send(my_vals, dest=0, tag=0)

# Make sure proc 0 collects all data
comm.Barrier()

# proc 0 saves the data
if rank == 0:
    # create data directory if it doesn't already exist
    datadir = clas0['datadir']
    if datadir is None:
        datadir = dirname + '/data/'
    if not os.path.isdir(datadir):
        os.makedirs(datadir)

    # set the save name by what we are saving (dataname)
    # and first and last iteration files for the trace
    iter1, iter2 = int_file_list[0], int_file_list[-1]
    savename = 'induct_spec-' + str(iter1).zfill(8) + '_' +\
            str(iter2).zfill(8) + '.pkl'
    savefile = datadir + savename

    # save the data
    f = open(savefile, 'wb')
    di_sav = {'vals': vals, 'lut': a0.lut, 'qv': a0.qv}
    if radtype in ['specav', 'ssav']:
        di_sav['rvals'] = a0.radius/rsun
    pickle.dump(di_sav, f, protocol=4)
    f.close()
    t2 = time.time()
    print (format_time(t2 - t1))
    print(make_bold(fill_str('total time', lent, char)), end='')
    print (make_bold(format_time(t2 - t1_glob)))
    print ('data saved at ')
    print (make_bold(savefile))
