##################################################################
# Routine to average Rayleigh data in time (generic)
# Author: Loren Matilsky
# Created: 04/08/2021
##################################################################
# This routine computes how long it takes to read Rayleigh data
# default data type is azav. to specify another 
# (specav, gav, shav, ssav, merav, eqav) use
# --radtype [radataname]
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
    print(fill_str('processes importing necessary modules', lent, char))

# modules needed by everyone (need to see the rayleigh_diagnostics
# module in post_processing)
import sys, os
sys.path.append(os.environ['rapp'])

# broadcast the desired datatype (azav by default)
if rank == 0:
    # get the CLAs
    args = sys.argv
    clas0, clas = read_clas(args)

    # overwrite defaults
    kw_default = dict({'radtype': 'azav'})
    kw = update_dict(kw_default, clas)
    radtype = kw.radtype

    # get reading function and dataname from the di_radtypes container
    reading_func = di_radtypes[radtype].reading_func
    dataname = di_radtypes[radtype].dataname
    meta = [radtype, reading_func, dataname]
else:
    meta = None
radtype, reading_func, dataname = comm.bcast(meta, root=0)


# proc 0 reads the file lists and distributes them
if rank == 0:
    # get dirname
    dirname = clas0['dirname']

    # Get the Rayleigh data directory
    radatadir = dirname + '/' + dataname + '/'

    # Get all the file names in datadir and their integer counterparts
    file_list, int_file_list, nfiles = get_file_lists(radatadir, clas)
    print (buff_line)
    print ('Considering %i %s files: %s through %s'\
        %(nfiles, dataname, file_list[0], file_list[-1]))
    print (buff_line)

    # Get the problem size
    nproc_min, nproc_max, n_per_proc_min, n_per_proc_max =\
            opt_workload(nfiles, nproc)

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
    meta = [dirname, radatadir]
else:
    meta = None
dirname, radatadir = comm.bcast(meta, root=0)

# now test if everything reads
my_nfiles = len(my_files)
for i in range(my_nfiles):
    try:
        a = reading_func(radatadir + str(my_files[i]).zfill(8), '')
    except:
        print ('rank = %i here' %rank)
        print ('could not read ' + radatadir + str(my_files[i]).zfill(8))

    if rank == 0:
        pcnt_done = (i + 1)/my_nfiles*100.
        print(fill_str('reading', lent, char) +\
                ('rank 0 %5.1f%% done' %pcnt_done), end='\r')

# Checkpoint and time
#comm.Barrier()
if rank == 0:
    t2 = time.time()
    print('\n' + fill_str('tot read time', lent, char), end='')
    print (format_time(t2 - t1))
    print (buff_line)
