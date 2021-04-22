##################################################################
# Routine to trace Rayleigh Shell_Slices data in time and longitude
# Author: Loren Matilsky
# Created: 04/01/2021
##################################################################
# This routine computes the trace in time/longitude of quantities in the 
# Shell_Slices data for a particular simulation. 
# By default, the [qvals] are computed at each time, in a latitude strip
# defined by clat = [central latitude for average] and
# dlat = [range of latitudes to average over]) 
#
# The strip range can be changed using the options --clat and --dlat, e.g., 
# --clat 60 --dlat 30, for a strip averaged between 45 and 75 degrees 
# (North)
# The final datacube output ('vals') will have shape
# (ntimes, nphi, nr, nq), where nr and nq are the attributes of Shell_Slices
############################################################################

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
from rayleigh_diagnostics import Shell_Slices
reading_func = Shell_Slices
dataname = 'Shell_Slices'

if rank == 0:
    # modules needed only by proc 0 
    import pickle
    from rayleigh_diagnostics import GridInfo

# Checkpoint and time
comm.Barrier()
if rank == 0:
    t2 = time.time()
    print (format_time(t2 - t1))
    t1 = time.time()

# proc 0 reads the file lists and distributes them, also the meta data
if rank == 0:
    # read the arguments
    dirname = sys.argv[1]
    args = sys.argv[2:]
    clas = read_clas(dirname, args)
    nargs = len(args)

    # Get the Rayleigh data directory
    radatadir = dirname + '/' + dataname + '/'

    # Get all the file names in datadir and their integer counterparts
    file_list, int_file_list, nfiles = get_file_lists(radatadir, args)

    # Set other defaults
    clat = 10.
    dlat = 0. # by default do NOT average over latitude

    # reset parameters with CLAs
    for i in range(nargs):
        arg = args[i]
        if arg == '--clat':
            clat = float(args[i+1])
        if arg == '--dlat':
            dlat = float(args[i+1])

    # Get metadata
    a0 = reading_func(radatadir + file_list[0], '')
    qvals = clas['qvals']
    nq = len(qvals)
    rvals = clas['rvals']
    if rvals is None:
        rvals = a0.radius/rsun
    nrvals = len(rvals)
    irvals = []
    for rval in rvals:
        irvals.append(np.argmin(np.abs(a0.radius/rsun - rval)))
    irvals = np.array(irvals)
    # recompute the actual sample values we get
    rvals = a0.radius[irvals]/rsun

    # didn't put this earlier because it got messed up by the message
    # saying which depths we were taking
    print(fill_str('proc 0 distributing the file lists', lent, char),\
            end='')
 
    # Desired latitude range + nphi
    di_grid = get_grid_info(dirname)
    ith1 = np.argmin(np.abs(di_grid['tt_lat'] - (clat - dlat/2.)))
    ith2 = np.argmin(np.abs(di_grid['tt_lat'] - (clat + dlat/2.)))
    tw_strip = di_grid['tw'][ith1:ith2+1]
    tw_strip /= np.sum(tw_strip)
    tw_strip = tw_strip.reshape((1, ith2 - ith1 + 1, 1, 1))
    nphi =  di_grid['nphi']

    # Get the problem size
    nproc_min, nproc_max, n_per_proc_min, n_per_proc_max =\
            opt_workload(nfiles, nproc)

    # Get metadata
    nq = len(qvals)

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

# Broadcast meta data
if rank == 0:
    meta = [dirname, radatadir, qvals, nq, irvals, nrvals, nphi,\
            ith1, ith2, tw_strip]
else:
    meta = None
dirname, radatadir, qvals, nq, irvals, nrvals, nphi,\
        ith1, ith2, tw_strip = comm.bcast(meta, root=0)

# Checkpoint and time
comm.Barrier()
if rank == 0:
    t2 = time.time()
    print (format_time(t2 - t1))
    print ('Considering %i %s files for the trace: %s through %s'\
        %(nfiles, dataname, file_list[0], file_list[-1]))
    print ("clat = %+.1f" %clat)
    if ith2 - ith1 + 1 > 1:
        print("dlat = %+.2f in latitude" %dlat)
    else:
        print("dlat = 0.00 (not averaging in latitude)")
    print ("lat1 = %+.1f" %(clat - dlat/2.0))
    print ("lat2 = %+.1f" %(clat + dlat/2.0))
    print ("sampling values:")
    print ("qvals = " + arr_to_str(qvals, "%i"))
    print ("sampling locations:")
    print ("rvals = " + arr_to_str(rvals, '%1.3f'))
    print(fill_str('computing', lent, char), end='\r')
    t1 = time.time()

# Now analyze the data
my_times = []
my_iters = []
my_vals = []

for i in range(my_nfiles):
    a = reading_func(radatadir + str(my_files[i]).zfill(8), '')
    for j in range(a.niter):
        # space in the arrays
        # get desired values in the strip, careful not to change 
        # the dimensionality of the array: (nphi, nlats, nrvals, nq)
        vals_strip = a.vals[:, ith1:ith2+1, :, :, j]
        vals_strip = vals_strip[:, :, irvals, :]
        vals_strip = vals_strip[:, :, :, a.lut[qvals]]
        my_vals.append(np.sum(vals_strip*tw_strip, axis=1))
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
            # Get my_times, my_iters, my_vals from rank j
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
    # and first and last iteration files for the trace (and optional tag)
    if clat >= 0.:
        hemisphere = 'N'
    else:
        hemisphere = 'S'
        
    savename = 'time_lon' + clas['tag'] + '-' +\
            ('clat%s%02.0f_dlat%03.0f' %(hemisphere, np.abs(clat), dlat)) +\
            '-' + file_list[0] + '_' + file_list[-1] + '.pkl'
    savefile = datadir + savename    

    # save the data
    f = open(savefile, 'wb')
    # convert everything to arrays
    vals = np.array(vals)
    times = np.array(times)
    iters = np.array(iters)
    pickle.dump({'vals': vals, 'times': times, 'iters': iters, 'rvals': rvals, 'qvals': qvals, 'clat': clat, 'dlat': dlat}, f, protocol=4)
    f.close()
    t2 = time.time()
    print (format_time(t2 - t1))
    print(make_bold(fill_str('total time', lent, char)), end='')
    print (make_bold(format_time(t2 - t1_glob)))
    print ('data saved at ')
    print (make_bold(savefile))
