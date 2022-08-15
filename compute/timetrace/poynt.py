##################################################################
# Routine to trace Rayleigh data in time, in different quadrants
# of the meridional plane
# Author: Loren Matilsky
# Created: 08/15/2022
##################################################################
# This routine computes the trace in time of quantities relevant for 
# mag. energy production in a spherical shell
# needed: 1102,1103,1104) and 1916, 1436, 2001
# on shell_avgs
#
# By default the "quadrant" is the entire shell
#
# if --nquadr or --rbounds or --irbounds is specified 
# quantities are computed over multiple shells
# (Poynting flux is evaluated at the radial boundaries, other quantities
# are integrated between these boundaries
#

# vals is shape
# [ntimes, nq=7, nquadr+1] 
# quantities: total me, radial me, theta me, phi me,
# -mag_work, -joule_heat, radial Poynting flux
# first six quantities are integrated over quadrants (last column in the
# radial dimension is zero)
# poynting flux is evaluated at the boundaries (increasing inner to outer)
# shells are ordered in increasing radius (inner shell to outer shell) 
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
sys.path.append(os.environ['raco'])
sys.path.append(os.environ['rapp'])
from rayleigh_diagnostics import Shell_Avgs

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
    args = sys.argv
    clas0, clas = read_clas(args)
    dirname = clas0['dirname']
    tag = clas0['tag']

    # get grid information
    di_grid = get_grid_info(dirname)
    nr = di_grid['nr']
    rw = di_grid['rw']
    rr = di_grid['rr']
    rmin, rmax = np.min(rr)/rsun, np.max(rr)/rsun

    # get grid info + default kwargs
    kwargs_default = dict({})
    kwargs_default['nquadr'] = None # can divide up the radial grid into nquadr equally spaced domains
    kwargs_default['rbounds'] = None # can specify radial domain boundaries directly (units of rsun, e.g., 0.721 0.863 0.92)
    kwargs_default['irbounds'] = None # can specify radial domain boundaries directly (radial index, e.g., 32 64 96

    # update these possibly
    kw = update_dict(kwargs_default, clas)

    # deal w/ radial boundaries
    dataname = 'Shell_Avgs'
    if not kw.nquadr is None: # equally spaced domain boundaries (not the default)
        kw.rbounds = np.linspace(rmax, rmin, kw.nquadr + 1) # remember: rr is DECREASING
    if kw.irbounds is None: # this is the default
        if kw.rbounds is None: # this is the default
            irbounds = [nr - 1, 0] # rr increases, r-inds decrease
        else:
            kw.rbounds = np.sort(kw.rbounds) # rr decreases
            irbounds = inds_from_vals(rr/rsun, kw.rbounds)
    else:
        irbounds = np.sort(kw.irbounds)[::-1] # r-inds derease

    # update the number of quadrants
    nquadr = len(irbounds) - 1

    # update the actual boundary vals
    rbounds = rr[irbounds]/rsun

    # compute the volumes of each quadrant
    volumes = np.zeros(nquadr)
    for ir in range(nquadr): # remember: rbounds increase but r-inds decrease
        ir1 = irbounds[ir + 1]
        ir2 = irbounds[ir]
        volumes[ir] = 4.*np.pi/3.*(rr[ir2]**3. - rr[ir1]**3.)

    # Get the Rayleigh data directory
    radatadir = dirname + '/' + dataname + '/'

    # Get all the file names in datadir and their integer counterparts
    file_list, int_file_list, nfiles = get_file_lists(radatadir, args)

    # Get the problem size
    nproc_min, nproc_max, n_per_proc_min, n_per_proc_max =\
            opt_workload(nfiles, nproc)

    # Distribute file_list to each process
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

        # send  my_files, my_nfiles if nproc > 1
        if k >= 1:
            comm.send([my_files, my_nfiles], dest=k)
else: # recieve my_files, my_nfiles
    my_files, my_nfiles = comm.recv(source=0)

# Broadcast dirname, radatadir, etc.
if rank == 0:
    meta = [\
dirname, dataname, radatadir, irbounds, rw]
else:
    meta = None
dirname, dataname, radatadir, irbounds, rw = comm.bcast(meta, root=0)

# figure out which reading_func to use
reading_func = Shell_Avgs

# Checkpoint and time
comm.Barrier()
if rank == 0:
    t2 = time.time()
    print (format_time(t2 - t1))
    print ("tracing over %i shellular quadrants" %nquadr)
    print ("rbounds/rsun = " + arr_to_str(rbounds, "%.3f"))
    print ('Considering %i %s files for the trace: %s through %s'\
        %(nfiles, dataname, file_list[0], file_list[-1]))
    print(fill_str('computing', lent, char), end='\r')
    t1 = time.time()

# Now analyze the data (some processes may not have nquadr)
nquadr = len(irbounds) - 1
my_times = []
my_iters = []
my_vals = []

for i in range(my_nfiles):
    a = reading_func(radatadir + str(my_files[i]).zfill(8), '')
    for j in range(a.niter):
        if dataname == 'Shell_Avgs':
            vals_loc = a.vals[:, 0, :, :]
        else:
            vals_loc = a.vals

        # get look up table, etc.
        lut = a.lut
        #nq = a.nq
        nq = 7

        # Get values in the separate quadrants, plus Poynting flux
        vals_gav = np.zeros((nq, nquadr))
        # note: the default is nquadr = 1, 1 (yes "extra" dimensions)
        # do volume-avg'd quantities first
        for ir in range(nquadr): # remember: rbounds increase but r-inds decrease
            ir1 = irbounds[ir + 1]
            ir2 = irbounds[ir]

            # radial integration weights
            rw_quad = (rw[ir1:ir2+1]/np.sum(rw[ir1:ir2+1]))

            count = 1
            for qval in [1102,1103,1104,1916,1436]:
                vals_quad = vals_loc[ir1:ir2+1, lut[qval], j]
                vals_gav[count,ir] = np.sum(rw_quad*vals_quad)
                count += 1

        # get total energy
        vals_gav[0,:] = vals_gav[1,:] + vals_gav[2,:] + vals_gav[3,:]

        # get Poynting flux
        for ir in range(nquadr+1): # remember: rbounds increase but r-inds decrease
            vals_gav[6,ir] = vals_loc[irbounds[ir], lut[2001],j]
            # not really "gav" for this guy, but anyway...

        # now we can append the array to the trace
        my_vals.append(vals_gav)
        my_times.append(a.time[j])
        my_iters.append(a.iters[j])
    if rank == 0:
        pcnt_done = i/my_nfiles*100.
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
    # convert everything to arrays
    vals = np.array(vals)
    times = np.array(times)
    iters = np.array(iters)
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
    basename = 'poynt_trace_nquadr%i' %nquadr
    savename = basename + tag + '-' + file_list[0] + '_' + file_list[-1] + '.pkl'
    savefile = datadir + savename

    # save the data
    # Get first and last iters of files
    iter1, iter2 = int_file_list[0], int_file_list[-1]
    f = open(savefile, 'wb')
    di_sav = {'vals': vals, 'times': times, 'iters': iters, 'volumes': volumes, 'volume_full': np.sum(volumes), 'rbounds': rbounds}
    pickle.dump(di_sav, f, protocol=4)
    f.close()
    t2 = time.time()
    print (format_time(t2 - t1))
    print(make_bold(fill_str('total time', lent, char)), end='')
    print (make_bold(format_time(t2 - t1_glob)))
    print ('data saved at ')
    print (make_bold(savefile))
