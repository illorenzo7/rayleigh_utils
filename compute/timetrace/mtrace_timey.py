##################################################################
# Routine to trace AZ_Avgs in time/latitude space (pick different radii)
# Author: Loren Matilsky
# Created: 02/27/2019
# Parallelized: 12/12/2020
##################################################################
# This routine computes the trace in time/latitude (by default) or
# time/radius of quantities in the 
# AZ_Avgs data for a particular simulation. 
##################################################################
#
# initialize communication
from mpi4py import MPI
import sys, os
sys.path.append(os.environ['raco'])

comm = MPI.COMM_WORLD
rank = comm.Get_rank()

# Start timing immediately
comm.Barrier()
if rank == 0:
    # timing module
    import time
    # info for print messages
    sys.path.append(os.environ['raco'])
    # import common here
    from common import *
    #from varprops import get_quantity_group
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
    print(fill_str('processes importing necessary modules'), end='')

# modules needed by everyone
import numpy as np
# data type and reading function
import sys, os
sys.path.append(os.environ['rapp'])
from rayleigh_diagnostics_alt import Shell_Slices
reading_func = Shell_Slices
dataname = 'Shell_Slices'

if rank == 0:
    # modules needed only by proc 0 
    import pickle

# Checkpoint and time
comm.Barrier()
if rank == 0:
    t2 = time.time()
    print (format_time(t2 - t1))
    print(fill_str('proc 0 distributing the file lists'), end='')
    t1 = time.time()

# proc 0 reads the file lists and distributes them, also the meta data
if rank == 0:
    # get CLAS
    args = sys.argv
    clas0, clas = read_clas(args)
    dirname = clas0['dirname']
    magnetism = clas0['magnetism']
    kwargs_default = dict({'rad': False, 'samplevals': None, 'rvals': None, 'qvals': None, 'groupname': 'v', 'mvals': np.array([1])})
    kw = update_dict(kwargs_default, clas)

    # can control samplevals with rvals directly:
    if not kw.rad and not kw.rvals is None:
        kw.samplevals = kw.rvals

    # get default sampling locations (samplevals)
    slice_info = get_sliceinfo(dirname)
    di_grid = get_grid_info(dirname)
    rr = slice_info.samplevals
    tt_lat = di_grid['tt_lat']

    if kw.rad:
        sampleaxis = tt_lat
    else:
        sampleaxis = rr

    if kw.samplevals is None:
        if kw.rad:
            samplevals = np.array([-85.] + np.linspace(-75., 75., 11).tolist() + [85.])
        else:
            samplevals = rr
        sampletag = ''
    else:
        samplevals = kw.samplevals
        sampletag = '_' + input("choose a tag name for your chosen sampling locations: ")

    # get sampling quantities we want
    if kw.qvals is None: # it's a quantity group
        groupname = kw.groupname
        qgroup = get_quantity_group(groupname, magnetism)
        qvals = qgroup['qvals']
    else:
        qvals = kw.qvals
        groupname = input("choose a groupname to save your data: ")

    # get mvals we want (default m = 1)
    mvals = make_array(kw.mvals)

    # get the Rayleigh data directory
    radatadir = dirname + '/' + dataname + '/'

    # get desired analysis range
    file_list, int_file_list, nfiles = get_file_lists(radatadir, clas)

    # get the problem size
    nproc_min, nproc_max, n_per_proc_min, n_per_proc_max =\
            opt_workload(nfiles, nproc)

    isamplevals = []
    for sampleval in samplevals:
        isamplevals.append(np.argmin(np.abs(sampleaxis - sampleval)))
    isamplevals = np.array(isamplevals)
    # recompute the actual sample values we get
    samplevals = sampleaxis[isamplevals]
    nsamplevals = len(samplevals)

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

        # send appropriate file info if nproc > 1
        if k >= 1:
            comm.send([my_files, my_nfiles], dest=k)
else: # recieve appropriate file info if rank > 1
    my_files, my_nfiles = comm.recv(source=0)

# broadcast meta data
if rank == 0:
    meta = [dirname, radatadir, qvals, kw.rad, isamplevals, nsamplevals, mvals]
else:
    meta = None

the_bcast = comm.bcast(meta, root=0)
dirname, radatadir, qvals, rad, isamplevals, nsamplevals, mvals = the_bcast

# Checkpoint and time
comm.Barrier()
if rank == 0:
    t2 = time.time()
    print (format_time(t2 - t1))
    print ('Considering %i %s files for the trace: %s through %s'\
        %(nfiles, dataname, file_list[0], file_list[-1]))
    print ("sampling values:")
    print ("qvals = " + arr_to_str(qvals, "%i"))
    print ("mvals = " + arr_to_str(mvals, "%i"))
    st = "sampling locations:"
    if rad:
        st2 = "lats = "
        fmt = '%.1f'
    else:
        st2 = "rvals = "
        fmt = "%1.3e"
    print(st)
    print (st2 + arr_to_str(samplevals, fmt))
    print(fill_str('computing'), end='\r')
    t1 = time.time()

# Now analyze the data
my_times = []
my_iters = []
my_vals = []

# no point in storing extraneous mvals
mmax = max(mvals)
for i in range(my_nfiles):
    a = reading_func(radatadir + str(my_files[i]).zfill(8), '', qvals=qvals)
    for j in range(a.niter):
        vals_loc = a.vals[..., j]
        vals_loc = np.fft.rfft(vals_loc, axis=0)
        vals_loc = vals_loc[:mmax+1, ...]

        if rad:
            my_vals.append(vals_loc[:, isamplevals, :,  :][:, :, :, a.lut[qvals]])
        else:
            my_vals.append(vals_loc[:, :, isamplevals, :][:, :, :, a.lut[qvals]])

        my_times.append(a.time[j])
        my_iters.append(a.iters[j])
    if rank == 0:
        pcnt_done = i/my_nfiles*100.0
        print(fill_str('computing') +\
                ('rank 0 %5.1f%% done' %pcnt_done), end='\r')

# Checkpoint and time
comm.Barrier()
if rank == 0:
    t2 = time.time()
    print('\n' + fill_str('computing time'), end='')
    print (format_time(t2 - t1))
    print(fill_str('rank 0 collecting and saving the results'), end='')
    t1 = time.time()

# proc 0 now collects the results from each process
if rank == 0:
    # initialize empty lists for the final arrays
    times = []
    iters = []
    vals = []

    # Gather the results into these "master" arrays
    istart = 0
    for j in range(nproc):
        if j >= 1:
            # Get my_my_times, my_iters, my_vals from rank j
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
    # Set the timetrace savename by the directory, what we are saving,
    # and first and last iteration files for the trace

    # convert everything to arrays
    vals = np.array(vals)
    times = np.array(times)
    iters = np.array(iters)

    # save the data for each mval individually
    for mval in mvals:
        # save the data in a sensible name
        if rad:
            basename = 'timerad'
        else:
            basename = 'timelat'

        # create data directory if it doesn't already exist
        datadir = clas0['datadir'] + basename + '/'
        if not os.path.isdir(datadir):
            os.makedirs(datadir)

        basename += ('_mval%03i' %mval) + '_' + groupname + sampletag

        savename = basename + clas0['tag'] + '-' +\
                file_list[0] + '_' + file_list[-1] + '.pkl'
        savefile = datadir + savename

        # reduce the vals array
        vals_sav = vals[:, mval, :, :]

        # save the data
        f = open(savefile, 'wb')
        di_sav = dict({'vals': vals_sav, 'times': times, 'iters': iters, 'qvals': qvals})
        di_sav['samplevals'] = samplevals
        pickle.dump(di_sav, f, protocol=4)
        f.close()
        print ('data saved at ')
        print (make_bold(savefile))

    t2 = time.time()
    print (format_time(t2 - t1))
    print (make_bold(fill_str('total time')), end='')
    print (make_bold(format_time(t2 - t1_glob)))
