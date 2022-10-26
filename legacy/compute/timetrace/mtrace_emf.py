##################################################################
# Routine to make m spectra (from Shell Slices), as a time trace
# Author: Loren Matilsky
# Created: 02/09/2022
##################################################################

# initialize communication
from mpi4py import MPI
import sys, os
sys.path.append(os.environ['raco'])
from cla_util import *
from common import rsun
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
    print(fill_str('proc 0 distributing the file lists', lent, char),\
            end='')
    t1 = time.time()

# proc 0 reads the file lists and distributes them
if rank == 0:
    # read the arguments
    args = sys.argv
    clas0, clas = read_clas(args)
    dirname = clas0['dirname']

    # Get the Rayleigh data directory
    radatadir = dirname + '/' + dataname + '/'

    # Get all the file names in datadir and their integer counterparts
    file_list, int_file_list, nfiles = get_file_lists(radatadir, args)

    # read first file for some metadata
    a0 = reading_func(radatadir + file_list[0], '')

    # set default values for irvals
    kwargs_default = dict({'irvals': np.array([0]), 'rvals': None,'mmax': None})

    # overwrite defaults
    kw = update_dict(kwargs_default, clas)

    # the grid
    gi = get_grid_info(dirname)
    nt = gi['nt']
    nphi = gi['nphi']
    
    # get the mvals
    # dealiased no. mvalues
    # remember to dealias---otherwise we are saving a bunch of
    # NOTHING!
    mmax = kw.mmax
    if mmax is None: # may manually strip even more m-values to
        # save space
        nm = int(np.floor(2./3.*nt))
    else:
        nm = mmax
    mvals = np.arange(nm)

    # get the rvals we want
    irvals = kw.irvals
    if not kw.rvals is None: # irvals haven't been set directly
        if isall(kw.rvals):
            irvals = np.arange(a0.nr)
        else:
            kw.rvals = make_array(kw.rvals)
            irvals = np.zeros_like(kw.rvals, dtype='int')
            for i in range(len(kw.rvals)):
                irvals[i] = np.argmin(np.abs(a0.radius/rsun - kw.rvals[i]))

   
    # everything must be array
    irvals = make_array(irvals)

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

# Broadcast meta data
if rank == 0:
    meta = [dirname, radatadir, irvals,  mmax, nm, nt]
else:
    meta = None
dirname, radatadir, irvals,  mmax, nm, nt = comm.bcast(meta, root=0)

# Checkpoint and time
comm.Barrier()
if rank == 0:
    t2 = time.time()
    print (format_time(t2 - t1))
    print (buff_line)
    print ("making %i data file(s)" %(len(irvals)*9))
    print ("irvals = ", irvals)
    print (buff_line)
    print ('Considering %i %s files for the trace: %s through %s'\
        %(nfiles, dataname, file_list[0], file_list[-1]))

count = 1
if rank == 0:
    totpklfiles = len(irvals)

for irval in irvals:
    if rank == 0:
        print (buff_line)
        fracpklfiles = count/totpklfiles*100
        print ("data file no. %04i of %04i (%5.1f%% done)" %(count, totpklfiles, fracpklfiles))
        print ("irval = ", irval)
        print (buff_line)
        print(fill_str('reading data', lent, char), end='\r')
        t1_thisfile = time.time()
        t1 = np.copy(t1_thisfile)

    # Now read the data
    my_times = []
    my_iters = []
    my_vals_r = []
    my_vals_r_nodr = []
    my_vals_t = []
    my_vals_t_nodr = []
    my_vals_p = []

    for i in range(my_nfiles):
        avr = reading_func(radatadir + str(my_files[i]).zfill(8), '', irvals=irval, qvals=1)
        avt = reading_func(radatadir + str(my_files[i]).zfill(8), '', irvals=irval, qvals=2)
        avp = reading_func(radatadir + str(my_files[i]).zfill(8), '', irvals=irval, qvals=3)

        abr = reading_func(radatadir + str(my_files[i]).zfill(8), '', irvals=irval, qvals=801)
        abt = reading_func(radatadir + str(my_files[i]).zfill(8), '', irvals=irval, qvals=802)
        abp = reading_func(radatadir + str(my_files[i]).zfill(8), '', irvals=irval, qvals=803)

        for j in range(avr.niter):
            vr = avr.vals[:, :, 0, 0, j]
            vt = avt.vals[:, :, 0, 0, j]
            vp = avp.vals[:, :, 0, 0, j]
            vp_nodr = vp - np.mean(vp, axis=0).reshape((1, nt))

            br = abr.vals[:, :, 0, 0, j]
            bt = abt.vals[:, :, 0, 0, j]
            bp = abp.vals[:, :, 0, 0, j]
           
            emfs_r_phys = [vt*bp - vp*bt, vt*bp, -vp*bt]
            emfs_r_phys_nodr = [vt*bp - vp_nodr*bt, vt*bp, -vp_nodr*bt]
            emfs_t_phys = [vp*br - vr*bp, vp*br, -vr*bp]
            emfs_t_phys_nodr = [vp_nodr*br - vr*bp, vp_nodr*br, -vr*bp]
            emfs_p_phys = [vr*bt - vt*br, vr*bt, -vt*br]

            emfs_r = []
            emfs_r_nodr = []
            emfs_t = []
            emfs_t_nodr = []
            emfs_p = []

            # throw away values above dealised max (nm)
            for iemf in range(3):
                emfs_r.append(np.fft.rfft(emfs_r_phys[iemf],\
                        axis=0)[:nm, :])
                emfs_r_nodr.append(np.fft.rfft(emfs_r_phys_nodr[iemf],\
                        axis=0)[:nm, :])
                emfs_t.append(np.fft.rfft(emfs_t_phys[iemf],\
                        axis=0)[:nm, :])
                emfs_t_nodr.append(np.fft.rfft(emfs_t_phys_nodr[iemf],\
                        axis=0)[:nm, :])
                emfs_p.append(np.fft.rfft(emfs_p_phys[iemf],\
                        axis=0)[:nm, :])

            my_vals_r.append(np.array(emfs_r))
            my_vals_r_nodr.append(np.array(emfs_r_nodr))
            my_vals_t.append(np.array(emfs_t))
            my_vals_t_nodr.append(np.array(emfs_t_nodr))
            my_vals_p.append(np.array(emfs_p))

            my_times.append(avr.time[j])
            my_iters.append(avr.iters[j])

        if rank == 0:
            if i == 0:
                t1_loc = time.time()
            t2_loc = time.time()
            pcnt_done = (i + 1)/my_nfiles*100.0
            print(fill_str('reading data', lent, char) +\
                    ('rank 0 %5.1f%% done' %pcnt_done) + ' ' +\
                    format_time(t2_loc - t1_loc) + 3*' ', end='\r')

    # make sure everybody does their part and "gets there"!
    # Checkpoint and time
    comm.Barrier()
    if rank == 0:
        t2 = time.time()

        print('\n' + fill_str('reading time', lent, char), end='')
        print (format_time(t2 - t1))
        print(fill_str('rank 0 collecting data', lent, char), end='\r')
        t1 = time.time()

    # proc 0 now collects the slices (trace in time)
    # from each process
    if rank == 0:
        # initialize empty lists for the final arrays
        times = []
        iters = []

        vals_r = []
        vals_r_nodr = []
        vals_t = []
        vals_t_nodr = []
        vals_p = []

        # Gather the results into these "master" arrays
        for j in range(nproc):
            if j >= 1:
                # Get data from rank j
                my_times, my_iters, my_vals_r, my_vals_r_nodr,\
                    my_vals_t, my_vals_t_nodr, my_vals_p =\
                    comm.recv(source=j)
            times += my_times
            iters += my_iters

            vals_r += my_vals_r
            vals_r_nodr += my_vals_r_nodr
            vals_t += my_vals_t
            vals_t_nodr += my_vals_t_nodr
            vals_p += my_vals_p

            pcnt_done = (j + 1)/nproc*100
            if j == 0:
                t1_loc = time.time()
            t2_loc = time.time()
            print(fill_str('rank 0 collecting data from rank %i' %j, lent, char) +\
                    ('%5.1f%% done' %pcnt_done) + ' ' +\
                    format_time(t2_loc - t1_loc) + 3*' ', end='\r')

    else: # other processes send their data
        comm.send([my_times, my_iters, my_vals_r, my_vals_r_nodr,\
                    my_vals_t, my_vals_t_nodr, my_vals_p], dest=0)

    # Checkpoint and time
    comm.Barrier()
    if rank == 0:
        t2 = time.time()

        print('\n' + fill_str('collection time', lent, char), end='')
        print (format_time(t2 - t1))

        print(fill_str('rank 0 saving data', lent, char), end='')
        t1 = time.time()

        # create data directory if it doesn't already exist
        if not mmax is None:
            datadir = clas0['datadir'] + ('mtrace_mmax%03i/' %mmax)
        else:
            datadir = clas0['datadir'] + 'mtrace/'

        if not os.path.isdir(datadir):
            os.makedirs(datadir)

        # Set the timetrace savename
        for emftag in ['emfr', 'emfr_nodr', 'emft', 'emft_nodr', 'emfp']:
            savename = ('mtrace_' + emftag + '_irval%02i' %irval) +\
                clas0['tag'] + '-' + file_list[0] + '_' + file_list[-1] + '.pkl'
            savefile = datadir + savename

            # save the data
            f = open(savefile, 'wb')
            # convert everything to arrays
            if emftag == 'emfr':
                vals = np.array(vals_r)
            elif emftag == 'emfr_nodr':
                vals = np.array(vals_r_nodr)
            elif emftag == 'emft':
                vals = np.array(vals_t)
            elif emftag == 'emft_nodr':
                vals = np.array(vals_r_nodr)
            elif emftag == 'emfp':
                vals = np.array(vals_p)
            times = np.array(times)
            iters = np.array(iters)
            pickle.dump({'vals': vals, 'times': times, 'iters': iters, 'mvals': mvals}, f, protocol=4)
            f.close()
            print ('data saved at ')
            print (make_bold(savefile))

            t2 = time.time()
            print (format_time(t2 - t1))
            print(make_bold(fill_str('time for file no. %02i' %count, lent, char)), end='')
            print (make_bold(format_time(t2 - t1_thisfile)))
            print (buff_line)

    count += 1

    # make sure everyone waits on proc 0 before starting the loop again
    comm.Barrier()

if rank == 0:
    t2 = time.time()
    print(make_bold(fill_str('total time', lent, char)), end='')
    print (make_bold(format_time(t2 - t1_glob)))
    print (buff_line)
