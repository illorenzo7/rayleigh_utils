##################################################################
# Routine to make temporal spectra
# Author: Loren Matilsky
# Created: 06/12/2021
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
#from nfft import nfft
# data type and reading function
import sys, os
from rayleigh_diagnostics_alt import Shell_Slices
reading_func = Shell_Slices
dataname = 'Shell_Slices'

# my own version of nfft (should be same --- good test)
def my_nfft(times, arr):
    # shift the times to lie in range -1/2, 1/2
    total_time = times[-1] - times[0]
    times = (times - times[0])/total_time - 1/2
    # get equally spaced times
    times_eq = np.linspace(times[0], times[-1], len(times))
    arr_eq = np.interp(times_eq, times, arr)
    return np.fft.fftshift(np.fft.fft(arr_eq))

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
    args = sys.argv
    clas0, clas = read_clas(args)
    dirname = clas0['dirname']

    # Get the Rayleigh data directory
    radatadir = dirname + '/' + dataname + '/'

    # Get all the file names in datadir and their integer counterparts
    file_list, int_file_list, nfiles = get_file_lists(radatadir, args)

    # read first file for some metadata
    a0 = reading_func(radatadir + file_list[0], '')

    # set default values for qval and irval
    kwargs_default = dict({'irvals': np.array([0]), 'rvals': None, 'qvals': np.array([1]), 'nonlin': False, 'mmax': None})

    # overwrite defaults
    kw = update_dict(kwargs_default, clas)

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

    # and the qvals
    qvals = kw.qvals
    if isall(qvals):
        qvals = np.sort(a0.qv)

    if qvals is None:
        qvals = np.array([1])
    
    # everything must be array
    irvals = make_array(irvals)
    qvals = make_array(qvals)

    # use nonlin fft or not
    nonlin = kw.nonlin

    # Get the problem size
    nproc_min, nproc_max, n_per_proc_min, n_per_proc_max =\
            opt_workload(nfiles, nproc)

    # mmax (if needed)
    mmax = kw.mmax

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
    meta = [dirname, radatadir, irvals, qvals, mmax]
else:
    meta = None
dirname, radatadir, irvals, qvals, mmax = comm.bcast(meta, root=0)

# Checkpoint and time
comm.Barrier()
if rank == 0:
    t2 = time.time()
    print (format_time(t2 - t1))
    print (buff_line)
    print ("making %i data file(s)" %(len(irvals)*len(qvals)))
    print ("irvals = ", irvals)
    print ("qvals = ", qvals)
    print (buff_line)
    print ('Considering %i %s files for the trace: %s through %s'\
        %(nfiles, dataname, file_list[0], file_list[-1]))

count = 1
if rank == 0:
    totpklfiles = len(irvals)*len(qvals)

for irval in irvals:
    for qval in qvals:
        if rank == 0:
            print (buff_line)
            fracpklfiles = count/totpklfiles*100
            print ("data file no. %04i of %04i (%5.1f%% done)" %(count, totpklfiles, fracpklfiles))
            print ("irval =", irval)
            print ("qval =", qval)
            if not mmax is None:
                print ("mmax =", mmax)
            print (buff_line)
            print(fill_str('reading data', lent, char), end='\r')
            t1_thisfile = time.time()
            t1 = np.copy(t1_thisfile)

        # Now read the data
        my_times = []
        my_iters = []
        my_vals = []

        for i in range(my_nfiles):
            a = reading_func(radatadir + str(my_files[i]).zfill(8), '', irvals=irval, qvals=qval)
            for j in range(a.niter):
                my_vals.append(a.vals[:, :, 0, 0, j])
                my_times.append(a.time[j])
                my_iters.append(a.iters[j])

            if rank == 0:
                if i == 0:
                    t1_loc = time.time()
                t2_loc = time.time()
                pcnt_done = (i + 1)/my_nfiles*100.0
                print(fill_str('reading data', lent, char) +\
                        ('rank 0 %5.1f%% done' %pcnt_done) + ' ' +\
                        format_time(t2_loc - t1_loc) + 3*' ', end='\r')

        #print ('mvasl =', my_vals)
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
            vals = []

            # Gather the results into these "master" arrays
            for j in range(nproc):
                if j >= 1:
                    # Get data from rank j
                    my_times, my_iters, my_vals = comm.recv(source=j)
                times += my_times
                iters += my_iters
                vals += my_vals

                pcnt_done = (j + 1)/nproc*100
                if j == 0:
                    t1_loc = time.time()
                t2_loc = time.time()
                print(fill_str('rank 0 collecting data from rank %i' %j, lent, char) +\
                        ('%5.1f%% done' %pcnt_done) + ' ' +\
                        format_time(t2_loc - t1_loc) + 3*' ', end='\r')

            vals = np.array(vals) # needs to be array
            # get shape here
            ntimes, nphi, nt = np.shape(vals)
            nphi = 2*nt
        else: # other processes send their data
            comm.send([my_times, my_iters, my_vals], dest=0)

        # Checkpoint and time
        comm.Barrier()
        if rank == 0:
            t2 = time.time()

            print('\n' + fill_str('collection time', lent, char), end='')
            print (format_time(t2 - t1))

            # get the frequencies and mvals
            delta_t = np.mean(np.diff(times))
            freq = np.fft.fftfreq(len(times), delta_t)
            freq = np.fft.fftshift(freq)

            # dealiased no. mvalues
            # remember to dealias---otherwise we are saving a bunch of
            if mmax is None: # may manually strip even more m-values to
                # save space
                nm = int(np.floor(2./3.*nt))
            else:
                nm = mmax
            mvals = np.arange(nm)

            # now redistribute the data by theta-val
            print(fill_str('proc 0 sending datacube by theta-val',\
                    lent, char), end='')
            t1 = time.time()

            nproc_min, nproc_max, n_per_proc_min, n_per_proc_max =\
                    opt_workload(nt, nproc)

            # Distribute itvals to each process
            for k in range(nproc - 1, -1, -1): # go in reverse so my_vals
                # for rank 0 is correct
                # distribute the partial file list to other procs 
                if k < nproc_max: # first processes analyzes more files
                    my_nt = np.copy(n_per_proc_max)
                    itstart = k*my_nt
                    itend = itstart + my_nt
                else: # last processes analyze fewer files
                    my_nt = np.copy(n_per_proc_min)
                    itstart = nproc_max*n_per_proc_max + (k - nproc_max)*my_nt
                    itend = itstart + my_nt

                # Get the portion of vals for rank k
                my_vals = np.copy(vals[:, :, itstart:itend])

                # send appropriate file info if not rank 0
                if k > 0:
                    comm.send([my_vals, my_nt, nonlin, nm, times], dest=k)
        else: # recieve appropriate file info if rank > 1
            my_vals, my_nt, nonlin, nm, times = comm.recv(source=0)

        # make sure everyone gets "their slice"
        comm.Barrier()
        if rank == 0:
            t2 = time.time()
            print (format_time(t2 - t1))

            # everyone do their FFT (or nonlin fft)
            print(fill_str('doing FFT', lent, char), end='\r')
            if nonlin:
                print ("using DFT for NONLINEARLY SPACED times")
            t1 = time.time()

        if rank == 0:
            t1_loc = time.time()
            print(fill_str('doing FFT', lent, char), end='\r')

        # do the phi FFT first
        # this is just a transform of real data
        my_vals = np.fft.rfft(my_vals, axis=1)
        # throw away values above dealised max
        my_vals = my_vals[:, :nm, :]

        # now the time FFT 
        if nonlin:
            for it in range(my_nt):
                for im in range(nm):
                    my_vals[:, im, it] = my_nfft(times, my_vals[:, im, it])
        else:
            my_vals = np.fft.fft(my_vals, axis=0)
            my_vals = np.fft.fftshift(my_vals, axes=0)

        # make it an array
        my_vals = np.array(my_vals)

        # make sure everybody does their part and "gets there"!
        comm.Barrier()

        # rank 0 collects Fourier transform
        if rank == 0:
            t2 = time.time()

            print('\n' + fill_str('FFT time', lent, char), end='')
            print (format_time(t2 - t1))
            print(fill_str('rank 0 collecting FFTs', lent, char), end='')
            t1 = time.time()

            # now overwrite vals with the transformed datacube
            vals = np.zeros((ntimes, nm, nt), 'complex')

            # Gather the FFT results back into vals array
            for k in range(nproc):
                if k < nproc_max: # first processes have more theta-values
                    my_nt = np.copy(n_per_proc_max)
                    itstart = k*my_nt
                    itend = itstart + my_nt
                else: # last processes analyze fewer files
                    my_nt = np.copy(n_per_proc_min)
                    itstart = nproc_max*n_per_proc_max + (k - nproc_max)*my_nt
                    itend = itstart + my_nt

                if k > 0:
                    # Get my_vals from rank k
                     my_vals = comm.recv(source=k, tag=0)
                vals[:, :, itstart:itend] = my_vals

        else: # other processes send their data
            comm.send(my_vals, dest=0, tag=0)

        # Make sure all data is collected
        comm.Barrier()
        if rank == 0:
            t2 = time.time()
            print (format_time(t2 - t1))
            print(fill_str('rank 0 saving data', lent, char), end='')
            t1 = time.time()

            # create data directory if it doesn't already exist
            if not mmax is None:
                datadir = clas0['datadir'] + ('tmspec_mmax%03i/' %mmax)
            else:
                datadir = clas0['datadir'] + 'tmspec/'
            if not os.path.isdir(datadir):
                os.makedirs(datadir)

            # Set the timetrace savename
            savename = ('tmspec_qval%04i_irval%02i' %(qval, irval)) +\
                    clas0['tag'] + '-' + file_list[0] + '_' + file_list[-1] + '.pkl'
            savefile = datadir + savename

            # save the data
            f = open(savefile, 'wb')
            # convert everything to arrays
            vals = np.array(vals)
            times = np.array(times)
            iters = np.array(iters)
            pickle.dump({'vals': vals, 'times': times, 'iters': iters, 'freq': freq, 'mvals': mvals}, f, protocol=4)
            f.close()
            if kw.nonlin:
                print ("the following need not be zero:")
            else:
                print ("the following should be zero:")
            print ("std(dt) = ", np.std(np.diff(times)))
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
