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
import time
from mpi4py import MPI
import sys, os
sys.path.append(os.environ['raco'])
comm = MPI.COMM_WORLD
rank = comm.Get_rank()

# Start timing immediately
comm.Barrier()
if rank == 0:
    # timing module
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
from rayleigh_diagnostics_alt import Shell_Slices
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
    # get CLAS
    args = sys.argv
    clas0, clas = read_clas(args)
    dirname = clas0['dirname']

    # Get the Rayleigh data directory
    radatadir = dirname + '/' + dataname + '/'

    # Get all the file names in datadir and their integer counterparts
    file_list, int_file_list, nfiles = get_file_lists(radatadir, args)

    # read first file for some metadata
    a0 = reading_func(radatadir + file_list[0], '')

    # Set other defaults
    kwargs_default = dict({'clat': 10, 'dlat': 0, 'irvals': np.array([0]), 'rvals': None, 'qvals': np.array([1])})

    # overwrite defaults
    kw = update_dict(kwargs_default, clas)

    # get the rvals we want
    rvals = kw.rvals
    irvals = kw.irvals
    if not kw.rvals is None: # irvals haven't been set directly
        if isall(rvals):
            irvals = np.arange(a0.nr)
        else:
            irvals = np.zeros_like(kw.rvals, dtype='int')
            for i in range(len(kw.rvals)):
                irvals[i] = np.argmin(np.abs(a0.radius/rsun - kw.rvals[i]))

    # and the qvals
    qvals = make_array(kw.qvals)

    # everything must be array
    irvals = make_array(irvals)

    # Get metadata
    rvals = a0.radius[irvals]/rsun
    clat = kw.clat
    dlat = kw.dlat

    # didn't put this earlier because it got messed up by the message
    # saying which depths we were taking
    print(fill_str('proc 0 distributing the file lists', lent, char),\
            end='')
 
    # Desired latitude range + nphi
    di_grid = get_grid_info(dirname)
    ith1 = np.argmin(np.abs(di_grid['tt_lat'] - (clat - dlat/2.)))
    ith2 = np.argmin(np.abs(di_grid['tt_lat'] - (clat + dlat/2.)))
    # recompute the actual clat and dlat
    clat = (di_grid['tt_lat'][ith2] + di_grid['tt_lat'][ith2])/2
    dlat = di_grid['tt_lat'][ith2] - di_grid['tt_lat'][ith2]
    tw_strip = di_grid['tw'][ith1:ith2+1]
    tw_strip /= np.sum(tw_strip)
    tw_strip = tw_strip.reshape((1, ith2 - ith1 + 1))
    nphi =  di_grid['nphi']

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

# Broadcast meta data
if rank == 0:
    meta = [dirname, radatadir, qvals, irvals, nphi,\
            ith1, ith2, tw_strip]
else:
    meta = None
dirname, radatadir, qvals, irvals, nphi,\
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
        print("dlat = %.2f in latitude" %dlat)
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

# loop over rvals and qvals, and make data
count = 1
if rank == 0:
    totpklfiles = len(irvals)*len(qvals)

for irval in irvals:
    for qval in qvals:
        if rank == 0:
            print (buff_line)
            fracpklfiles = count/totpklfiles*100
            print ("data file no. %04i of %04i (%5.1f%% done)" %(count, totpklfiles, fracpklfiles))
            print ("irval = ", irval)
            print ("qval = ", qval)
            print (buff_line)
            print(fill_str('reading data', lent, char), end='\r')
            t1_thisfile = time.time()
            t1 = np.copy(t1_thisfile)

        # Now analyze the data, for each rval, qval
        my_times = []
        my_iters = []
        my_vals = []

        for i in range(my_nfiles):
            a = reading_func(radatadir + str(my_files[i]).zfill(8), '', irvals=irval, qvals=qval)
            for j in range(a.niter):
                # space in the arrays
                # get desired values in the strip
                vals_strip = a.vals[:, ith1:ith2+1, 0, 0, j]
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
            datadir = clas0['datadir'] + 'timelon/'
            if not os.path.isdir(datadir):
                os.makedirs(datadir)

            # Set the timetrace savename
            savename = 'timelon_clat' + lat_format(clat) + '_dlat%03.0f' %dlat + ('_qval%04i_irval%02i' %(qval, irval)) +\
                    clas0['tag'] + '-' + file_list[0] + '_' + file_list[-1] + '.pkl'
            savefile = datadir + savename

            # save the data
            f = open(savefile, 'wb')
            # convert everything to arrays
            vals = np.array(vals)
            times = np.array(times)
            iters = np.array(iters)
            pickle.dump({'vals': vals, 'times': times, 'iters': iters, 'clat': clat, 'dlat': dlat}, f, protocol=4)
            f.close()
            print ('data saved at ')
            print (make_bold(savefile))

            t2 = time.time()
            print (format_time(t2 - t1))
            print(make_bold(fill_str('time for file no. %02i' %count, lent, char)), end='')
            print (make_bold(format_time(t2 - t1_thisfile)))
            print (buff_line)
        count += 1

if rank == 0:
    t2 = time.time()
    print(make_bold(fill_str('total time', lent, char)), end='')
    print (make_bold(format_time(t2 - t1_glob)))
    print (buff_line)
