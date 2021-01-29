# Routine to trace Rayleigh G_Avgs data in time
# Created by: Loren Matilsky
# On: 01/28/2021
##################################################################
# This routine computes the trace in time of the values in the G_Avgs data 
# for a particular simulation. 
# traces separtely over different "quadrants" (4-8 in total, depending on
# if rbcz is specified to separate the domains in radius)

# By default, the routine traces over the last 100 files of datadir, though
# the user can specify a different range in sevaral ways:
# -n 10 (last 10 files)
# -range iter1 iter2 (names of start and stop data files; iter2 can be 'last')
# -centerrange iter0 nfiles (trace about central file iter0 over nfiles)

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
    print ('%8.2e s' %(t2 - t1))
    print(fill_str('proc 0 distributing the file lists', lent, char),\
            end='')
    t1 = time.time()

# proc 0 reads the file lists and distributes them
if rank == 0:
    # Get the name of the run directory
    dirname = sys.argv[1]
    args = sys.argv[2:]
    nargs = len(args)

    # Get the Rayleigh data directory
    radatadir = dirname + '/' + dataname + '/'

    # Get all the file names in datadir and their integer counterparts
    file_list, int_file_list, nfiles = get_file_lists(radatadir)
    
    # get desired range for analysis
    the_tuple = get_desired_range(int_file_list, args)
    if the_tuple is None:
        index_first, index_last = nfiles - 101, nfiles - 1  
        # By default trace over the last 100 files
    else:
        index_first, index_last = the_tuple

    # Remove parts of file lists we don't need
    file_list = file_list[index_first:index_last + 1]
    int_file_list = int_file_list[index_first:index_last + 1]
    nfiles = index_last - index_first + 1

    # Get the problem size
    nproc_min, nproc_max, n_per_proc_min, n_per_proc_max =\
            opt_workload(nfiles, nproc)

    # Will need the first data file for a number of things
    a0 = reading_func(radatadir + file_list[0], '')
    nrec_full = a0.niter
    nq = a0.nq + 3 # will add the three internal energies after

    # Will need nrec (last niter) to get proper time axis size
    af = reading_func(radatadir + file_list[-1], '')
    nrec_last = af.niter
    ntimes = (nfiles - 1)*nrec_full + nrec_last

    # get grid information
    gi = GridInfo(dirname + '/grid_info', '')
    rr = gi.radius
    cost = gi.costheta
    sint = gi.sintheta
    tt = np.arccos(cost)
    tt_lat = 180./np.pi*(np.pi/2. - tt)
    rw = gi.rweights
    nr = gi.nr
    tw = gi.tweights
    nt = gi.ntheta

    # domain bounds
    ncheby, domain_bounds = get_domain_bounds(dirname)
    ri = np.min(domain_bounds)
    ro = np.max(domain_bounds)
    d = ro - ri

    # get zone separators, defaulting to 
    # the domain bounds in radius and
    # [-45, 0, 45] in latitude

    # Read in CLAs
    latvals = [-45., 0., 45.]
    rvals = []
    for i in range(nargs):
        arg = args[i]
        if arg == '-rbcz':
            rvals = [float(args[i+1])*rsun]
        elif arg == '-rbczcm':
            rvals = [float(args[i+1])]
        elif arg == '-depths':
            strings = args[i+1].split()
            for st in strings:
                rval = ro - float(st)*d
                rvals.append(rval)
        elif arg == '-depthscz':
            rm = domain_bounds[1]
            dcz = ro - rm
            strings = args[i+1].split()
            for st in strings:
                rval = ro - float(st)*dcz
                rvals.append(rval)
        elif arg == '-depthsrz':
            rm = domain_bounds[1]
            drz = rm - ri
            strings = args[i+1].split()
            for st in strings:
                rval = rm - float(st)*drz
                rvals.append(rval)
        elif arg == '-rvals':
            rvals = []
            strings = args[i+1].split()
            for st in strings:
                rval = float(st)*rsun
                rvals.append(rval)
        elif arg == '-rvalscm':
            rvals = []
            strings = args[i+1].split()
            for st in strings:
                rval = float(st)
                rvals.append(rval)
        elif arg == '-rrange':
            r1 = float(args[i+1])
            r2 = float(args[i+2])
            n = int(args[i+3])
            rvals = np.linspace(r1, r2, n)*rsun
        elif arg == '-lats':
            latvals = []
            strings = args[i+1].split()
            for st in strings:
                latval = float(st)
                latvals.append(latval)
            latvals = latvals[::-1] + [0.] + latvals
            # indicators separtors in only one hemisphere

    ndomains = len(ncheby)
    if rvals == [] and ndomains > 1:
        rvals = list(domain_bounds[1:-1][::-1])
    rbounds = [ro] + rvals + [ri]
    latbounds = [tt_lat[0]] + latvals + [tt_lat[-1]]

    # now get the separator indices
    ir_sep = []
    for rval in rvals:
        ir_sep.append(np.argmin(np.abs(rr - rval)))
    it_sep = []
    for lat in latvals:
        it_sep.append(np.argmin(np.abs(tt_lat - lat)))

    nsep_r = len(ir_sep)
    nsep_t = len(it_sep)
    nquad = (nsep_r + 1)*(nsep_t + 1)

    # compute the volumes of each quadrant
    volumes = np.zeros((nsep_t + 1, nsep_r + 1))
    for it in range(nsep_t + 1):
        if it == 0:
            it1 = 0
        else:
            it1 = it_sep[it - 1]
        if it == nsep_t:
            it2 = nt - 1
        else:
            it2 = it_sep[it]
        for ir in range(nsep_r + 1):
            if ir == 0:
                ir1 = 0
            else:
                ir1 = ir_sep[ir - 1]
            if ir == nsep_r:
                ir2 = nr - 1
            else:
                ir2 = ir_sep[ir]
            volumes[it, ir] = 2.*np.pi/3.*(rr[ir1]**3. - rr[ir2]*3.)*\
                    np.abs(cost[it1] - cost[it2])

    # Get rho*T
    eq = get_eq(dirname)
    rhot = (eq.density*eq.temperature).reshape((1, nr))

    # Distribute file_list and my_ntimes to each process
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

        if k == nproc - 1: # last process may have nrec_last != nrec_full
            my_ntimes = (my_nfiles - 1)*nrec_full + nrec_last
        else:
            my_ntimes = my_nfiles*nrec_full

        # Get the file list portion for rank k
        my_files = np.copy(int_file_list[istart:iend])

        # send  my_files, my_nfiles, my_ntimes if nproc > 1
        if k >= 1:
            comm.send([my_files, my_nfiles, my_ntimes], dest=k)
else: # recieve my_files, my_nfiles, my_ntimes
    my_files, my_nfiles, my_ntimes = comm.recv(source=0)

# Broadcast dirname, radatadir, nq
if rank == 0:
    meta = [dirname, radatadir, nq, nt, nr, ir_sep, it_sep, rhot, nsep_r, nsep_t, tw, rw]
else:
    meta = None
dirname, radatadir, nq, nt, nr, ir_sep, it_sep, rhot, nsep_r, nsep_t, tw, rw = comm.bcast(meta, root=0)

# Checkpoint and time
comm.Barrier()
if rank == 0:
    t2 = time.time()
    print ('%8.2e s' %(t2 - t1))
    print ("tracing over %i quadrants" %nquad)
    print (("nsep_t = %i" %nsep_r), " latbounds = ", latbounds)
    print (("nsep_r = %i" %nsep_r), " rbounds = ", rbounds)
    print ('Considering %i %s files for the trace: %s through %s'\
        %(nfiles, dataname, file_list[0], file_list[-1]))
    print ("ntimes for trace = %i" %ntimes)
    print(fill_str('computing', lent, char), end='\r')
    t1 = time.time()

# Now analyze the data
my_times = np.zeros(my_ntimes)
my_iters = np.zeros(my_ntimes, dtype='int')
my_vals = np.zeros((my_ntimes, nq, nsep_t + 1, nsep_r + 1))

my_count = 0
for i in range(my_nfiles):
    if rank == 0 and i == 0:
        a = a0
    else:   
        a = reading_func(radatadir + str(my_files[i]).zfill(8), '')
    for j in range(a.niter):
        if my_count < my_ntimes: # make sure we don't go over the allotted
            # space in the arrays
            vals_loc = np.copy(a.vals[:, :, :, j])
            # add in internal energy
            inte_loc = rhot*vals_loc[:, :, a.lut[501]]
            # top S subtracted
            topS = np.sum(tw*vals_loc[:, 0, a.lut[501]])
            inte_loc_subt = rhot*(vals_loc[:, :, a.lut[501]] - topS)
            # bottom S subtracted
            botS = np.sum(tw*vals_loc[:, -1, a.lut[501]])
            inte_loc_subb = rhot*(vals_loc[:, :, a.lut[501]] - botS)

            # add in the three energies
            vals_loc = np.concatenate((vals_loc, inte_loc.reshape((nt, nr, 1)), inte_loc_subt.reshape((nt, nr, 1)),\
                    inte_loc_subb.reshape((nt, nr, 1))), axis=2)

            # Get the values in the separate quadrants
            for it in range(nsep_t + 1):
                if it == 0:
                    it1 = 0
                else:
                    it1 = it_sep[it - 1]
                if it == nsep_t:
                    it2 = nt - 1
                else:
                    it2 = it_sep[it]
                for ir in range(nsep_r + 1):
                    if ir == 0:
                        ir1 = 0
                    else:
                        ir1 = ir_sep[ir - 1]
                    if ir == nsep_r:
                        ir2 = nr - 1
                    else:
                        ir2 = ir_sep[ir]

                    vals_quad = vals_loc[it1:it2+1, ir1:ir2+1, :]
                    tw_quad = (tw[it1:it2+1]/np.sum(tw[it1:it2+1])).\
                            reshape((it2 - it1 + 1, 1, 1))
                    rw_quad = (rw[ir1:ir2+1]/np.sum(rw[ir1:ir2+1])).\
                            reshape((ir2 - ir1 + 1, 1))


                    gav = np.sum(tw_quad*vals_quad, axis=0)
                    gav = np.sum(rw_quad*gav, axis=0)
                    my_vals[my_count, :, it, ir] = gav
            my_times[my_count] = a.time[j] 
            my_iters[my_count] = a.iters[j]
        my_count += 1
        if rank == 0:
            pcnt_done = my_count/my_ntimes*100.
            print(fill_str('computing', lent, char) +\
                    ('rank 0 %5.1f%% done' %pcnt_done), end='\r')
# Checkpoint and time
comm.Barrier()
if rank == 0:
    t2 = time.time()
    print(fill_str('\ncomputing time', lent, char), end='')
    print ('%8.2e s                                 ' %(t2 - t1))
    print(fill_str('rank 0 collecting and saving the results',\
            lent, char), end='')
    t1 = time.time()

# proc 0 now collects the results from each process
if rank == 0:
    # Initialize zero-filled 'vals/times/iters' arrays to store the data
    times = np.zeros(ntimes)
    iters = np.zeros(ntimes, dtype='int')
    vals = np.zeros((ntimes, nq, nsep_t + 1, nsep_r + 1))

    # Gather the results into these "master" arrays
    istart = 0
    for j in range(nproc):
        if j >= 1:
            # Get my_ntimes, my_times, my_iters, my_vals from rank j
            my_ntimes, my_times, my_iters, my_vals = comm.recv(source=j)
        times[istart:istart+my_ntimes] = my_times
        iters[istart:istart+my_ntimes] = my_iters
        vals[istart:istart+my_ntimes, :, :, :] = my_vals
        istart += my_ntimes
else: # other processes send their data
    comm.send([my_ntimes, my_times, my_iters, my_vals], dest=0)

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
    savename = dirname_stripped + '_trace_quad_G_Avgs_' +\
            file_list[0] + '_' + file_list[-1] + '.pkl'
    savefile = datadir + savename

    # compute the full average over the whole volume
    vals_full = np.zeros((ntimes, nq))
    for it in range(nsep_t + 1):
        for ir in range(nsep_r + 1):
            vals_full += volumes[it, ir]*vals[:, :, it, ir]
    vals_full /= np.sum(volumes)

    # save the data
    # Get first and last iters of files
    iter1, iter2 = int_file_list[0], int_file_list[-1]
    # append the lut (making the inte, inte_subt, and inte_subb quantities
    # (4000, 4001, 4002)
    lut_app = np.array([a0.nq, a0.nq + 1, a0.nq + 2])
    lut = np.hstack((a0.lut, lut_app))
    qv_app = np.array([4000, 4001, 4002])
    qv = np.hstack((a0.qv, qv_app))
    f = open(savefile, 'wb')
    pickle.dump({'vals': vals, 'vals_full': vals_full,'times': times, 'iters': iters, 'lut': lut, 'ntimes': ntimes, 'iter1': iter1, 'iter2': iter2, 'rr': a0.radius, 'nr': nr, 'nt': nt, 'qv': qv, 'nq': nq, 'nsep_r': nsep_r, 'nsep_t': nsep_t, 'it_sep': it_sep, 'ir_sep': ir_sep, 'volumes': volumes, 'latvals_sep': latvals, 'rvals_sep': rvals, 'rbounds': rbounds, 'latbounds': latbounds, 'nquad': nquad}, f, protocol=4)
    f.close()
    t2 = time.time()
    print ('%8.2e s' %(t2 - t1))
    print(fill_str('total time', lent, char), end='')
    print ('%8.2e s' %(t2 - t1_glob))
