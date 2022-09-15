##################################################################
# Author: Loren Matilsky
# Created: 06/17/2021
##################################################################
# This routine computes the time average of the Ferraro magnetic tension
# force:
# (1/4pi) div ( T r^2 sin^2(theta) B_P B_P dot grad B_P
# breaking down by:
# nonzero terms in triple product: tot, pmp, ppm, mmm, mpp, ppp
# production of B_[r, theta, phi]
# induct, shear, advec, comp
# 6 x 3 x 4 = 72 terms
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
sys.path.append(os.environ['raco'])
from rayleigh_diagnostics import AZ_Avgs, Meridional_Slices
from common import drad, dth, get_eq
reading_func1 = AZ_Avgs
reading_func2 = Meridional_Slices
dataname1 = 'AZ_Avgs'
dataname2 = 'Meridional_Slices'

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
    # get the name of the run directory + CLAs
    args = sys.argv
    clas0, clas = read_clas(args)
    if 'instant' in clas:
        instantshear = True
    else:
        instantshear = False # the default

    dirname = clas0['dirname']

    # Get the Rayleigh data directory
    radatadir1 = dirname + '/' + dataname1 + '/'
    radatadir2 = dirname + '/' + dataname2 + '/'

    # Get all the file names in datadir and their integer counterparts
    file_list, int_file_list, nfiles = get_file_lists(radatadir1, args)

    # Get the problem size
    nproc_min, nproc_max, n_per_proc_min, n_per_proc_max =\
            opt_workload(nfiles, nproc)

    # get grid information
    di_grid = get_grid_info(dirname)
    nt = di_grid['nt']
    nr = di_grid['nr']
    rr = di_grid['rr']
    tt = di_grid['tt']
    rr_2d = di_grid['rr_2d']
    sint_2d = di_grid['sint_2d']
    cott_2d = di_grid['cott_2d']
    rr_3d = di_grid['rr_3d']
    sint_3d = di_grid['sint_3d']

    # will need rho and nu to get angular velocity gradients
    eq = get_eq(dirname)
    nu = eq.nu.reshape((1, 1, nr))
    rho = eq.rho.reshape((1, 1, nr))

    # get time averaged shear from longest AZ_Avgs file if required
    if instantshear:
        print ("Will get instantaneaous shear from the AZ_Avgs at each timestep")
        dOmdr_ms = 0.0 # just placeholders...won't be used
        dOmdt_ms = 0.0
    else:
        the_file = get_widest_range_file(clas0['datadir'], 'AZ_Avgs')
        print ('Getting mean shear from ' + the_file)
        di = get_dict(the_file)
        amom_vflux_r = di['vals'][:, :, di['lut'][1813]].reshape((1, nt, nr))
        amom_vflux_t = di['vals'][:, :, di['lut'][1814]].reshape((1, nt, nr))

        # this is the "mean shear" grad Omega
        dOmdr_ms = -amom_vflux_r/(rho*nu*rr_3d**2*sint_3d**2)
        dOmdt_ms = -amom_vflux_t/(rho*nu*rr_3d**2*sint_3d**2)

    # Distribute file lists to each process
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
    meta =[\
dirname, radatadir1, radatadir2, nt, nr, rr, tt, rr_2d, sint_2d, cott_2d, rr_3d, sint_3d, nfiles, instantshear, rho, nu, dOmdr_ms, dOmdt_ms]
else:
    meta = None
dirname, radatadir1, radatadir2, nt, nr, rr, tt, rr_2d, sint_2d, cott_2d, rr_3d, sint_3d, nfiles, instantshear, rho, nu, dOmdr_ms, dOmdt_ms = comm.bcast(meta, root=0)

# Checkpoint and time
comm.Barrier()
if rank == 0:
    t2 = time.time()
    print (format_time(t2 - t1))
    print ('Considering %i %s/%s files for the trace: %s through %s'\
        %(nfiles, dataname1, dataname2, file_list[0], file_list[-1]))
    print(fill_str('computing', lent, char), end='\r')
    t1 = time.time()

# Now analyze the data
nset = 5 # total in a "set" of values
nq = 3*nset # Reynolds dcomp.
my_vals = np.zeros((nt, nr, nq))

my_nfiles = len(my_files)
for i in range(my_nfiles):
    a = reading_func1(radatadir1 + str(my_files[i]).zfill(8), '')
    mer = reading_func2(radatadir2 + str(my_files[i]).zfill(8), '')

    # take mean along the time axis;
    niter = min(a.niter, mer.niter)
    my_weight = 1.0/(nfiles*niter)
    for j in range(niter -1, -1, -1): # go last to first in case "niters" 
        # don't agree
        # full B_pol
        br = mer.vals[:, :, :, mer.lut[801], j]
        bt = mer.vals[:, :, :, mer.lut[802], j]

        # mean B_pol
        br_m = a.vals[:, :, a.lut[801], j].reshape((1, nt, nr))
        bt_m = a.vals[:, :, a.lut[802], j].reshape((1, nt, nr))

        # mean Omega-gradients (instantaneous)
        amom_vflux_r = a.vals[:, :, a.lut[1813], j].reshape((1, nt, nr))
        amom_vflux_t = a.vals[:, :, a.lut[1814], j].reshape((1, nt, nr))

        if instantshear:
            # need rho and nu to get angular velocity gradients
            dOmdr = -amom_vflux_r/(rho*nu*rr_3d**2*sint_3d**2)
            dOmdt = -amom_vflux_t/(rho*nu*rr_3d**2*sint_3d**2)
        else:
            dOmdr = dOmdr_ms
            dOmdt = dOmdt_ms

        # calculate B_phi terms from mean shear
        Bp_fromshear = (br*dOmdr + bt*dOmdt)  # full mean shear
        Bp_fromshear_m = (br_m*dOmdr + bt_m*dOmdt) 
        # mean shear from mean fields
        # still needs to be multiplied by a time scale...

        # next get the angular momentum fluxes from magnetic tension
        factor = rr_3d*sint_3d/(4*np.pi) 
        amom_magflux_r = -np.mean(factor*Bp_fromshear*br, axis=0)
        amom_magflux_t = -np.mean(factor*Bp_fromshear*bt, axis=0)

        amom_magflux_r_m = -np.mean(factor*Bp_fromshear_m*br_m, axis=0)
        amom_magflux_t_m = -np.mean(factor*Bp_fromshear_m*bt_m, axis=0)

        torque_r = - ( drad(amom_magflux_r, rr) + (2/rr_2d)*amom_magflux_r )
        torque_t = - ( (1/rr_2d)*(dth(amom_magflux_t, tt) + cott_2d*amom_magflux_t) )
        torque = torque_r + torque_t

        torque_r_m = - ( drad(amom_magflux_r_m, rr) + (2/rr_2d)*amom_magflux_r_m )
        torque_t_m = - ( (1/rr_2d)*(dth(amom_magflux_t_m, tt) + cott_2d*amom_magflux_t_m) )
        torque_m = torque_r_m + torque_t_m

        # add to my_vals array
        # first the full terms
        indstart = 0
        my_vals[:, :, indstart + 0] += amom_magflux_r*my_weight
        my_vals[:, :, indstart + 1] += amom_magflux_t*my_weight
        my_vals[:, :, indstart + 2] += torque_r*my_weight
        my_vals[:, :, indstart + 3] += torque_t*my_weight
        my_vals[:, :, indstart + 4] += torque*my_weight

        # mean terms
        indstart += nset
        my_vals[:, :, indstart + 0] += amom_magflux_r_m*my_weight
        my_vals[:, :, indstart + 1] += amom_magflux_t_m*my_weight
        my_vals[:, :, indstart + 2] += torque_r_m*my_weight
        my_vals[:, :, indstart + 3] += torque_t_m*my_weight
        my_vals[:, :, indstart + 4] += torque_m*my_weight

    if rank == 0:
        pcnt_done = (i + 1)/my_nfiles*100.
        print(fill_str('computing', lent, char) +\
                ('rank 0 %5.1f%% done' %pcnt_done), end='\r')

# fluc terms
indstart = 2*nset
for k in range(nset):
    my_vals[:, :, indstart + k] +=\
            my_vals[:, :, k] - my_vals[:, :, k + nset]

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
    vals = np.zeros((nt, nr, nq))

    # Gather the results into this "master" array
    for j in range(nproc):
        if j >= 1:
            # Get vals from rank j
            my_vals = comm.recv(source=j)
        # "my_vals" are all weighted: their sum equals the overall average
        vals += my_vals 
else: # other processes send their data
    comm.send(my_vals, dest=0)

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
    if instantshear:
        savename = 'ferraroinst'
    else:
        savename = 'ferraro' # the default
    savename += '-' + file_list[0] + '_' + file_list[-1] + '.pkl'
    savefile = datadir + savename

    # save the data
    # Get first and last iters of files
    f = open(savefile, 'wb')
    pickle.dump({'vals': vals}, f, protocol=4)
    f.close()
    t2 = time.time()
    print (format_time(t2 - t1))
    print(make_bold(fill_str('total time', lent, char)), end='')
    print (make_bold(format_time(t2 - t1_glob)))
    print ('data saved at ')
    print (make_bold(savefile))
