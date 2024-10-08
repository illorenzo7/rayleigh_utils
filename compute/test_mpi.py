##################################################################
# Routine to test MPI, making sure mpi4py works and all processes
# respond
##################################################################
##################################################################

# initialize communication and import modules
from mpi4py import MPI
import numpy as np
import sys
comm = MPI.COMM_WORLD
rank = comm.Get_rank()

comm.Barrier() # sync processes

if rank == 0:
    print ("=============")
    nproc = comm.Get_size()
    if nproc > 1:
        print ('testing parallel computation with %i ranks' %nproc)
        print ('communication initialized')
    else:
        print ('testing serial computation with 1 rank')

# modules needed by everyone
import numpy as np

comm.Barrier() # sync processes

a = np.random.rand(100)
mean = np.mean(a)
print ("rank %i checking in. I computed the mean of 100 random numbers: %.3f" %(rank, mean))

sys.stdout.flush()
comm.Barrier() # sync processes

# proc 0 now collects the results from each process
if rank == 0:
    # initialize empty lists for the final arrays
    totmean = 0.
    for j in range(nproc):
        if j >= 1:
            mean = comm.recv(source=j)
        totmean += mean
    totmean /= nproc

    print("=============")
    print("rank 0 here: the mean of %i means = %.3f" %(nproc,totmean))

else: # other processes send their data
    comm.send(mean, dest=0)
