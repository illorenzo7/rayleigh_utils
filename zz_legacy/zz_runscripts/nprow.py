import numpy as np
import sys

ncpus = int(sys.argv[1])

root = np.sqrt(ncpus)

# find the next highest power of 2
exp = int(np.ceil(np.log2(root)))

nprow = 2**exp
npcol = ncpus // nprow

print (nprow)
