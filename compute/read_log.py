# Author: Loren Matilsky
# Created: 02/24/2020
# Reads in the info in a Rayleigh log file (.txt) line by line, using numpy
# Outputs numpy arrays of "iters", "delta_t", and "iters_per_sec" (if 
# available), storing them in a dictionary
import numpy as np

def read_log(fname):
    f = open(fname, 'r')
    lines = f.readlines()
    iters = []
    delta_t = []
    iters_per_sec = []
    for line in lines:
        if "Iteration:" in line:
            split = line.split()
            for i in range(len(split)):
                if split[i] == 'Iteration:':
                    iters.append(int(split[i+1]))
                elif split[i] == 'DeltaT:':
                    delta_t.append(float(split[i+1]))
                elif split[i] == 'Iter/sec:':
                    iters_per_sec.append(float(split[i+1]))
    if len(iters_per_sec) > 0:
        di = dict({'iters': np.array(iters), 'delta_t': np.array(delta_t),\
                'iters_per_sec': np.array(iters_per_sec)})
    else:
        di = dict({'iters': np.array(iters), 'delta_t': np.array(delta_t)})

    return di
