# Author: Loren Matilsky
# Created: 01/10/2020
# Reads in kappa00 file byte by byte, using numpy
import numpy as np

def read_eq_vp(the_file, nt, nr):
    eq_vp = np.zeros((nt, nr))
    f = open(the_file, 'rb')
    endian_char = '>' # assume big-endian a priori
    pi_int = int(np.fromfile(f, endian_char + 'i4', 1))
    if pi_int != 314:
        f.close
        f = open(the_file, 'rb')
        pi_int = np.fromfile(f, '<i4', 1)
        endian_char = '<'
        if pi_int != 314:
            sys.exit("In read_inner_vp: Python is confused by the endianness of your data in " + the_file)
    for ir in range(nr):
        for it in range(nt):
            eq_vp[it, ir] = np.fromfile(f, endian_char + 'f8', count=1) 
    f.close()
    return eq_vp
