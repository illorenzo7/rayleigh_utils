# Author: Loren Matilsky
# Created: 01/10/2020
# Reads in kappa00 file byte by byte, using numpy
import numpy as np

def read_inner_vp(the_file):
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
#    nr = int(np.fromfile(f, endian_char + 'i4', 1))
#    kappa00 = np.fromfile(f, endian_char + 'f8', nr)
    inner_vp = np.fromfile(f, endian_char + 'f8') # default is to read
        # all items
    f.close()
    return inner_vp
