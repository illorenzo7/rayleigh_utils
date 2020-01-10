# Author: Loren Matilsky
# Created: 01/10/2020
# Reads in kappa00 file byte by byte, using numpy
import numpy as np

def read_kappa00(kappa00_file):
    f = open(kappa00_file, 'rb')
    endian_char = '>' # assume big-endian a priori
    pi_int = int(np.fromfile(f, endian_char + 'i4', 1))
    if pi_int != 314:
        f.close
        f = open(kappa00_file, 'rb')
        pi_int = np.fromfile(f, '<i4', 1)
        endian_char = '<'
        if pi_int != 314:
            sys.exit("In read_kappa00: Python is confused by the endianness of your data in " + kappa00_file)
    nr = int(np.fromfile(f, endian_char + 'i4', 1))
    rr = np.fromfile(f, endian_char + 'f8', nr)
    kappa00 = np.fromfile(f, endian_char + 'f8', nr)
    dlnkappa00 = np.fromfile(f, endian_char + 'f8', nr)
    f.close()
    return nr, rr, kappa00, dlnkappa00
