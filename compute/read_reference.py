# Author: Loren Matilsky
# Created: 05/06/2019
# Reads in a custom reference state file byte by byte, using numpy
import numpy as np

def read_reference(ref_file):
    f = open(ref_file, 'rb')
    endian_char = '>' # assume big-endian a priori
    pi_int = int(np.fromfile(f, endian_char + 'i4', 1))
    if pi_int != 314:
        f.close
        f = open(ref_file, 'rb')
        pi_int = np.fromfile(f, '<i4', 1)
        endian_char = '<'
        if pi_int != 314:
            sys.exit("In read_reference: Python is confused by the endianness of your data in " + ref_file)
    nr = int(np.fromfile(f, endian_char + 'i4', 1))
    rr = np.fromfile(f, endian_char + 'f8', nr)
    rho = np.fromfile(f, endian_char + 'f8', nr)
    dlnrho = np.fromfile(f, endian_char + 'f8', nr)
    d2lnrho = np.fromfile(f, endian_char + 'f8', nr)
    p = np.fromfile(f, endian_char + 'f8', nr)
    T = np.fromfile(f, endian_char + 'f8', nr)
    dlnT = np.fromfile(f, endian_char + 'f8', nr)
    dsdr = np.fromfile(f, endian_char + 'f8', nr)
    s = np.fromfile(f, endian_char + 'f8', nr)
    g = np.fromfile(f, endian_char + 'f8', nr)
    f.close()
    return nr, rr, rho, dlnrho, d2lnrho, p, T, dlnT, dsdr, s, g
