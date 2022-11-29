##################################################################
# Routine to make mtrace to tmspec (transform along time axis)
# Author: Loren Matilsky
# Created: 05/01/2022
##################################################################

# import
import sys, os
sys.path.append(os.environ['raco'])
from cla_util import *
from common import *

import numpy as np
import pickle

# CLAs
args = sys.argv
clas0, clas = read_clas(args)
dirname = clas0['dirname']

# Get the Rayleigh data directory
radatadir = dirname + '/Shell_Slices/'

# Get all the file names in datadir and their integer counterparts
file_list, int_file_list, nfiles = get_file_lists(radatadir, clas)

# read first file for some metadata
a0 = Shell_Slices(radatadir + file_list[0], '')

# set default values for qval and irval
kwargs_default = dict({'irvals': np.array([0]), 'rvals': None, 'qvals': None, 'mmax': None, 'nonlin': None})

# overwrite defaults
kw = update_dict(kwargs_default, clas)
nonlin = kw.nonlin

# the grid
gi = get_grid_info(dirname)
nt = gi['nt']
nphi = gi['nphi']

# get the mvals
mmax = kw.mmax
if mmax is None: # may manually strip even more m-values to
    # save space
    nm = int(np.floor(2./3.*nt))
else:
    nm = mmax
mvals = np.arange(nm)

# get the rvals we want
irvals = kw.irvals
if not kw.rvals is None: # irvals haven't been set directly
    if isall(kw.rvals):
        irvals = np.arange(a0.nr)
    else:
        kw.rvals = make_array(kw.rvals)
        irvals = np.zeros_like(kw.rvals, dtype='int')
        for i in range(len(kw.rvals)):
            irvals[i] = np.argmin(np.abs(a0.radius/rsun - kw.rvals[i]))

# and the qvals
qvals = kw.qvals
if isall(qvals):
    qvals = np.sort(a0.qv)

if qvals is None:
    qvals = np.array([1])

# everything must be array
irvals = make_array(irvals)
qvals = make_array(qvals)

# get mtrace directory
if not mmax is None:
    datadir_mtrace = clas0['datadir'] + ('mtrace_mmax%03i/' %mmax)
else:
    datadir_mtrace = clas0['datadir'] + 'mtrace/'

# tmspec directory
datadir = datadir_mtrace.replace('mtrace', 'tmspec')

# create data directory if it doesn't already exist
if not os.path.isdir(datadir):
    os.makedirs(datadir)

print (buff_line)
print ("converting %i data file(s)" %(len(irvals)*len(qvals)))
print ("irvals = ", irvals)
print ("qvals = ", qvals)

# loop over rvals and qvals and convert data
for irval in irvals:
    for qval in qvals:

        # names of datafiles
        dataname_mtrace = ('mtrace_qval%04i_irval%02i' %(qval, irval)) +\
                clas0['tag']
        dataname = dataname_mtrace.replace('mtrace', 'tmspec')

        # get mtrace data
        the_file = get_widest_range_file(datadir_mtrace, dataname_mtrace)
        print (buff_line)
        print ('reading ' + the_file)
        di = get_dict(the_file)
        iter1, iter2 = get_iters_from_file(the_file)
        vals = di['vals']

        # get the times 
        times = di['times']
        delta_t = np.mean(np.diff(times))
        if kw.nonlin is None: # by default, determine "nonlin" from dispersion
            tol = 1e-6
            # of times
            disp_t = np.sqrt(np.mean((np.diff(times) - delta_t)**2))
            if disp_t/delta_t > tol:
                nonlin = True
            else:
                nonlin = False

        # Fourier transform the vals
        print (buff_line)
        print ('doing Fourier transform along time axis')
        if nonlin:
            print ("using DFT for NONLINEARLY SPACED times")
            vals_fft, freq = my_nfft(times, vals)
        else:
            vals_fft = np.fft.fft(vals, axis=0)
            vals_fft = np.fft.fftshift(vals_fft, axes=0)
            # get the frequencies
            delta_t = np.mean(np.diff(times))
            freq = np.fft.fftfreq(len(times), delta_t)
            freq = np.fft.fftshift(freq)

        # Set the tmspec savename
        savename = ('tmspec_qval%04i_irval%02i' %(qval, irval)) +\
                clas0['tag'] + '-' + str(iter1).zfill(8) + '_' +\
                str(iter2).zfill(8) + '.pkl'
        savefile = datadir + savename

        # save the data
        print (buff_line)
        print ('saving ' + savefile)
        f = open(savefile, 'wb')

        pickle.dump({'vals': vals_fft, 'times': times, 'iters': di['iters'], 'freq': freq, 'mvals': mvals}, f, protocol=4)

        f.close()
        print (buff_line)
        if nonlin:
            print ("the following need not be zero:")
        else:
            print ("the following should be zero:")
        print ("std(dt) = ", np.std(np.diff(times)))
        print (buff_line)
        print ('data saved at ')
        print (make_bold(savefile))

        print (buff_line)
        print (buff_line)
