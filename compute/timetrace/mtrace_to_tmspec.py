##################################################################
# Routine to make mtrace to tmspec (transform along time axis)
# Author: Loren Matilsky
# Created: 05/01/2022
##################################################################
# import
import numpy as np
import pickle
import sys, os
sys.path.append(os.environ['rapp'])
sys.path.append(os.environ['raco'])
from common import *
from cla_util import *

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
kwargs_default = dict({'rad': False, 'groupname': 'b', 'mvals': 1, 'nonlin': None, 'files': None})

# overwrite defaults
kw = update_dict(kwargs_default, clas)

# get data directory
if kw.rad:
    basename_mtrace = 'timerad'
else:
    basename_mtrace = 'timelat'
datadir_mtrace = clas0['datadir'] + basename_mtrace + '/'

# tmspec directory (output data directory)
if kw.rad:
    basename = 'tmspecrad'
else:
    basename = 'tmspeclat'
datadir = clas0['datadir'] + basename + '/'

# create data directory if it doesn't already exist
if not os.path.isdir(datadir):
    os.makedirs(datadir)

# determine files to transform
kw.mvals = make_array(kw.mvals)
kw.groupname = make_array(kw.groupname)

if kw.files is None:
    kw.files = []
    for mval in kw.mvals:
        for groupname in kw.groupname:
            dataname = basename_mtrace + ('_mval%03i' %mval) + '_' + groupname + clas0['tag'] 
            kw.files.append(get_widest_range_file(datadir_mtrace, dataname))
elif isall(kw.files):
    allfiles = os.listdir(datadir_mtrace)
    kw.files = []
    for the_file in allfiles:
        if 'mval' in the_file:
            kw.files.append(datadir_mtrace + '/' + the_file)
kw.files = np.sort(kw.files)

print (buff_line)
print ("converting %i data file(s)" %(len(kw.files)))

# loop over files and convert data
for the_file in kw.files:
    # get mtrace data
    print ('reading ' + the_file)
    di = get_dict(the_file)
    iter1, iter2 = get_iters_from_file(the_file)
    vals = di['vals']

    # get the times 
    times = di['times']
    mean_dt = np.mean(np.diff(times))
    std_dt = np.std(np.diff(times))
    tol = 1.0e-6
    if kw.nonlin is None: # by default, determine "kw.nonlin" from dispersion
        # of times
        if std_dt/mean_dt > tol:
            kw.nonlin = True
        else:
            kw.nonlin = False

    # Fourier transform the vals
    print (buff_line)
    print ('doing Fourier transform along time axis')
    if kw.nonlin:
        print ("using DFT for NONLINEARLY SPACED times")
        vals_fft, freq = my_nfft(times, vals)
    else:
        vals_fft = np.fft.fft(vals, axis=0)
        vals_fft = np.fft.fftshift(vals_fft, axes=0)
        # get the frequencies
        freq = np.fft.fftfreq(len(times), mean_dt)
        freq = np.fft.fftshift(freq)

    # Set the tmspec savename
    savefile = the_file.replace(basename_mtrace, basename)
    # save the data
    print (buff_line)
    print ('saving ' + savefile)
    f = open(savefile, 'wb')
    di_sav = dict(di)
    di_sav['vals'] = vals_fft
    di_sav['freq'] = freq
    pickle.dump(di_sav, f, protocol=4)
    f.close()

    print (buff_line)
    if kw.nonlin:
        print ("the following need not be zero:")
    else:
        print ("the following should be zero:")
    print ("std(dt)/mean(dt) = %1.3e" %(std_dt/mean_dt))
    print (buff_line)
    print ('data saved at ')
    print (make_bold(savefile))
