# Author: Loren Matilsky
# Created: 11/10/2018
# Common routines for routines for post-processing Rayleigh data 

import numpy as np
from scipy.interpolate import interp1d
import sys, os, pickle
from string_to_num import string_to_number_or_array
sys.path.append(os.environ['rapp'])
from reference_tools import equation_coefficients
from rayleigh_diagnostics import G_Avgs, Shell_Slices, ReferenceState,\
    TransportCoeffs, GridInfo
from rayleigh_diagnostics_alt import sliceinfo
from compute_grid_info import compute_grid_info, compute_theta_grid,\
        compute_r_grid

# handy class for making dictionaries "dot-accessible" "key-accessible" and vice versa
class dotdict(dict):
    """dot.notation access to dictionary attributes"""
    __getattr__ = dict.get
    __setattr__ = dict.__setitem__
    __delattr__ = dict.__delitem__


###############################
# BASIC ASTROPHYSICAL CONSTANTS
###############################

# universal gravitational constant
guniv= 6.67e-8      # Rounded to two decimals in Rayleigh...

# solar stuff
sun = dotdict(dict({}))

# Solar radius, luminosity, and mass (as we have been assuming in Rayleigh)
rsun = sun.r = 6.957e10  # value taken from IAU recommendation: arxiv, 1510.07674
                 # should probably figure out how Nick chose 
                 # 5.000 and 6.586209 as the base of the CZ and location 
                 # of 3rd density scale height (Comparing polytrope to 
                 # model S?)
lsun = sun.l = 3.846e33  # Used in Rayleigh: disagrees with IAU recommended value 
                 # of 3.828e33
msun = sun.m = 1.98891e33 # FROM WIKIPEDIA: 1.98847 \pm 0.00007
                  # From IAU recommendation: 1.9885, with 
                  # G = 6.67408 \pm 0.00031 (10^-8 c.g.s.)
                # NOTE: ALL THESE QUANTITIES CHANGE IN TIME (except G, if
                # the cosmologists are right...)

#Thermodyanmic variables
sun.cp = 3.5e8
gammaideal = 5./3.
sun.thermor = sun.cp*(1. - 1./gammaideal)

# solar base of the convection zone stuff
sun.rhobcz = 0.18053428
sun.tempbcz = 2111256.4
sun.rbcz = 5.0e10
sun.rbcznond = sun.rbcz/sun.r
# radius for the third density scale height (according to model S)
sun.rnrho3 = 6.5860209e10 

#######################################
# RANDOM CONSTANTS FOR COMPUTE ROUTINES
#######################################

# width of print messages in parallel routines
lent = 50
buff_frac = 0.05 # default buffer to set axes limits

# character to make text bold
bold_char_begin = "\033[1m"
bold_char_end = "\033[0m"

letters = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l',\
    'm', 'n', 'o', 'p', 'q', 'r', 's', 't', 'u', 'v', 'w', 'x', 'y', 'z']

buff_line = buffer_line = "=============================================="

# default latitudes to sample Rayleigh sim
default_latvals = np.array([-85., -75., -60., -45., -30., -15., 0., 15., 30., 45., 60., 75., 85.])
# default radial levels within a given zone 
# (samples basically every 12.5% of the way through, and at the top)
basedepths = np.array([0.05, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 0.95, 1.0])

######################################
# FORMATTING ROUTINES FOR PRINT OUTPUT
######################################

def make_bold(st):
    return bold_char_begin + st + bold_char_end

def format_time(seconds):
    if seconds < 1: # output in milliseconds
        return "%5i ms" %(round(seconds*1000))
    elif seconds < 10:
        # out put in decimals
        return "%1.2f s" %seconds
    elif seconds < 3600: # output as MM:SS
        seconds = int(seconds)
        minutes = seconds // 60
        seconds %= 60
        return "%02i:%02i" %(minutes, seconds)
    else: # output as HH:MM:SS
        seconds = int(seconds)
        hours = seconds // 3600
        seconds %= 3600
        minutes = seconds // 60
        seconds %= 60
        return "%02i:%02i:%02i" %(hours, minutes, seconds)

def format_size(nbytes):
    if nbytes < 1024: # output in B
        return "%i B" %nbytes
    elif nbytes < 1024**2:
        # output K
        return "%1.1f K" %(nbytes/1024)
    elif nbytes < 1024**3:
        # output M
        return "%1.1f M" %(nbytes/1024**2)
    else:
        # output G
        return "%1.1f G" %(nbytes/1024**3)

def fill_str(stri, lent, char):
    len_loc = len(stri)
    nfill = lent - len_loc
    return stri + char*nfill

def make_array(arr, tolist=False, length=None): 
    # arr is scalar, list, or array
    # convert everything to a list
    if arr is None:
        return None
    elif np.isscalar(arr):
        out = [arr]
    else:
        out = list(arr)

    # if length was specified (only good if OG length was 1
    # make the list that length
    if not length is None and len(out) == 1:
        out = out*length

    if tolist:
        return out
    else:
        return (np.array(out))

def lat_format(latval):
    if latval < 0:
        hemisphere = 'S'
    else:
        hemisphere = 'N'
    return hemisphere + '%02.0f' %np.abs(latval)

def array_of_strings(arr):
    li = []
    for ele in arr:
        li.append(str(ele))
    return np.array(li)

def arr_to_str(a, fmt):
    st = ''
    for ele in a:
        st += (fmt + ' ') %ele
    return '[' + st[:-1] + ']'

def update_dict(dict_orig, dict_update):
    dict_out = {**dict_orig} # start with original dict
    for key, val in dict_update.items():
        if key in dict_orig: # only update arguments that are
            # specified in the default set
            dict_out[key] = val
    return dotdict(dict_out)

def print_dict(di):
    for key, value in di.items():
        print (key, '\t\t', value)

def find_bad_keys(dict_orig, dict_update, funcname, justwarn=False):
    bad_keys = []
    for key, val in dict_update.items():
        if not key in dict_orig: 
            bad_keys.append(key)
    if len (bad_keys) > 0:
        if justwarn:
            print ("WARNING!")
        else:
            print ("ERROR!")
        print (funcname + '() got bad keyword args:')
        print (bad_keys)
        if not justwarn:
            print ("exiting")
            sys.exit()

########################
# BASIC COMPUTE ROUTINES
########################

def inds_from_vals(arr, arrvals):
    nind = len(arrvals)
    indarr = np.zeros(nind, 'int')
    for i in range(nind):
        indarr[i] = np.argmin(np.abs(arr - arrvals[i]))
    return indarr

def is_an_int(string):
    # obviously, first check if it's actually an int
    if isinstance(string, int) or isinstance(string, np.int64):
        return True
    len_str = len(string)
    bool_val = True
    for i in range(len_str):
        char = string[i]
        bool_val *= (char >= '0' and char <= '9')
    return(bool(bool_val))

def isall(arg):
    if np.isscalar(arg):
        isitall = arg
    else:
        isitall = args[0]

    if arg == 'all':
        return True
    else:
        return False

def frac_nonzero(arr):
    num_nonzero = len(np.where(arr != 0)[0])
    num_total = np.size(arr)
    return (num_nonzero/num_total)

def rms(array):
    if np.size(array) == 0:
        return 0
    else:
        return np.sqrt(np.mean(array**2))

# derivative routines
def drad(arr, rr): # this works for any dimension array, as long as
    # the radial index is the last one
    # make the output array
    darr = np.zeros_like(arr)

    # check for repeated radii; they are bad!
    nr = len(rr)
    ir_rep = []
    for ir in range(1, nr):
        if rr[ir] == rr[ir-1]:
            ir_rep.append(ir)

    ir1 = 0
    ir2 = nr - 1
    for ir_rep_loc in ir_rep:
        darr[..., ir1:ir_rep_loc] = np.gradient(arr[..., ir1:ir_rep_loc],\
                rr[ir1:ir_rep_loc], axis=arr.ndim-1)
        darr[..., ir_rep_loc] = darr[..., ir_rep_loc-1]
        ir1 = ir_rep_loc + 1
    # then do the last segment:
    darr[..., ir1:ir2+1] = np.gradient(arr[..., ir1:ir2+1],\
            rr[ir1:ir2+1], axis=arr.ndim-1)
    return darr

def dth(arr, tt): # assumes theta falls along the second-to-last axis
    # if array is > 1D
    if arr.ndim < 2:
        return np.gradient(arr, tt)
    else:
        return np.gradient(arr, tt, axis=arr.ndim-2)

def dph(arr): # assumes phi falls along first axis
    nphi = np.shape(arr)[0]
    dphi = 2.*np.pi/nphi
    return np.gradient(arr, dphi, axis=0)

def opt_workload(n, nproc):
    # optimally distributes workload (n tasks) over processes (n workers)
    n_per_proc_min = np.int(np.floor(n/nproc)) # min workload
    n_per_proc_max = np.int(np.ceil(n/nproc)) # max workload
    # min/max workloads differ by 1
    r = n/nproc - n_per_proc_min # remainder: r sets optimal number of processes
    # to perform max workload
    nproc_max = np.int(np.floor(nproc*r))
    nproc_min = nproc - nproc_max # there are total nproc processes

    # "optimal choice" assumes partial processes; but processes are whole
    # correct nproc_max/min to make sure all n tasks are perofrmed
    n_real_life = nproc_min*n_per_proc_min + nproc_max*n_per_proc_max 
    diff = n - n_real_life
    if diff > 0:
        nproc_max += diff
        nproc_min -= diff
    else:
        nproc_max -= diff
        nproc_min += diff
    return (nproc_min, nproc_max, n_per_proc_min, n_per_proc_max)

# Thin out the arrays to not deal obscene quantities of data 
# (and unreadable "curves")
def thin_data(vals, ntot=None):
    if ntot is None: # do nothing
        return vals
    else:
        nx = np.shape(vals)[0]
        nskip = nx//ntot
        if not nskip in [0, 1]: #for ntot < 2*nx, do nothing
            vals_new = vals[::nskip]
        else:
            vals_new = vals
        return vals_new

# Nonlinear Fourier transforms
def my_nfft(times, arr, axis=0):
    # shift the times to lie in range -1/2, 1/2
    total_time = times[-1] - times[0]
    times_shift = (times - times[0])/total_time - 1/2
    # get equally spaced times
    times_eq = np.linspace(-1/2, 1/2, len(times))
    interpolant = interp1d(times_shift, arr, axis=axis)
    arr_interp = interpolant(times_eq)
    arr_fft = np.fft.fft(arr_interp, axis=axis)
    arr_fft = np.fft.fftshift(arr_fft, axes=axis)

    # may as well get frequencies here too
    delta_t = np.mean(np.diff(times))
    freq = np.fft.fftfreq(len(times), delta_t)
    freq = np.fft.fftshift(freq)
    # return everything
    return arr_fft, freq

def my_infft(times, arr_fft, axis=0):
    # undo all the good work of the FFT

    # undo the frequency shift
    arr_fft = np.fft.ifftshift(arr_fft, axes=axis)

    # undo the fft
    arr_interp = np.fft.ifft(arr_fft, axis=axis)

    #  calculate shifted times (to -1/2, 1/2 interval)
    total_time = times[-1] - times[0]
    times_shift = (times - times[0])/total_time - 1/2

    # get equally spaced times
    times_eq = np.linspace(-1/2, 1/2, len(times))

    # get the "interpolant of the undoing"
    interpolant = interp1d(times_eq, arr_interp, axis=axis)

    # get the original array --- backward interpolated onto
    # the original times
    arr_orig = interpolant(times_shift)

    # return the original array
    return arr_orig

#############################################################
# ROUTINES FOR CHARACTERIZING/SORTING RAYLEIGH RAW DATA FILES
#############################################################

# options for selecting a range of files (to average over, trace over, etc.)
range_options = ['--range', '--centerrange', '--leftrange', '--rightrange', '--iter', '--n', '--f', '--all']

def my_mkdir(dirname):
    if not os.path.isdir(dirname):
        os.makedirs(dirname)
    return dirname

def get_file_lists_all(radatadir):
    # Get all the file names in datadir and their integer counterparts
    try:
        if 'Spherical_3D' in strip_dirname(radatadir):
            file_list_long = os.listdir(radatadir)
            file_list_long.sort()
            file_list = []
            for fname in file_list_long:
                if '_0001' in fname:
                    file_list.append(fname[:-5])
        else:
            file_list = os.listdir(radatadir)
            file_list.sort()
        nfiles = len(file_list)
        file_list = np.array(file_list)
    except: # if this fails, the directory must not have existed
        nfiles = 0
        file_list = np.array([])
    
    # Create an integer file list
    int_file_list = np.zeros(len(file_list), dtype=int)
    for i in range(nfiles):
        int_file_list[i] = int(file_list[i])
    
    return file_list, int_file_list, nfiles

def get_file_lists(radatadir, args):
    # Get file names in datadir and their integer counterparts
    # (only the ones in the desired range determined by args)
    # all the "action" occurs in get_desired_range() function below

    # get all files
    file_list, int_file_list, nfiles = get_file_lists_all(radatadir)
    # get the desired range
    index_first, index_last = get_desired_range(int_file_list, args)
    # Remove parts of file lists we don't need
    file_list = file_list[index_first:index_last + 1]
    int_file_list = int_file_list[index_first:index_last + 1]
    nfiles = index_last - index_first + 1
    return file_list, int_file_list, nfiles

def get_desired_range(int_file_list, args):
    # Get first and last index (within the int_file_list) associated with the desired range
    nargs = len(args)
    nfiles = len(int_file_list)
    # By default, the range will always be the last 100 files:
    index_first, index_last = nfiles - 100, nfiles - 1

    # user can modify this default in a number of ways
    for i in range(nargs):
        arg = args[i]
        if arg in ['--range', '--centerrange', '--leftrange',\
                '--rightrange', '--iter']: # first arg will be iter no.
            # 'first' means first available file.
            # 'last' means last available file.
            desired_iter = args[i+1]
            if desired_iter == 'first':
                desired_iter = int_file_list[0]
            elif desired_iter == 'last':
                desired__iter = int_file_list[-1]
            else:
                desired_iter = int(desired_iter)
            index = np.argmin(np.abs(int_file_list - desired_iter))
        if arg in ['--centerrange', '--rightrange', '--leftrange']:
            # many options include an "ndatafiles" argument
            ndatafiles = int(args[i+2])
        if arg in ['--n', '--f']:
            ndatafiles = int(args[i+1])
        if arg == '--range': # average between two specific files
            index_first = index # first arg is first desired iter
            # also need last iter
            desired_iter = args[i+2]
            if desired_iter == 'first':
                desired_iter = int_file_list[0]
            elif desired_iter == 'last':
                desired_iter = int_file_list[-1]
            else:
                desired_iter = int(desired_iter)
            index_last = np.argmin(np.abs(int_file_list - desired_iter))
        elif arg == '--centerrange': #range centered around specific file
            if ndatafiles % 2 == 0: #ndatafiles is even
                index_first = index - ndatafiles//2 + 1
                index_last = index + ndatafiles//2
            else:  #ndatafiles is odd
                index_first = index - ndatafiles//2
                index_last = index + ndatafiles//2
        elif arg == '--leftrange': # range with specific file first
            index_first = index
            index_last = index + ndatafiles - 1
        elif arg == '--rightrange': # range with specific file last
            index_last = index
            index_first = index - ndatafiles + 1
        elif arg == '--n': 
            # range with certain no. files ending with the last
            index_last = nfiles - 1
            index_first = nfiles - ndatafiles
        elif arg == '--f': 
            # range with certain no. files starting with the first
            index_first = 0
            index_last = ndatafiles - 1
        elif arg == '--all': # all files
            index_first = 0
            index_last = nfiles - 1
        elif arg == '--iter': # just get 1 iter
            index_first = index_last = index
    # Check to see if either of the indices fall "out of bounds"
    # and if they do replace them with the first or last index
    if index_first < 0: 
        index_first = 0
    if index_last > nfiles - 1: 
        index_last = nfiles - 1
    # Return the desired indices
    return index_first, index_last

########################################################################
# ROUTINES FOR CHARACTERIZING/READING RAYLEIGH POST-PROCESSED DATA FILES
########################################################################

def get_widest_range_file(datadir, dataname):
    # Find the desired post-processed file(s) in the data directory. If there are 
    # multiple, by default choose the one with widest range as far
    # as the input raw data files are concerned (i.e., the largest last_iter - first_iter)
    # If there is no matching file (i.e., the data product has not been computed yet), 
    # return None
    if os.path.isdir(datadir):
        datafiles = os.listdir(datadir)
        specific_files = []
        for i in range(len(datafiles)):
            datafile = datafiles[i]
            if dataname == datafile.split('-')[0]:
                specific_files.append(datafile)

        ranges = []
        iters1 = []
        iters2 = []
        if len(specific_files) > 0:
            for specific_file in specific_files:
                iter1, iter2 = get_iters_from_file(specific_file)
                ranges.append(iter2 - iter1)
                iters1.append(iter1)
                iters2.append(iter2)
            
            ranges = np.array(ranges)
            iters1 = np.array(iters1)
            iters2 = np.array(iters2)
            
            inds_max_range = np.where(ranges == np.max(ranges))
            iters2_maxrange = iters2[inds_max_range]
            # By default, use the file closest to the end of the simulation
            ind = inds_max_range[0][np.argmax(iters2_maxrange)]
            return datadir + specific_files[ind]
    # if we reached this point, no file can be found
    return None

def get_dict(fname):
    # basic routine for reading post-processed data (either .pkl or .npy format)
    # update 10/07/22: maybe I don't need the "npy" format anymore?
    if fname[-3:] == 'npy':
        di = np.load(fname, encoding='latin1', allow_pickle=True).item()
    elif fname[-3:] == 'pkl':
        f = open(fname, 'rb')
        di = pickle.load(f)
        f.close()
    return di

def read_log(fname):
    # routine for reading Rayleigh logfiles (the files I redirect Rayleigh output to)
    f = open(fname, 'r')
    lines = f.readlines()
    f.close()
    iters = []
    iters_io = []
    delta_t = []
    delta_t_io = []
    iters_per_sec = []
    iters_per_sec_io = []
    ncpu = -1
    for i in range(len(lines)):
        line = lines[i]
        if i > 0:
            lastline = lines[i - 1]
        if 'NCPU' in line:
            split = line.split()
            ncpu = int(split[-1])
        elif "teration" in line:
            split = line.split()
            lensplit = len(split)
            if lensplit == 6:
                for i in range(lensplit):
                    if split[i] == 'Iteration:':
                        iters.append(int(split[i+1]))
                    elif split[i] == 'DeltaT:':
                        delta_t.append(float(split[i+1]))
                    elif split[i] == 'Iter/sec:':
                        iters_per_sec.append(float(split[i+1]))
            elif lensplit == 7:
                for i in range(lensplit):
                    if split[i] == 'Iteration':
                        iters.append(int(split[i+2]))
                    elif split[i] == 'DeltaT':
                        delta_t.append(float(split[i+2]))
            else:
                print("read_log(%s): unrecognized text file format")
            if "Creating" in lastline: # this was a line right after output
                iters_io.append(iters[-1])
                delta_t_io.append(delta_t[-1])
                if lensplit == 6:
                    iters_per_sec_io.append(iters_per_sec[-1])

    di_out = dict({'iters': np.array(iters), 'delta_t': np.array(delta_t), 'ncpu': ncpu, 'iters_io': np.array(iters_io), 'delta_t_io': np.array(delta_t_io)})
    if len(iters_per_sec) > 0:
        di_out['iters_per_sec'] = np.array(iters_per_sec)
        di_out['iters_per_sec_io'] = np.array(iters_per_sec_io)

    return di_out

def get_dataname_from_file(filename):
    # route does what it's named to do...gets the dataname associated with post-processed 
    # data product name
    just_file = filename.split('/')[-1] #strip the possible full path info
    return just_file.split('-')[0]

def get_datadir_from_file(filename):
    # route does what it's named to do...
    just_dir = filename.split('/')[:-1] #strip the filename away from full path
    datadir_rel = ''
    start_adding = False
    for subdir in just_dir:
        if subdir == 'data':
            start_adding = True

        if start_adding:
            datadir_rel += subdir + '/'
    return datadir_rel

def get_iters_from_file(filename):
    # route does what it's named to do...
    filename_end = filename.split('-')[-1][:-4] 
    # (gets the [iter1]_[iter2].pkl and removes the trailing .ext)
    iters_st = filename_end.split('_')
    iter1, iter2 = int(iters_st[0]), int(iters_st[1])
    return iter1, iter2
        
##################################################################################
# ROUTINES TO GET BASIC INFO ABOUT A RAYLEIGH SIMULATION, GIVEN THE SIM. DIRECTORY
##################################################################################

###################################
# Start with basic input parameters
###################################

def strip_dirname(dirname, wrap=False):
    dirname_stripped = dirname.split('/')[-1]
    if dirname_stripped == '':
        dirname_stripped = dirname.split('/')[-2]
    if dirname == '.':
        dirname_stripped = os.getcwd().split('/')[-1]
    if dirname == '..':
        orig_dir = os.getcwd()
        os.chdir('..')
        dirname_stripped = os.getcwd().split('/')[-1]
        os.chdir(orig_dir)

    ncut = 20
    if wrap and len(dirname_stripped) > ncut:
        # Split dirname_stripped into two lines if it is very long
        dirname_stripped = dirname_stripped[:ncut] + '\n' +\
                dirname_stripped[ncut:]
    return dirname_stripped

def get_parameter(dirname, parameter):
    # read a paremter from main_input
    f = open(dirname + '/main_input')
    lines = f.readlines()
    n = len(lines)

    # search for parameter, line-by-line
    for i in range(n):
        # process each line
        line = lines[i]

        # make lower case
        line = line.lower()

        # remove spaces and newline character (at the end of each line)
        line = line.replace(' ', '')
        line = line.replace('\n', '')

        if line == '':
            continue
        
        # remove possible trailing comma from line
        if line[-1] == ',':
            line = line[:-1]

        # line ready to process
        # only lines with "=" are relevant
        # ignore the ones commented out with !
        if '=' in line and line[0] != '!':
            if '!' in line:
                # there was a comment after the code statement
                # throw it away!
                excl_index = line.index('!')
                line = line[:excl_index]

            lhs, rhs = line.split('=')
            if parameter == lhs: # found the parameter!
                num_string = rhs
                return (string_to_number_or_array(num_string))

    # if we reached this point, nothing was returned
    if parameter in ['magnetism', 'use_extrema', 'rotation']:
        return False # if these weren't specified, they are false
    else:
        raise Exception('The parameter ' + parameter + ' was not\n' +\
                        'specified in run: ' + dirname + '. \n' +\
                        'exiting NOW\n')

#########################################################
# some parameters (like luminosity and domain_bounds) can 
# be specified using multiple keywords
# in main_input, so need special routines to extract them
#########################################################

def get_lum(dirname):
    # Make lstar = lsun unless otherwise specified in main_input
    try:
        # First see if we can get c_10 from equation_coefficients:
        try:
            eq = equation_coefficients()
            eq.read(dirname + '/equation_coefficients')
            lstar = eq.constants[9]
            print("Got luminosity from 'equation_coefficients' file")
        except: # otherwise get "luminosity" from main_input
            lstar = get_parameter(dirname, 'luminosity')
            print ("Got luminosity from 'main_input' file")
    except:
        lstar = lsun
        print ("Cannot find luminosity in either 'equation_coefficients'")
        print("or 'main_input' files. Setting luminosity to lsun.")
    return lstar

def compute_Prot(dirname):
    try:
        Om0 = get_parameter(dirname, 'angular_velocity')
    except:
        eq = equation_coefficients()
        eq.read(dirname + '/equation_coefficients')
        Om0 = eq.constants[0]/2.
    return 2*np.pi/Om0

def get_domain_bounds(dirname):
    try:
        # if one-domain, should be able to get rmin, rmax somehow
        # (no way I can think of to determine if it's one domain
        # without just trying it, seeing if I get error)
        try: # can set boundaries via rmin, rmax directly
            rmin, rmax = get_parameter(dirname, 'rmin'),\
                    get_parameter(dirname, 'rmax')
        except: # ... or can set boundaries via aspect ratio and
            # shell depth
            aspect_ratio = get_parameter(dirname, 'aspect_ratio')
            shell_depth = get_parameter(dirname, 'shell_depth')
            rmin = shell_depth/(1/aspect_ratio - 1)
            rmax = shell_depth/(1 - aspect_ratio)

        nr = get_parameter(dirname, 'n_r')
        domain_bounds = np.array([rmin, rmax])
        ncheby = np.array([nr])
    except: # otherwise things must be multi-domain
        domain_bounds = get_parameter(dirname, 'domain_bounds')
        ncheby = get_parameter(dirname, 'ncheby')
    return ncheby, domain_bounds

def compute_tdt(dirname, mag=False, visc=False, tach=False):
    # Returns computed diffusion time (in sec) across whole layer
    # If tach=True, return diffusion time across whole layer,
    # across CZ and across RZ (tuple of 3)
    # Read in the diffusion profile
    eq = get_eq(dirname)
    rr = eq.radius
    if mag:
        diff = eq.eta
    elif visc:
        diff = eq.nu
    else:
        diff = eq.kappa

    # Compute and return the diffusion time
    if tach:
        domain_bounds = get_parameter(dirname, 'domain_bounds')
        ri, rm, ro = domain_bounds
        rmid = 0.5*(ri + ro)
        rmidrz = 0.5*(ri + rm)
        rmidcz = 0.5*(rm + ro)

        irmidrz = np.argmin(np.abs(rr - rmidrz))
        irmidcz = np.argmin(np.abs(rr - rmidcz))
        irmid = np.argmin(np.abs(rr - rmid))

        diff_midrz = diff[irmidrz]
        diff_midcz = diff[irmidcz]
        diff_mid = diff[irmid]

        Hrz = rm - ri
        Hcz = ro - rm
        H = ro - ri

        return Hrz**2.0/diff_midrz, Hcz**2.0/diff_midcz, H**2.0/diff_mid
    else:
        ri, ro = np.min(rr), np.max(rr)
        rmid = 0.5*(ri + ro)
        irmid = np.argmin(np.abs(rr - rmid))
        diff_mid = diff[irmid]
        H = ro - ri
        return H**2.0/diff_mid

#################################################################################
# Get basic radial coefficients (grid info, reference state) associated with sim.
#################################################################################

def get_eq(dirname, fname='equation_coefficients'): 
    # return a human readable version of equation_coefficients
    # [dirname], either using equation_coefficients or 
    # transport/reference files
    if os.path.exists(dirname + '/' + fname):
        # by default, get info from equation_coefficients (if file exists)
        eq = equation_coefficients()
        eq.read(dirname + '/' + fname)
        #eq_hr = eq_human_readable(eq.nr)
        eq_hr = dotdict(dict({}))
        eq_hr.radius = eq_hr.rr = eq.radius
        eq_hr.density = eq.functions[0]
        eq_hr.rho = eq_hr.density
        eq_hr.dlnrho = eq.functions[7]
        eq_hr.d2lnrho = eq.functions[8]
        eq_hr.temperature = eq.functions[3]
        eq_hr.T = eq_hr.temperature
        eq_hr.pressure = sun.thermor*eq_hr.rho*eq_hr.T
        eq_hr.P = eq_hr.pressure
        eq_hr.dlnT = eq.functions[9]
        eq_hr.gravity = eq.functions[1]/eq_hr.rho*sun.cp
        eq_hr.g = eq_hr.gravity
        eq_hr.dsdr = eq.functions[13]
        eq_hr.Nsq = (eq_hr.g/sun.cp)*eq_hr.dsdr
        eq_hr.heating = eq.constants[9]*eq.functions[5]
        eq_hr.Q = eq_hr.heating
        eq_hr.nu = eq.constants[4]*eq.functions[2]
        eq_hr.dlnu = eq.functions[10]
        eq_hr.kappa = eq.constants[5]*eq.functions[4]
        eq_hr.dlnkappa = eq.functions[11]
        eq_hr.eta = eq.constants[6]*eq.functions[6] # these are built-in to
        eq_hr.dlneta = eq.functions[12] # equation_coefficients as "zero"
        eq_hr.lum = eq.constants[9]
    else:
        ref = ReferenceState(dirname + '/reference')
        eq_hr = eq_human_readable(ref.nr)

        eq_hr.radius = eq_hr.rr = ref.radius
        eq_hr.density = ref.density
        eq_hr.rho = eq_hr.density
        eq_hr.dlnrho = ref.dlnrho
        eq_hr.d2lnrho = ref.d2lnrho
        eq_hr.temperature = ref.temperature
        eq_hr.T = eq_hr.temperature
        eq_hr.dlnT = ref.dlnt
        eq_hr.pressure = sun.thermor*eq_hr.rho*eq_hr.T
        eq_hr.P = eq_hr.pressure
        eq_hr.gravity = ref.gravity
        eq_hr.g = eq_hr.gravity
        eq_hr.dsdr = ref.dsdr
        eq_hr.heating = eq_hr.rho*eq_hr.T*ref.heating
        eq_hr.Q = eq_hr.heating
        # 'transport' didn't always used to exist, so only read it if possible
        if os.path.exists(dirname + '/transport'):
            trans = TransportCoeffs(dirname + '/transport')
            eq_hr.nu = trans.nu
            eq_hr.dlnu = trans.dlnu
            eq_hr.kappa = trans.kappa
            eq_hr.dlnkappa = trans.dlnkappa
            try:
                eq_hr.eta = trans.eta
                eq_hr.dlneta = dlneta # this will fail for hydro cases
                # "trans" will not have attributes eta, dlneta
            except: # if it failed, just keep the arrays zero             
                pass # (magnetism = False)
        else:
            print ("get_eq(): neither 'equation_coefficients' nor 'transport' found")
            print ("nu, dlnu, etc. will be zero")
            eq_hr.nu = np.zeros_like(eq_hr.rr)
            eq_hr.dlnu = np.zeros_like(eq_hr.rr)
            eq_hr.kappa = np.zeros_like(eq_hr.rr)
            eq_hr.dlnkappa = np.zeros_like(eq_hr.rr)
        eq_hr.lum = get_parameter(dirname, 'luminosity')
    return eq_hr

def get_grid_info(dirname):
    # get basic grid info
    di_out = dict({})
    gi = GridInfo(dirname + '/grid_info', '')
    # 1D arrays
    di_out['rr'] = gi.radius
    di_out['shell_depth'] = np.max(gi.radius) - np.min(gi.radius)
    di_out['tt'] = gi.theta
    di_out['cost'] = gi.costheta
    di_out['sint'] = gi.sintheta
    di_out['cott'] = di_out['cost']/di_out['sint']
    di_out['tt_lat'] = (np.pi/2 - di_out['tt'])*180/np.pi
    di_out['phi'] = gi.phi
    di_out['lons'] = gi.phi*180./np.pi
    di_out['rw'] = gi.rweights
    di_out['tw'] = gi.tweights
    di_out['pw'] = gi.pweights
    # grid dimensions
    di_out['nr'] = gi.nr
    di_out['nt'] = gi.ntheta
    di_out['nphi'] = gi.nphi
    # 2D arrays (theta, r)
    di_out['tt_2d'] = di_out['tt'].reshape((di_out['nt'], 1))
    di_out['sint_2d'] = np.sin(di_out['tt_2d'])
    di_out['cost_2d'] = np.cos(di_out['tt_2d'])
    di_out['cott_2d'] = di_out['cost_2d']/di_out['sint_2d']
    di_out['tw_2d'] = di_out['tw'].reshape((di_out['nt'], 1))
    di_out['rr_2d'] = di_out['rr'].reshape((1, di_out['nr']))
    di_out['rw_2d'] = di_out['rw'].reshape((1, di_out['nr']))
    di_out['xx'] = di_out['rr_2d']*di_out['sint_2d']
    di_out['zz'] = di_out['rr_2d']*di_out['cost_2d']
    # 3D arrays (phi, theta, r)
    di_out['phi_3d'] = di_out['phi'].reshape((di_out['nphi'], 1, 1))
    di_out['pw_3d'] = di_out['pw'].reshape((di_out['nphi'], 1, 1))
    di_out['tt_3d'] = di_out['tt'].reshape((1, di_out['nt'], 1))
    di_out['sint_3d'] = np.sin(di_out['tt_3d'])
    di_out['cost_3d'] = np.cos(di_out['tt_3d'])
    di_out['cott_3d'] = di_out['cost_3d']/di_out['sint_3d']
    di_out['tw_3d'] = di_out['tw'].reshape((1, di_out['nt'], 1))
    di_out['rr_3d'] = di_out['rr'].reshape((1, 1, di_out['nr']))
    di_out['rw_3d'] = di_out['rw'].reshape((1, 1, di_out['nr']))
    di_out['xx_3d'] = di_out['rr_3d']*di_out['sint_3d']
    di_out['zz_3d'] = di_out['rr_3d']*di_out['cost_3d']
    return dotdict(di_out)


####################################
# Routines associated with grid info
####################################

def integrate_in_r(arr, dirname):
    # routine to integrate in radius (over each domain separately)
    di = get_grid_info(dirname)
    rw = di['rw']
    nr = len(rw)
    ndim = arr.ndim
    newshape = [nr]
    for i in range(ndim - 1):
        newshape = [1] + newshape
    rw_nd = rw.reshape(newshape)
    ncheby, domain_bounds = get_domain_bounds(dirname)
    ndomains = len(ncheby)
   
    # start the return with the fully averaged arr
    li = [np.sum(arr*rw_nd, axis=ndim-1)]

    # loop over the domains and integrate
    if ndomains > 1:
        ir2 = nr
        ir1 = ir2 - ncheby[0]
        for idom in range(ndomains):
            if idom > 0:
                ir2 -= ncheby[idom - 1]
                ir1 -= ncheby[idom]
            rw_loc = rw_nd[..., ir1:ir2]
            li.append(np.sum(arr[..., ir1:ir2]*rw_loc, axis=ndim-1))
    return li

def get_volumes(dirname):
    # routine to integrate in radius (over each domain separately)
    di = get_grid_info(dirname)
    rr = di['rr']
    nr = di['nr']
    fact = 4.0*np.pi/3.0
    ncheby, domain_bounds = get_domain_bounds(dirname)
    ndomains = len(ncheby)
   
    # start the return with the full volume
    li = [fact*(np.max(rr)**3.0 - np.min(rr)**3.0)]

    # loop over the domains and integrate
    if ndomains > 1:
        ir2 = nr
        ir1 = ir2 - ncheby[0]
        for idom in range(ndomains):
            if idom > 0:
                ir2 -= ncheby[idom - 1]
                ir1 -= ncheby[idom]
            r1 = rr[ir1]
            r2 = rr[ir2 - 1]
            li.append(fact*(r2**3.0 - r1**3.0))
    return li

def interpret_rvals(dirname, rvals):
    # interpret array of rvals (array of strings), some could have the special keywords, rmin, rmid, rmax
    # but otherwise assumed to be float

    # get the actual min/max/mid of the full domain:
    ncheby, domain_bounds = get_domain_bounds(dirname)
    rmin, rmax = np.min(domain_bounds), np.max(domain_bounds)
    rmid = 0.5*(rmin + rmax)

    gi = get_grid_info(dirname)
    rr = gi.rr
    rvals_out = []
    rvals = make_array(rvals)
    for rval in rvals:
        if rval == 'rmin':
            rvals_out.append(rmin)
        elif rval == 'rmid':
            rvals_out.append(rmid)
        elif rval == 'rmax':
            rvals_out.append(rmax)
        else:
            # get the closest r value to the one specified
            ind = np.argmin(np.abs(rr - float(rval)))
            rvals_out.append(rr[ind])   
    return np.array(rvals_out)

def get_default_rvals(dirname, rvals=None):
    # default sampling locations
    # to specify multiple domains, need to specify boundary points (rvals, in cm)
    ncheby, domain_bounds = get_domain_bounds(dirname)
    rmin, rmax = np.min(domain_bounds), np.max(domain_bounds)

    if rvals is None:
        rvals = np.array([rmin, rmax])
    else:
        rvals = interpret_rvals(dirname, rvals)
        rvals = np.sort(rvals) # easier if rvals are in ascending order

    rvals_out = np.array([])
    for ir in range(len(rvals) - 1): # rvals define len(rvals) - 1 domains
        rbot = rvals[ir]
        rtop = rvals[ir+1]
        if ir == 0: # for the first domain, include the bottom boundary
            rvals_to_add = np.hstack((rbot, rbot + (rtop - rbot)*basedepths))
        else:
            rvals_to_add = rbot + (rtop - rbot)*basedepths
        rvals_out = np.hstack((rvals_out, rvals_to_add))
    return rvals_out

def get_sliceinfo(dirname, datatype='Shell_Slices', fname=None):
    radatadir = dirname + '/' + datatype + '/'
    file_list, int_file_list, nfiles = get_file_lists_all(radatadir)
    if fname is None:
        fname = file_list[0]
    return sliceinfo(fname, path=radatadir)

############################################
# ROUTINES FOR TIME PARAMETERS OF SIMULATION
############################################

def get_time_unit(dirname):
    # get basic time unit of simulation (rotation period or diffusion time)
    rotation = get_parameter(dirname, 'rotation')
    if rotation:
        time_unit = compute_Prot(dirname)
        time_label = r'${\rm{P_{rot}}}$'
        simple_label = 'rotations'
    else:
        time_unit = compute_tdt(dirname)
        time_label = r'${\rm{TDT}}$'
        simple_label = 'TDT'
    return time_unit, time_label, rotation, simple_label

def translate_times(time, dirname, translate_from='iter'):
    # TO USE MUST HAVE G_Avgs_trace file
    # change between different time units (can translate from: iter, unit (default time unit), 
    # rotation period, thermal diffusion time, seconds

    # Get the baseline time unit
    time_unit, time_label, rotation, simple_label = get_time_unit(dirname)

    # Get the G_Avgs trace
    datadir = dirname + '/data/'
    the_file = get_widest_range_file(datadir, 'G_Avgs_trace')
    if the_file is None:
        print ("translate_times(): you need to have G_Avgs_trace file")
        print ("to use me! Exiting.")
        sys.exit()
    else:
        di = get_dict(the_file)

    # Get times and iters from trace file
    times = di['times']
    iters = di['iters']

    if translate_from == 'iter':
        ind = np.argmin(np.abs(iters - time))
    elif translate_from in ['unit', 'prot', 'tdt']:
        ind = np.argmin(np.abs(times/time_unit - time))
    elif translate_from == 'sec':
        ind = np.argmin(np.abs(times - time))

    val_sec = times[ind]
    val_iter = iters[ind]
    val_unit = times[ind]/time_unit

    return dotdict(dict({'val_sec': val_sec,'val_iter': val_iter,\
            'val_unit': val_unit, 'simple_label': simple_label}))

def get_time_string(dirname, iter1, iter2=None, oneline=False):
    # Get the time range in sec
    t1 = translate_times(iter1, dirname, translate_from='iter')['val_sec']
    if not iter2 == None:
        t2 = translate_times(iter2, dirname, translate_from='iter')['val_sec']

    # Get the baseline time unit
    time_unit, time_label, rotation, simple_label = get_time_unit(dirname)

    # set the averaging-interval label
    if rotation:
        fmt = '%.0f' # measure rotations
    else:
        fmt = '%.3f'

    if not iter2 is None:
        if oneline:
            time_string = (('t = ' + fmt + ' to ' + fmt) %(t1/time_unit, t2/time_unit)) + ' ' + time_label + ' ' + r'$(\Delta t$' + ' = ' + (fmt %((t2 - t1)/time_unit)) + ' ' + time_label + ')'
        else:
            time_string = (('t = ' + fmt + ' to ' + fmt) %(t1/time_unit, t2/time_unit)) + ' ' + time_label + '\n' + r'$\Delta t$' + ' = ' + (fmt %((t2 - t1)/time_unit)) + ' ' + time_label
    else:
        time_string = (('t = ' + fmt + ' ') %(t1/time_unit)) + time_label

    return time_string

# default variable names to plot on slices 
def get_default_varnames(dirname):
    magnetism = get_parameter(dirname, 'magnetism')
    varnames_default = ['vr', 'vt', 'vp', 'omr',\
                'omt', 'omp', 'sprime', 'pprime', 'ssph', 'psph']
    if magnetism:
        varnames_default += ['br', 'bt', 'bp', 'jr', 'jt', 'jp']
    varnames_default = np.array(varnames_default)    
    return varnames_default

##############################################################
# ROUTINES FOR FIELD AMPLITUDES, LENGTH SCALES AND NON-D NUMBERS
##############################################################

def field_amp(dirname, the_file=None):
    # Make empty dictionary for field-amplitude arrays
    di_out = dotdict(dict([]))

    # See if run is magnetic
    magnetism = get_parameter(dirname, 'magnetism')

    # First get density
    eq = get_eq(dirname)
    rho = eq.rho
    rr = eq.radius
    nr = len(rr)
    di_out['rr'] = rr
    di_out['nr'] = nr

    # Get data directory
    datadir = dirname + '/data/'

    # Try to read in the Shell_Avgs data
    try:
        if the_file is None: # default
            the_file = get_widest_range_file(datadir, 'Shell_Avgs')
        print ('field_amp(): reading ' + the_file)
        di = get_dict(the_file)
        di_out['iter1'], di_out['iter2'] = get_iters_from_file(the_file)
        vals = di['vals']
        lut = di['lut']

        # Read in velocity-squared of flows
        # get this from kinetic energy
        vsqr = 2.0*vals[:, 0, lut[402]]/eq.rho
        vsqt = 2.0*vals[:, 0, lut[403]]/eq.rho
        vsqp = 2.0*vals[:, 0, lut[404]]/eq.rho

        vsqr_fluc = 2.0*vals[:, 0, lut[410]]/eq.rho
        vsqt_fluc = 2.0*vals[:, 0, lut[411]]/eq.rho
        vsqp_fluc = 2.0*vals[:, 0, lut[412]]/eq.rho

        # inherently positive, but avoid slightly negative
        # (machine error) when 0
        vsqr_mean = np.abs(vsqr - vsqr_fluc)
        vsqt_mean = np.abs(vsqt - vsqt_fluc)
        vsqp_mean = np.abs(vsqp - vsqp_fluc)

        # get the velocity amplitudes
        di_out['vr'] = np.sqrt(vsqr)
        di_out['vt'] = np.sqrt(vsqt)
        di_out['vp'] = np.sqrt(vsqp)
        di_out['vpol'] = np.sqrt(vsqr + vsqt)
        di_out['vhor'] = np.sqrt(vsqt + vsqp)
        di_out['v'] = np.sqrt(vsqr + vsqt + vsqp)

        di_out['vrfluc'] = np.sqrt(vsqr_fluc)
        di_out['vtfluc'] = np.sqrt(vsqt_fluc)
        di_out['vpfluc'] = np.sqrt(vsqp_fluc)
        di_out['vpolfluc'] = np.sqrt(vsqr_fluc + vsqt_fluc)
        di_out['vhorfluc'] = np.sqrt(vsqt_fluc + vsqp_fluc)
        di_out['vfluc'] = np.sqrt(vsqr_fluc + vsqt_fluc + vsqp_fluc)

        di_out['vrmean'] = np.sqrt(vsqr_mean)
        di_out['vtmean'] = np.sqrt(vsqt_mean)
        di_out['vpmean'] = np.sqrt(vsqp_mean)
        di_out['vpolmean'] = np.sqrt(vsqr_mean + vsqt_mean)
        di_out['vhormean'] = np.sqrt(vsqt_mean + vsqp_mean)
        di_out['vmean'] = np.sqrt(vsqr_mean + vsqt_mean + vsqp_mean)

        # Read in enstrophy of convective flows
        omsqr = vals[:, 0, lut[314]] 
        omsqt = vals[:, 0, lut[315]]
        omsqp = vals[:, 0, lut[316]]

        omsqr_fluc = vals[:, 0, lut[317]] 
        omsqt_fluc = vals[:, 0, lut[318]]
        omsqp_fluc = vals[:, 0, lut[319]]

        # inherently positive, but avoid slightly negative
        # (machine error) when 0
        omsqr_mean = np.abs(omsqr - omsqr_fluc)
        omsqt_mean = np.abs(omsqt - omsqt_fluc)
        omsqp_mean = np.abs(omsqp - omsqp_fluc)

        # get the vorticity amplitudes
        di_out['omr'] = np.sqrt(omsqr)
        di_out['omt'] = np.sqrt(omsqt)
        di_out['omp'] = np.sqrt(omsqp)
        di_out['ompol'] = np.sqrt(omsqr + omsqt)
        di_out['omhor'] = np.sqrt(omsqt + omsqp)
        di_out['om'] = np.sqrt(omsqr + omsqt + omsqp)

        di_out['omrfluc'] = np.sqrt(omsqr_fluc)
        di_out['omtfluc'] = np.sqrt(omsqt_fluc)
        di_out['ompfluc'] = np.sqrt(omsqp_fluc)
        di_out['ompolfluc'] = np.sqrt(omsqr_fluc + omsqt_fluc)
        di_out['omhorfluc'] = np.sqrt(omsqt_fluc + omsqp_fluc)
        di_out['omfluc'] = np.sqrt(omsqr_fluc + omsqt_fluc + omsqp_fluc)

        di_out['omrmean'] = np.sqrt(omsqr_mean)
        di_out['omtmean'] = np.sqrt(omsqt_mean)
        di_out['ompmean'] = np.sqrt(omsqp_mean)
        di_out['ompolmean'] = np.sqrt(omsqr_mean + omsqt_mean)
        di_out['omhormean'] = np.sqrt(omsqt_mean + omsqp_mean)
        di_out['ommean'] = np.sqrt(omsqr_mean + omsqt_mean + omsqp_mean)

        # get some thermo amplitudes

        # full spherical moments
        ssq = vals[:, 1, lut[501]]
        dsdrsq = vals[:, 1, lut[507]]
        psq = vals[:, 1, lut[502]]
        dpdrsq = vals[:, 1, lut[508]]

        # spherical means (square them)
        ssq_sphmean = vals[:, 0, lut[501]]**2
        psq_sphmean = vals[:, 0, lut[502]]**2
        dsdrsq_sphmean = vals[:, 0, lut[507]]**2
        dpdrsq_sphmean = vals[:, 0, lut[508]]**2

        # spherical fluctuations
        # inherently positive, but avoid slightly negative
        # (machine error) when 0
        ssq_sphfluc = np.abs(ssq - ssq_sphmean)
        dsdrsq_sphfluc = np.abs(dsdrsq - dsdrsq_sphmean)
        psq_sphfluc = np.abs(psq - psq_sphmean)
        dpdrsq_sphfluc = np.abs(dpdrsq - dpdrsq_sphmean)

        # add thermo amplitudes to dictionary        
        di_out['s'] = np.sqrt(ssq)
        di_out['ssphfluc'] = np.sqrt(ssq_sphfluc)
        di_out['ssphmean'] = np.sqrt(ssq_sphmean)

        di_out['p'] = np.sqrt(psq)
        di_out['psphfluc'] = np.sqrt(psq_sphfluc)
        di_out['psphmean'] = np.sqrt(psq_sphmean)
        
        di_out['dsdr'] = np.sqrt(dsdrsq)
        di_out['dsdrsphfluc'] = np.sqrt(dsdrsq_sphfluc)
        di_out['dsdrsphmean'] = np.sqrt(dsdrsq_sphmean)

        di_out['dpdr'] = np.sqrt(dpdrsq)
        di_out['dpdrsphfluc'] = np.sqrt(dpdrsq_sphfluc)
        di_out['dpdrsphmean'] = np.sqrt(dpdrsq_sphmean)

        if magnetism:
            # Read in squared B-fields
            # get this from magnetic energy
            eightpi = 8.*np.pi
            bsqr = eightpi*vals[:, 0, lut[1102]]
            bsqt = eightpi*vals[:, 0, lut[1103]]
            bsqp = eightpi*vals[:, 0, lut[1104]]

            bsqr_fluc = eightpi*vals[:, 0, lut[1110]]
            bsqt_fluc = eightpi*vals[:, 0, lut[1111]]
            bsqp_fluc = eightpi*vals[:, 0, lut[1112]]

            # inherently positive, but avoid slightly negative
            # (machine error) when 0
            bsqr_mean = np.abs(bsqr - bsqr_fluc)
            bsqt_mean = np.abs(bsqt - bsqt_fluc)
            bsqp_mean = np.abs(bsqp - bsqp_fluc)

            di_out['br'] = np.sqrt(bsqr)
            di_out['bt'] = np.sqrt(bsqt)
            di_out['bp'] = np.sqrt(bsqp)
            di_out['bpol'] = np.sqrt(bsqr + bsqt)
            di_out['bhor'] = np.sqrt(bsqt + bsqp)
            di_out['b'] = np.sqrt(bsqr + bsqt + bsqp)

            di_out['brfluc'] = np.sqrt(bsqr_fluc)
            di_out['btfluc'] = np.sqrt(bsqt_fluc)
            di_out['bpfluc'] = np.sqrt(bsqp_fluc)
            di_out['bpolfluc'] = np.sqrt(bsqr_fluc + bsqt_fluc)
            di_out['bhorfluc'] = np.sqrt(bsqt_fluc + bsqp_fluc)
            di_out['bfluc'] = np.sqrt(bsqr_fluc + bsqt_fluc + bsqp_fluc)

            di_out['brmean'] = np.sqrt(bsqr_mean)
            di_out['btmean'] = np.sqrt(bsqt_mean)
            di_out['bpmean'] = np.sqrt(bsqp_mean)
            di_out['bpolmean'] = np.sqrt(bsqr_mean + bsqt_mean)
            di_out['bhormean'] = np.sqrt(bsqt_mean + bsqp_mean)
            di_out['bmean'] = np.sqrt(bsqr_mean + bsqt_mean + bsqp_mean)

            # Read in squared current densities
            jsqr = vals[:, 0, lut[1014]]
            jsqt = vals[:, 0, lut[1017]]
            jsqp = vals[:, 0, lut[1020]]

            jsqr_fluc = vals[:, 0, lut[1015]]
            jsqt_fluc = vals[:, 0, lut[1018]]
            jsqp_fluc = vals[:, 0, lut[1021]]

            # inherently positive, but avoid slightly negative
            # (machine error) when 0
            jsqr_mean = np.abs(jsqr - jsqr_fluc)
            jsqt_mean = np.abs(jsqt - jsqt_fluc)
            jsqp_mean = np.abs(jsqp - jsqp_fluc)

            di_out['jr'] = np.sqrt(jsqr)
            di_out['jt'] = np.sqrt(jsqt)
            di_out['jp'] = np.sqrt(jsqp)
            di_out['jpol'] = np.sqrt(jsqr + jsqt)
            di_out['jhor'] = np.sqrt(jsqt + jsqp)
            di_out['j'] = np.sqrt(jsqr + jsqt + jsqp)

            di_out['jrfluc'] = np.sqrt(jsqr_fluc)
            di_out['jtfluc'] = np.sqrt(jsqt_fluc)
            di_out['jpfluc'] = np.sqrt(jsqp_fluc)
            di_out['jpolfluc'] = np.sqrt(jsqr_fluc + jsqt_fluc)
            di_out['jhorfluc'] = np.sqrt(jsqt_fluc + jsqp_fluc)
            di_out['jfluc'] = np.sqrt(jsqr_fluc + jsqt_fluc + jsqp_fluc)

            di_out['jrmean'] = np.sqrt(jsqr_mean)
            di_out['jtmean'] = np.sqrt(jsqt_mean)
            di_out['jpmean'] = np.sqrt(jsqp_mean)
            di_out['jpolmean'] = np.sqrt(jsqr_mean + jsqt_mean)
            di_out['jhormean'] = np.sqrt(jsqt_mean + jsqp_mean)
            di_out['jmean'] = np.sqrt(jsqr_mean + jsqt_mean + jsqp_mean)

    except:
        print ("field_amplitudes(): need to compute Shell_Avgs time avg")
        print ("returning empty dictionary")

    # Return the dictionary 
    return di_out

def length_scales(dirname, the_file=None):
    # Make empty dictionary for length_scale arrays
    di_out = dotdict(dict([]))

    # See if run is magnetic
    magnetism = get_parameter(dirname, 'magnetism')

    # First get mixing length scale + grid
    eq = get_eq(dirname)
    hrho = -1./eq.dlnrho
    rr = eq.radius
    nr = len(rr)

    di_out['rr'] = rr
    di_out['nr'] = nr
    di_out.hrho = hrho

    # Get data directory
    datadir = dirname + '/data/'

    # Get field amplitudes
    print ('length_scales(): ', end='')
    fa = field_amp(dirname, the_file=the_file)

    # Compute lengthscales (from flows) and put them in dictionary

    # full flows
    di_out.v = fa.v/fa.om
    di_out.vhor = fa.vhor/fa.omr
    di_out.vpol = fa.vpol/fa.omp
    # NOTE: the following will only be accurate if v_phi is the strongest component of v
    di_out.vp = fa.vp/fa.ompol

    # mean flows
    di_out.vmean = fa.vmean/fa.ommean
    di_out.vhormean = fa.vhormean/fa.omrmean
    di_out.vpolmean = fa.vpolmean/fa.ompmean
    # NOTE: the following will only be accurate if v_phi is the strongest component of v
    di_out.vpmean = fa.vpmean/fa.ompolmean

    # fluc flows
    di_out.vfluc = fa.vfluc/fa.omfluc
    di_out.vhorfluc = fa.vhorfluc/fa.omrfluc
    di_out.vpolfluc = fa.vpolfluc/fa.ompfluc
    # NOTE: the following will only be accurate if v_phi is the strongest component of v
    di_out.vpfluc = fa.vpfluc/fa.ompolfluc

    # Compute lengthscales (from fields) and put them in dictionary
    if magnetism:
        # full fields
        di_out.b = fa.b/fa.j
        di_out.bhor = fa.bhor/fa.jr
        di_out.bpol = fa.bpol/fa.jp
        # NOTE: the following will only be accurate if b_phi is the strongest component of b
        di_out.bp = fa.bp/fa.jpol

        # mean fields
        di_out.bmean = fa.bmean/fa.jmean
        di_out.bhormean = fa.bhormean/fa.jrmean
        di_out.bpolmean = fa.bpolmean/fa.jpmean
        # NOTE: the following will only be accurate if b_phi is the strongest component of b
        di_out.bpmean = fa.bpmean/fa.jpolmean

        # fluc fields
        di_out.bfluc = fa.bfluc/fa.jfluc
        di_out.bhorfluc = fa.bhorfluc/fa.jrfluc
        di_out.bpolfluc = fa.bpolfluc/fa.jpfluc
        # NOTE: the following will only be accurate if b_phi is the strongest component of b
        di_out.bpfluc = fa.bpfluc/fa.jpolfluc

    # Return the dictionary 
    return di_out

def get_numbers(dirname, the_file=None, shell_depth=None):
    # output dictionary
    di = dotdict(dict({}))

    # See if run is rotating and/or magnetic 
    rotation = get_parameter(dirname, 'rotation')
    magnetism = get_parameter(dirname, 'magnetism')
    if rotation:
        Om0 = 2*np.pi/compute_Prot(dirname)

    # get equation coefs
    eq = get_eq(dirname)

    # get grid info
    gi = get_grid_info(dirname)
    if shell_depth is None:
        shell_depth = gi.shell_depth

    # get field amplitudes and length scales
    print ("get_numbers(): ", end='')
    fa = field_amp(dirname, the_file=the_file)
    print ("get_numbers(): ", end='')
    ls = length_scales(dirname, the_file=the_file)

    # get some numbers!
    # (first, prognostic ones)

    # Prandtl
    di.pr = eq.nu/eq.kappa
    if magnetism:
        di.prm = eq.nu/eq.eta

    if rotation:
        # Ekman
        di.ek = eq.nu/(2*eq.shell_depth**2*Om0)

        # Taylor
        di.ta = 1/di.ek**2

        # Buoyancy parameter (zero for convection only)
        di.bref = eq.Nsq/Om0**2

    # (then, diagnostic ones)

    return 

def nonD_numbers(dirname, rbcz=None):
    # all the nonD numbers (as functions of radius and in different zones)
    # we could ever want

    # Make empty dictionary for length_scale arrays
    di_out = dict([])

    # See if run is magnetic
    magnetism = get_parameter(dirname, 'magnetism')
    rotation = get_parameter(dirname, 'rotation')

    # get reference state
    eq = get_eq(dirname)
    rr = eq.radius
    nr = len(rr)
    #di_out['rr'] = rr
    #di_out['nr'] = nr

    di_amp = field_amp(dirname)
    di_len = length_scales(dirname)

    # get the reference state
    eq = get_eq(dirname)

    # get the Reynolds numbers
    shell_depth = di_len['shell_depth']
    hrho = di_len['L_rho']

    di_out['Re'] = di_amp['vamp']*shell_depth/eq.nu
    di_out['Re_fluc'] = di_amp['vfluc']*shell_depth/eq.nu
    di_out['Re_mean'] = di_amp['vmean']*shell_depth/eq.nu

    di_out['Rehrho'] = di_amp['vamp']*hrho/eq.nu
    di_out['Rehrho_fluc'] = di_amp['vfluc']*hrho/eq.nu
    di_out['Rehrho_mean'] = di_amp['vfluc']*hrho/eq.nu

    L_om = di_len['L_om']
    di_out['Revort'] = di_amp['vamp']*L_om/eq.nu
    di_out['Revort_fluc'] = di_amp['vfluc']*L_om/eq.nu
    di_out['Revort_mean'] = di_amp['vmean']*L_om/eq.nu

    # Read in the Shell_Spectra data
    datadir = dirname + '/data/'
    the_file = get_widest_range_file(datadir, 'Shell_Spectra')
    if the_file == '':
        have_spec = False
    else: 
        have_spec = True

    if have_spec:
        ir_spec = di_len['ir_spec']
        rr_spec = di_len['rr_spec']
        L_v = di_len['L_v']
        di_out['Respec'] = (di_amp['vamp']/eq.nu)[ir_spec]*L_v
        di_out['Respec_fluc'] = (di_amp['vfluc']/eq.nu)[ir_spec]*L_v
        di_out['Respec_mean'] = (di_amp['vmean']/eq.nu)[ir_spec]*L_v

    if magnetism: # magnetic Reynolds numbers Rm
        L_J = di_len['L_J']
        di_out['Rm'] = di_amp['vamp']*L_J/eq.eta
        di_out['Rm_fluc'] = di_amp['vfluc']*L_J/eq.eta
        di_out['Rm_mean'] = di_amp['vmean']*L_J/eq.eta

        if have_spec:
            L_B = di_len['L_B']
            di_out['Rmspec'] = (di_amp['vamp']/eq.eta)[ir_spec]*L_B
            di_out['Rmspec_fluc'] = (di_amp['vfluc']/eq.eta)[ir_spec]*L_B
            di_out['Rmspec_mean'] = (di_amp['vmean']/eq.eta)[ir_spec]*L_B

    if rotation: # Rossby numbers
        Om0 = 2*np.pi/compute_Prot(dirname)
        di_out['Ro'] = di_amp['vamp']/(2.0*Om0*shell_depth)
        di_out['Ro_fluc'] = di_amp['vfluc']/(2.0*Om0*shell_depth)
        di_out['Ro_mean'] = di_amp['vmean']/(2.0*Om0*shell_depth)

        di_out['Rohrho'] = di_amp['vamp']/(2.0*Om0*hrho)
        di_out['Rohrho_fluc'] = di_amp['vfluc']/(2.0*Om0*hrho)
        di_out['Rohrho_mean'] = di_amp['vmean']/(2.0*Om0*hrho)

        di_out['Rovort'] = di_amp['vamp']/(2.0*Om0*L_om)
        di_out['Rovort_fluc'] = di_amp['vfluc']/(2.0*Om0*L_om)
        di_out['Rovort_mean'] = di_amp['vmean']/(2.0*Om0*L_om)

        if have_spec:
            di_out['Rospec'] = (di_amp['vamp']/eq.eta)[ir_spec]/(2.0*Om0*L_v)
            di_out['Rospec_fluc'] = (di_amp['vfluc']/eq.eta)[ir_spec]/(2.0*Om0*L_v)
            di_out['Rospec_mean'] = (di_amp['vmean']/eq.eta)[ir_spec]/(2.0*Om0*L_v)

    # now compute the global average of all numbers
    gi = GridInfo(dirname + '/grid_info', '')
    rw = gi.rweights
    if not rbcz is None:
        irbcz = np.argmin(np.abs(rr/rsun - rbcz))
        if have_spec:
            irbcz_spec = np.argmin(np.abs(rr_spec/rsun - rbcz))
        if not (irbcz == 0 or irbcz == nr - 1):
            rwcz = rw[:irbcz+1]/np.sum(rw[:irbcz+1])
            rwrz = rw[irbcz+1:]/np.sum(rw[irbcz+1:])
        else:
            print ('nonD_numbers(): dude, you entered a stupid value for')
            print ('rbcz. you set rbcz = %1.3e' %rbcz)
            print ('it needs be in the range [%.3f, %.3f]' %(np.min(rr)/rsun, np.max(rr)/rsun))
            print ('resetting rbcz = None')
            rbcz = None

    all_keys = list(di_out.keys())
    for key in all_keys:
        if 'spec' in key:
            di_out[key + '_gav'] = np.mean(di_out[key])
        else:
            di_out[key + '_gav'] = np.sum(di_out[key]*rw)
        if not rbcz is None:
            if 'spec' in key:
                if not (irbcz_spec == 0 or irbcz_spec == len(rr_spec) - 1):
                    di_out[key + '_cz'] = np.mean(di_out[key][:irbcz_spec+1])
                    di_out[key + '_rz'] = np.mean(di_out[key][irbcz_spec+1:])
                else:
                    di_out[key + '_cz'] = di_out[key]
                    di_out[key + '_rz'] = di_out[key]
            else:
                di_out[key + '_cz'] = np.sum(di_out[key][:irbcz+1]*rwcz)
                di_out[key + '_rz'] = np.sum(di_out[key][irbcz+1:]*rwrz)
    # I think we got it all!
    return di_out

