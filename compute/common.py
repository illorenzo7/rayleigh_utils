# Author: Loren Matilsky
# Created: 11/10/2018
# Common routines for routines for post-processing Rayleigh data 

import numpy as np
from scipy.interpolate import interp1d
import sys, os, pickle
from string_to_num import string_to_number_or_array
sys.path.append(os.environ['rapp'])

from reference_tools import equation_coefficients
from rayleigh_diagnostics import G_Avgs, Shell_Slices, GridInfo
from rayleigh_diagnostics_alt import sliceinfo
from grid_info import compute_grid_info

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
g_univ= 6.67e-8      # Rounded to two decimals in Rayleigh...

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
                 # Update 10/24/2022: Camissasa & Featherstone 2022 uses 3.839e33...revisit this
msun = sun.m = 1.98891e33 # FROM WIKIPEDIA: 1.98847 \pm 0.00007
                  # From IAU recommendation: 1.9885, with 
                  # G = 6.67408 \pm 0.00031 (10^-8 c.g.s.)
                # NOTE: ALL THESE QUANTITIES CHANGE IN TIME (except G, if
                # the cosmologists are right...)

#Thermodyanmic variables
sun.c_p = 3.5e8
gamma_ideal = 5./3.
sun.gas_constant = sun.c_p*(1. - 1./gamma_ideal)

# solar base of the convection zone stuff
sun.rho_bcz = 0.18053428
sun.temp_bcz = 2111256.4
sun.rbcz = 5.0e10
sun.rbcz_nond = sun.rbcz/sun.r
# radius for the third density scale height (according to model S)
sun.r_nrho3 = 6.5860209e10 

#######################################
# RANDOM CONSTANTS FOR COMPUTE ROUTINES
#######################################

# width of print messages in parallel routines
lenfill = 50
charfill = "."
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
base_depths = np.array([0.05, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 0.95, 1.0])

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

def fill_str(st, lenfill=lenfill, charfill=charfill):
    len_st = len(st)
    nfill = lenfill - len_st
    return st + charfill*nfill

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
    # read a parameter from main_input

    # first read in main_input
    f = open(dirname + '/main_input')
    lines = f.readlines()

    # pare down this down a bit
    lines_new = []
    for i in range(len(lines)):
        # process each line
        line = lines[i]

        # make lower case
        line = line.lower()

        # remove spaces and newline/tab characters
        line = line.replace(' ', '')
        line = line.replace('\n', '')
        line = line.replace('\t', '')

        # some lines will be blank:
        if line != '':
            # only keep line if it's not a comment
            if line[0] != "!":
                # remove trailing comments
                i_exclamation = line.find('!')
                if i_exclamation != -1: # find returns -1 if character 
                        # wasn't found (i.e. there were no comments)
                    line = line[:i_exclamation]
                lines_new.append(line)

    # see if parameter is there --- search line by line
    st_param = ''
    for i in range(len(lines_new)):
        line = lines_new[i]
        if line.split('=')[0] == parameter: # found the parameter
            st_param += line.split('=')[1]
            # sometimes values continue on next line(s)
            keep_searching = True
            for j in range(i+1, len(lines_new)):
                if "=" in lines_new[j] or lines_new[j] == '/':
                    keep_searching = False
                if keep_searching:
                    st_param += lines_new[j]
    if st_param != '': # we successfully read the parameter
        return string_to_number_or_array(st_param)
    else: # we couldn't find parameter explicitly

        # these parameters still have a value (False)
        # if they weren't specified
        if parameter in ['magnetism', 'rotation']:
            return False # if these weren't specified, they are false
        else: # or finally, the parameter might just not be there
            return None

#########################################################
# some parameters (like luminosity and domain_bounds) can 
# be specified using multiple keywords
# in main_input, so need special routines to extract them
#########################################################

def get_domain_bounds(dirname):
    # get radial domain bounds associated with simulation

    # see if things are multi-domain
    domain_bounds = get_parameter(dirname, 'domain_bounds')
    ncheby = get_parameter(dirname, 'ncheby')

    if domain_bounds is None or ncheby is None: # need to determine these
        # a different way

        # if one-domain, should be able to get rmin, rmax somehow
        rmin, rmax = get_parameter(dirname, 'rmin'),\
                get_parameter(dirname, 'rmax')
        if rmin is None or rmax is None: # if still None, one more option
            # ... or can set boundaries via aspect ratio shell depth
            aspect_ratio = get_parameter(dirname, 'aspect_ratio')
            shell_depth = get_parameter(dirname, 'shell_depth')
            rmin = shell_depth/(1/aspect_ratio - 1)
            rmax = shell_depth/(1 - aspect_ratio)

        nr = get_parameter(dirname, 'n_r')
        domain_bounds = np.array([rmin, rmax])
        ncheby = np.array([nr])
    return ncheby, domain_bounds

#################################################################################
# Get basic radial coefficients (grid info, reference state) associated with sim.
#################################################################################

####################################
# Routines associated with grid info
####################################
def get_grid_info(dirname):
    # get basic grid info; try to read from grid_info file
    # or directly from main_input if grid_info doesn't exist
    di = dotdict()

    # get basic grid (colocation points and weights)
    if os.path.exists(dirname + '/grid_info'):
        gi = GridInfo(dirname + '/grid_info', '')
        # 1D arrays
        di.rr = rr = gi.radius
        di.rw = rr = gi.rweights
        di.tt = tt = gi.theta
        di.tw = tw = gi.tweights
    else:
        ncheby, domain_bounds = get_domain_bounds(dirname)
        nt = get_parameter(dirname, 'n_theta')
        out = compute_grid_info(ncheby, domain_bounds, nt)
        di.rr = rr = out[0]
        di.rw = rr = out[1]
        di.tt = tt = out[2]
        di.tw = tw = out[3]

    # some derivative theta (tt) quantities 
    di['cost'] = np.cos(tt)
    di['sint'] = np.sin(tt)
    di['cott'] = di['cost']/di['sint']
    di['tt_lat'] = (np.pi/2 - di['tt'])*180/np.pi

    # grid dimensions
    di['nr'] = nr =  len(rr)
    di['nt'] = nt = len(tt)
    di['nphi'] = nphi = 2*nt

    # phi stuff
    di['phi'] = phi = np.linspace(0, 2*np.pi, nphi, endpoint=False)
    di['lons'] = di.phi*180./np.pi
    di['pw'] = 1./nphi + np.zeros(nphi)

    # 2D arrays (theta, r)
    di['tt_2d'] = di['tt'].reshape((di['nt'], 1))
    di['sint_2d'] = np.sin(di['tt_2d'])
    di['cost_2d'] = np.cos(di['tt_2d'])
    di['cott_2d'] = di['cost_2d']/di['sint_2d']
    di['tw_2d'] = di['tw'].reshape((di['nt'], 1))
    di['rr_2d'] = di['rr'].reshape((1, di['nr']))
    di['rw_2d'] = di['rw'].reshape((1, di['nr']))
    di['xx'] = di['rr_2d']*di['sint_2d']
    di['zz'] = di['rr_2d']*di['cost_2d']
    # 3D arrays (phi, theta, r)
    di['phi_3d'] = di['phi'].reshape((di['nphi'], 1, 1))
    di['pw_3d'] = di['pw'].reshape((di['nphi'], 1, 1))
    di['tt_3d'] = di['tt'].reshape((1, di['nt'], 1))
    di['sint_3d'] = np.sin(di['tt_3d'])
    di['cost_3d'] = np.cos(di['tt_3d'])
    di['cott_3d'] = di['cost_3d']/di['sint_3d']
    di['tw_3d'] = di['tw'].reshape((1, di['nt'], 1))
    di['rr_3d'] = di['rr'].reshape((1, 1, di['nr']))
    di['rw_3d'] = di['rw'].reshape((1, 1, di['nr']))
    di['xx_3d'] = di['rr_3d']*di['sint_3d']
    di['zz_3d'] = di['rr_3d']*di['cost_3d']
    return di


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
    return np.sort(np.array(rvals_out))

def get_sliceinfo(dirname, datatype='Shell_Slices', fname=None):
    radatadir = dirname + '/' + datatype + '/'
    file_list, int_file_list, nfiles = get_file_lists_all(radatadir)
    if fname is None:
        fname = file_list[0]
    return sliceinfo(fname, path=radatadir)

def get_vol(dirname, r1='rmin', r2='rmax'):
    # get the shell volume in the range (r1, r2)
    r1, r2 = interpret_rvals(dirname, np.array([r1, r2]))
    return 4./3.*np.pi*(r2**3 - r1**3)

def volav_in_radius(dirname, arr, r1='rmin', r2='rmax'):
    gi = get_grid_info(dirname)
    rr, rw = gi.rr, gi.rw
    r1, r2 = interpret_rvals(dirname, np.array([r1, r2]))
    ir1, ir2 = np.argmin(np.abs(rr - r1)), np.argmin(np.abs(rr - r2))
    return np.sum((arr*rw)[ir2:ir1+1])/np.sum(rw[ir2:ir1+1])

##########################################
# Routines associated with reference state
##########################################

def compute_polytrope(rmin=sun.rbcz, rmax=sun.r_nrho3, nr=500, rr=None, poly_nrho=3.0, poly_n=1.5, poly_rho_i=sun.rho_bcz, poly_mass=sun.m, c_p=sun.c_p):
    # compute polytrope from N_rho and inner density, following 
    # Jones et al. (2011)

    di = dotdict()

    if rr is None:
        di.rr = rr = np.linspace(rmax, rmin, nr)
    else:
        di.rr = rr
        rmin, rmax = np.min(rr), np.max(rr)
        nr = len(rr)

    d = rmax - rmin
    beta = rmin/rmax
    gas_constant = c_p/(poly_n + 1)
    poly_gamma = (poly_n + 1)/poly_n

    exp = np.exp(poly_nrho/poly_n)

    zeta_o = (beta + 1)/(beta*exp + 1)
    zeta_i = (1 + beta - zeta_o)/beta
    c0 = (2*zeta_o - beta - 1)/(1 - beta)
    c1 = (1 + beta)*(1 - zeta_o)/(1 - beta)**2
    zeta = c0 + c1*d/rr

    tmp_c = g_univ*poly_mass/(c_p*c1*d)
    rho_c = poly_rho_i/zeta_i**poly_n
    prs_c = gas_constant*rho_c*tmp_c

    di.rho = rho_c*zeta**(poly_n)
    di.prs = prs_c*zeta**(poly_n+1.)
    di.tmp = tmp_c*zeta
    di.grav = g_univ*poly_mass/rr**2
    
    dzeta = -c1*d/rr**2
    d2zeta = 2*c1*d/rr**3
    dlnzeta = dzeta/zeta
    d2lnzeta = d2zeta/zeta - dzeta**2/zeta**2 

    di.dlnrho = poly_n*dlnzeta
    di.d2lnrho = poly_n*d2lnzeta
    di.dlntmp = dlnzeta

    dlnprs = di.dlnrho + di.dlntmp
    di.dsdr = c_p*(dlnprs/poly_gamma - di.dlnrho)

    return di

def get_eq(dirname, fname=None): 
    # return a human readable version of equation_coefficients
    # for [dirname], using either (in order of priority)

    # 1. binary file specified by fname
    # 2. equation_coefficients file
    # 3. custom_reference_binary file
    # 4. polytrope + transport coefs defined by main_input

    # human readable equation coefficients object to store reference state
    eq_hr = dotdict()

    # need this no matter what -- not sure if it works for non-D anelastic
    # (i.e. may set dissipation number instead)
    eq_hr.c_p = get_parameter(dirname, 'pressure_specific_heat')

    # and this no matter what
    magnetism = get_parameter(dirname, 'magnetism')

    # read reference state from binary file or main_input
    if fname is None:
        if os.path.exists(dirname + '/' + 'equation_coefficients'):
            fname = 'equation_coefficients'
        elif os.path.exists(dirname + '/' + 'custom_reference_binary'):
            fname = 'custom_reference_binary'

    if fname is None: # no binary file; get everything from main_input
        gi = get_grid_info(dirname)
        eq_hr.rr = gi.rr 
        zero = np.zeros_like(eq_hr.rr)

        poly_nrho = get_parameter(dirname, 'poly_nrho')
        poly_n = get_parameter(dirname, 'poly_n')
        poly_rho_i = get_parameter(dirname, 'poly_rho_i')
        poly_mass = get_parameter(dirname, 'poly_mass')
        poly = compute_polytrope(rr=eq_hr.rr, poly_nrho=poly_nrho, poly_n=poly_n, poly_rho_i=poly_rho_i, poly_mass=poly_mass, c_p=eq_hr.c_p)

        # get the background reference state
        eq_hr.rho = poly.rho 
        eq_hr.dlnrho = poly.dlnrho
        eq_hr.d2lnrho = poly.d2lnrho
        eq_hr.tmp = poly.tmp
        eq_hr.prs = poly.prs
        eq_hr.dlntmp = poly.dlntmp
        eq_hr.grav = poly.grav
        eq_hr.dsdr = poly.dsdr
        eq_hr.nsq = poly.dsdr*eq_hr.grav/eq_hr.c_p

        # get heating
        eq_hr.lum = get_parameter(dirname, 'luminosity')
        heating_type = get_parameter(dirname, 'heating_type')
        if heating_type is None: # default heating_type = 0
            eq_hr.heat = zero
        elif heating_type == 1: # heating proportional to pressure
                                # "consant entropy heating"
            eq_hr.heat = eq_hr.prs - eq_hr.prs[0]
            integral = volav_in_radius(dirname, eq_hr.heat)
            integral *= get_vol(dirname)
            eq_hr.heat = eq_hr.heat*eq_hr.lum/integral
        elif heating_type == 4: # constant energy heating
            eq_hr.heat = zero + eq_hr.lum/get_vol(dirname)

        # get the transport coefficients
        nu_top = get_parameter(dirname, 'nu_top')
        kappa_top = get_parameter(dirname, 'kappa_top')
        nu_type = get_parameter(dirname, 'nu_type')
        kappa_type = get_parameter(dirname, 'kappa_type')

        if nu_type is None or nu_type == 1:
            eq_hr.nu = zero + nu_top
            eq_hr.dlnu = zero
        elif nu_type == 2:
            nu_power = get_parameter(dirname, 'nu_power')
            eq_hr.nu = nu_top*(eq_hr.rho/eq_hr.rho[0])**nu_power
            eq_hr.dlnu = nu_power*eq_hr.dlnrho

        if kappa_type is None or kappa_type == 1:
            eq_hr.kappa = zero + kappa_top
            eq_hr.dlnkappa = zero
        elif kappa_type == 2:
            kappa_power = get_parameter(dirname, 'nu_power')
            eq_hr.kappa = kappa_top*(eq_hr.rho/eq_hr.rho[0])**kappa_power
            eq_hr.dlnkappa = kappa_power*eq_hr.dlnkappa
        
        if magnetism:
            if eta_type is None or eta_type == 1:
                eq_hr.eta = zero + eta_top
                eq_hr.dlneta = zero
            elif eta_type == 2:
                eta_power = get_parameter(dirname, 'eta_power')
                eq_hr.eta = eta_top*(eq_hr.rho/eq_hr.rho[0])**eta_power
                eq_hr.dlneta = eta_power*eq_hr.dlnrho

        # finally, get the rotation rate
        eq_hr.om0 = get_parameter(dirname, 'angular_velocity')
        eq_hr.prot = 2*np.pi/eq_hr.om0

    else:
        # by default, get info from equation_coefficients (if file exists)
        # get the background reference state
        eq = equation_coefficients()
        eq.read(dirname + '/' + fname)
        eq_hr.rr = eq.radius
        eq_hr.rho = eq.functions[0]
        eq_hr.dlnrho = eq.functions[7]
        eq_hr.d2lnrho = eq.functions[8]
        eq_hr.tmp = eq.functions[3]
        eq_hr.dlntmp = eq.functions[9]
        eq_hr.grav = eq.functions[1]/eq_hr.rho*eq_hr.c_p
        eq_hr.dsdr = eq.functions[13]
        eq_hr.heat = eq.constants[9]*eq.functions[5]
        eq_hr.lum = eq.constants[9]

        # assume gas is ideal and get pressure
        gas_constant = (gamma_ideal-1)*eq_hr.c_p/gamma_ideal
        eq_hr.prs = gas_constant*eq_hr.rho*eq_hr.tmp

        # get the transport coefficients
        eq_hr.nu = eq.constants[4]*eq.functions[2]
        eq_hr.dlnu = eq.functions[10]
        eq_hr.kappa = eq.constants[5]*eq.functions[4]
        eq_hr.dlnkappa = eq.functions[11]
        if magnetism:
            eq_hr.eta = eq.constants[6]*eq.functions[6] # these are built-in to
            eq_hr.dlneta = eq.functions[12] # equation_coefficients as "zero"

        # finally, get the rotation rate
        eq_hr.om0 = eq.constants[0]/2
        eq_hr.prot = 2*np.pi/eq_hr.om0

    # some derivative quantities

    # buoyancy frequency
    eq_hr.nsq = (eq_hr.grav/eq_hr.c_p)*eq_hr.dsdr

    # density scale height
    eq_hr.hrho = -1/eq_hr.dlnrho

    # thermal diffusion time
    rmin, rmax = np.min(eq_hr.rr), np.max(eq_hr.rr)
    irmid = np.argmin(np.abs(eq_hr.rr - (rmin + rmax)/2))
    eq_hr.tdt = (rmax - rmin)**2/eq_hr.kappa[irmid]

    return eq_hr


############################################
# ROUTINES FOR TIME PARAMETERS OF SIMULATION
############################################

# for backwards compatibility (remove references as I see them)
def compute_Prot(dirname):
    eq = get_eq(dirname)
    return eq.prot

def compute_tdt(dirname):
    eq = get_eq(dirname)
    return eq.tdt

def get_time_unit(dirname):
    # get basic time unit of simulation (rotation period or diffusion time)
    rotation = get_parameter(dirname, 'rotation')
    eq = get_eq(dirname)
    if rotation:
        time_unit = eq.prot
        time_label = r'${\rm{P_{rot}}}$'
        simple_label = 'rotations'
    else:
        time_unit = eq.tdt
        time_label = simple_label = 'TDT'
    return time_unit, time_label, rotation, simple_label

def translate_times(time, dirname, translate_from='iter'):
    # TO USE MUST HAVE G_Avgs_trace file
    # change between different time units (can translate from: 
    # iter, prot, tdt, sec

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

    # get the equation_coefficients file
    eq = get_eq(dirname)

    # translate the time
    if translate_from == 'iter':
        ind = np.argmin(np.abs(iters - time))
    elif translate_from == 'prot':
        ind = np.argmin(np.abs(times/eq.prot - time))
    elif translate_from == 'tdt':
        ind = np.argmin(np.abs(times/eq.tdt - time))
    elif translate_from == 'sec':
        ind = np.argmin(np.abs(times - time))

    # prepare the dictionary to return
    di = dotdict()
    di.val_sec = times[ind]
    di.val_iter = iters[ind]
    di.val_tdt = times[ind]/eq.tdt
    rotation = get_parameter(dirname, 'rotation')
    if rotation:
        di.val_prot = times[ind]/eq.prot
    return di

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

##############################################################
# ROUTINES FOR FIELD AMPLITUDES, LENGTH SCALES AND NON-D NUMBERS
##############################################################

def field_amp(dirname, the_file=None):
    # Make empty dictionary for field-amplitude arrays
    di_out = dotdict()

    # See if run is magnetic
    magnetism = get_parameter(dirname, 'magnetism')

    # First get density
    eq = get_eq(dirname)
    rho = eq.rho
    rr = eq.rr
    di_out['rr'] = rr

    # Get data directory
    datadir = dirname + '/data/'

    # Read in the Shell_Avgs data
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
    rr = eq.rr

    di_out['rr'] = rr
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
