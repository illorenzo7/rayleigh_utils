# Author: Loren Matilsky
# Created: 11/10/2018
# Common routines for routines for post-processing Rayleigh data 

import numpy as np
import os, pickle

# Solar radius, luminosity, and mass (as we have been assuming in Rayleigh)
rsun = 6.957e10  # value taken from IAU recommendation: arxiv, 1510.07674
                 # should probably figure out how Nick chose 
                 # 5.000 and 6.586209 as the base of the CZ and location 
                 # of 3rd density scale height (Comparing polytrope to 
                 # model S?)
G = 6.67e-8      # Rounded to two decimals in Rayleigh...
lsun = 3.846e33  # Used in Rayleigh: disagrees with IAU recommended value 
                 # of 3.828e33
msun = 1.98891e33 # FROM WIKIPEDIA: 1.98847 \pm 0.00007
                  # From IAU recommendation: 1.9885, with 
                  # G = 6.67408 \pm 0.00031 (10^-8 c.g.s.)
                # NOTE: ALL THESE QUANTITIES CHANGE IN TIME (except G, if
                # the cosmologists are right...)
# Read in all files from the Rayleigh data directory and sort them by name (number)

# I am now calling r_m the base of the convection zone, 
# while r_i (the inner shell radius) can vary
rhom = 0.18053428
Tm = 2111256.4
rm = 5.0e10
ro = 6.5860209e10 # Radii consistent with the bottom 3 density scale 
        # heights in the Sun rho_i above corresponds to the density
        # at the base of the convection zone

def get_file_lists(radatadir):
    file_list = os.listdir(radatadir)
    file_list.sort()
    nfiles = len(file_list)
    file_list = np.array(file_list)
    
    # Create an integer file list
    int_file_list = np.zeros(len(file_list), dtype=int)
    for i in range(nfiles):
        int_file_list[i] = int(file_list[i])
    
    return file_list, int_file_list, nfiles

range_options = ['-range', '-centerrange', '-leftrange', '-rightrange', '-n',\
               '-f', '-all', '-iter']

def get_desired_range(int_file_list, args):
    nargs = len(args)
    nfiles = len(int_file_list)
    # By default, the range will always be the last 100 files:
    index_first, index_last = nfiles - 100, nfiles - 1
    for i in range(nargs):
        arg = args[i]
        if arg == '-range':   
            # Give user option to specify range of iterations to
            # average over. 'first' means first available file.
            # 'last' means last available file.
            desired_first_iter = args[i+1]
            if desired_first_iter == 'first':
                desired_first_iter = int_file_list[0]
            elif desired_first_iter == 'last':
                # The only consistent option here would be
                # -range last last (to average over ONLY the last file)
                desired_first_iter = int_file_list[-1]
            else:
                desired_first_iter = int(desired_first_iter) 
                
            desired_last_iter = args[i+2]
            if desired_last_iter == 'last':
                desired_last_iter = int_file_list[-1]
            elif desired_last_iter == 'first':
                # only consistent option -range first first
                desired_last_iter = int_file_list[0]
            else:
                desired_last_iter = int(desired_last_iter)

            # Find the index in the int_file_list closest 
            # to the user's preferences
            index_first = np.argmin(np.abs(desired_first_iter -\
                    int_file_list))
            index_last = np.argmin(np.abs(desired_last_iter -\
                    int_file_list))
        elif arg == '-centerrange':
            desired_central_iter = args[i+1]
            if desired_central_iter == 'first':
                desired_central_iter = int_file_list[0]
            elif desired_central_iter == 'last':
                desired_central_iter = int_file_list[-1]
            else:
                desired_central_iter = int(desired_central_iter)
            central_index = np.argmin(np.abs(desired_central_iter -\
                    int_file_list))
            ndatafiles = int(args[i+2])
            if ndatafiles % 2 == 0: #ndatafiles is even
                index_first = central_index - ndatafiles//2 + 1
                index_last = central_index + ndatafiles//2
            else:  #ndatafiles is odd
                index_first = central_index - ndatafiles//2
                index_last = central_index + ndatafiles//2
        elif arg == '-leftrange':
            desired_left_iter = args[i+1]
            if desired_left_iter == 'first':
                desired_left_iter = int_file_list[0]
            elif desired_left_iter == 'last':
                # Only viable option here is if ndatafiles is 1
                desired_left_iter = int_file_list[-1]
            else:
                desired_left_iter = int(desired_left_iter)
            left_index = np.argmin(np.abs(desired_left_iter -\
                    int_file_list))
            ndatafiles = int(args[i+2])
            index_first = left_index
            index_last = left_index + ndatafiles - 1
        elif arg == '-rightrange':
            desired_right_iter = args[i+1]
            if desired_right_iter == 'first':
                desired_right_iter = int_file_list[0]
            elif desired_right_iter == 'last':
                # Only viable option here is if ndatafiles is 1
                desired_right_iter = int_file_list[-1]
            else:
                desired_right_iter = int(desired_right_iter)
            right_index = np.argmin(np.abs(desired_right_iter -\
                    int_file_list))
            ndatafiles = int(args[i+2])
            index_first = right_index - ndatafiles + 1
            index_last = right_index
        elif arg == '-n': 
            # allow the user to specify a desired number of iterations
            # to average over, ending with the last data file
            index_last = nfiles - 1
            number_to_average = int(args[i+1])
            index_first = nfiles - number_to_average
        elif arg == '-f': 
            # allow the user to specify a desired number of iterations
            # to average over, starting with the first data file
            index_first = 0
            number_to_average = int(args[i+1])
            index_last = number_to_average - 1
        elif arg == '-all':
            index_first = 0
            index_last = nfiles - 1
        elif arg == '-iter': # just get 1 iter
            desired_iter = args[i+1]
            if desired_iter == 'first':
                desired_iter = int_file_list[0]
            elif desired_iter == 'last':
                desired_iter = int_file_list[-1]
            else:
                desired_iter = int(desired_iter)
            index_first = np.argmin(np.abs(desired_iter - int_file_list))
            index_last = index_first
    # Check to see if either of the indices fall "out of bounds"
    # and if they do replace them with the first or last index
    if index_first < 0: 
        index_first = 0
    if index_last > nfiles - 1: 
        index_last = nfiles - 1

    # Return the desired indices
    return index_first, index_last
            
def strip_dirname(dirname):
    dirname_stripped = dirname.split('/')[-1]
    if (dirname_stripped == ''):
        dirname_stripped = dirname.split('/')[-2]
    if (dirname == '.'):
        dirname_stripped = os.getcwd().split('/')[-1]
    if (dirname == '..'):
        orig_dir = os.getcwd()
        os.chdir('..')
        dirname_stripped = os.getcwd().split('/')[-1]
        os.chdir(orig_dir)
    return dirname_stripped

def strip_filename(filename):
    filename_stripped = filename[:-4] # strip off the .npy
    filename_split = filename_stripped.split('_')
    return filename_split[0]

def get_iters_from_file(filename):
    filename_stripped = filename[:-4] # strip off the .npy
    filename_split = filename_stripped.split('_')
    iter1, iter2 = int(filename_split[-2]), int(filename_split[-1])
    return iter1, iter2

def is_an_int(string):
    len_str = len(string)
    bool_val = True
    for i in range(len_str):
        char = string[i]
        bool_val *= (char >= '0' and char <= '9')
    return(bool(bool_val))
        
def get_widest_range_file(datadir, data_name):
    # Find the desired file(s) in the data directory. If there are 
    # multiple, by default choose the one with widest range in the
    # trace/average/distribution
    # If there is no matching file, return the empty string
    datafiles = os.listdir(datadir)
    len_name = len(data_name)
    specific_files = []
    for i in range(len(datafiles)):
        datafile = datafiles[i]
        if data_name in datafile:
            istart = datafile.find(data_name)
            num = 1 # iterations should usually start right after the 
                    # data type name
            if data_name == 'time-longitude': # except for time-latitude
                num = 17
            possible_iter = datafile[istart + len_name + num:istart + len_name + num + 8]
            if is_an_int(possible_iter):
                specific_files.append(datafile)

    ranges = []
    iters1 = []
    iters2 = []
    if len(specific_files) > 0:
        for specific_file in specific_files:
            specific_file_stripped = specific_file[:-4] 
                # get rid of '.npy'...
            li2 = specific_file_stripped.split('_')
            iter1, iter2 = int(li2[-2]), int(li2[-1])
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
        return specific_files[ind]
    else:
        return ''

def frac_nonzero(arr):
    num_nonzero = len(np.where(arr != 0)[0])
    num_total = np.size(arr)
    return (num_nonzero/num_total)

def interpy(x1, y1, x2, y2, x):
    ''' Given two points (x1, y1), (x2, y2) and an abscissa x, interpy 
    interpolates (or extrapolates) to find y(x)'''
    slope = (y2 - y1)/(x2 - x1)
    y = y1 + slope*(x - x1)
    return (y)

def interpx(x1, y1, x2, y2, y):
    ''' Given two points (x1, y1), (x2, y2) and an ordinate y, interpx 
    interpolates (or extrapolates) to find x(y)'''
    slope = (y2 - y1)/(x2 - x1)
    x = x1 + (y - y1)/slope
    return (x)

def get_dict(fname):
    if fname[-3:] == 'npy':
        di = np.load(fname, encoding='latin1', allow_pickle=True).item()
    elif fname[-3:] == 'pkl':
        f = open(fname, 'rb')
        di = pickle.load(f)
        f.close()
    return di

def rms(array):
    if np.size(array) == 0:
        return 0
    else:
        return np.sqrt(np.mean(array**2))

def get_satvals(field, posdef=False, logscale=False, symlog=False):
    # Get good place to saturate array [field], assuming either
    # posdef (True or False) and/or logscale (True or False)
    # and/or symlog (True or False)
    if logscale:
        logfield = np.log(field)
        medlog = np.median(logfield)
        shiftlog = logfield - medlog
        std_plus =\
            np.std(shiftlog[np.where(shiftlog > 0.)].flatten())
        std_minus =\
            np.std(shiftlog[np.where(shiftlog <= 0.)].flatten())
        av_std = (std_plus + std_minus)/2.

        minexp = medlog - 5.*av_std
        maxexp = medlog + 5.*av_std
        minmax = np.exp(minexp), np.exp(maxexp)        
    elif posdef:
        sig = rms(field)
        minmax = 0., 3.*sig        
    elif symlog:
        maxabs = np.max(np.abs(field))
        minmax = -maxabs, maxabs       
    else:
        sig = np.std(field)
        minmax = -3.*sig, 3.*sig
    return minmax

def get_symlog_params(field, field_max=None):
    if field_max is None:
        maxabs = np.max(np.abs(field))
        maxabs_exp = np.floor(np.log10(maxabs))
        field_max = 10.**maxabs_exp
    sig = np.std(field)
    linthresh = 0.15*sig
    dynamic_range = field_max/linthresh
    dynamic_range_decades = np.log10(dynamic_range)
    linscale = dynamic_range_decades
    return linthresh, linscale
    
def saturate_array(arr, my_min, my_max):
    arr[np.where(arr < my_min)] = my_min
    arr[np.where(arr > my_max)] = my_max

def get_exp(num):
    if num != 0.:
        return int(np.floor(np.log10(np.abs(num))))
    else:
        return 1

def sci_format(num):
    exponent = get_exp(num)
    mantissa = num/10.**exponent
    return (r'$%1.1f\times10^{%i}$' %(mantissa, exponent))


def my_bool(x):
    x = x.lower()
    if x == 't' or x == 'true' or x == '1':
        return True
    elif x == 'f' or x == 'false' or x == '0':
        return False

def trim_field(field, rr, cost):
    # "Trim" a field--i.e., return the field with the values within
    # 5% of the boundary by radius, and within 15 degrees of the poles, 
    # all removed.

    # Get rid of the data close to poles
    lats = (np.pi/2. - np.arccos(cost))*180./np.pi
    lat_cutoff = 75.
    it_cutm = np.argmin(np.abs(lats + lat_cutoff))
    it_cutp = np.argmin(np.abs(lats - lat_cutoff))

    # Also stay away from within 5 percent of top and bottom!
    ri, ro = np.min(rr), np.max(rr)
    shell_depth = ro - ri
    rr_depth = (ro - rr)/shell_depth
    ir_cuttop = np.argmin(np.abs(rr_depth - 0.05))
    ir_cutbot = np.argmin(np.abs(rr_depth - 0.95))
    field_cut = field[it_cutm:it_cutp+1, ir_cuttop:ir_cutbot + 1]
    return field_cut

def append_logfile(logfile, message):
    f = open(logfile, 'a')
    f.write(message)
    f.close()