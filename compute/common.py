# Author: Loren Matilsky
# Created: 11/10/2018
# Common routines for routines for post-processing Rayleigh data 

import numpy as np
from scipy.interpolate import interp1d
from scipy.integrate import simps
import sys, os, pickle
from string_to_num import string_to_number_or_array
sys.path.append(os.environ['rapp'])

from reference_tools import equation_coefficients
from rayleigh_diagnostics import G_Avgs, AZ_Avgs, Shell_Avgs, Shell_Slices, Shell_Spectra, Meridional_Slices, Equatorial_Slices, GridInfo
from rayleigh_diagnostics_alt import sliceinfo
from grid_util import compute_grid_info, compute_theta_grid

# handy class for making dictionaries "dot-accessible" "key-accessible" and vice versa
class dotdict(dict):
    """dot.notation access to dictionary attributes"""
    __getattr__ = dict.get
    __setattr__ = dict.__setitem__
    __delattr__ = dict.__delitem__

# container holding basic info about the Rayleigh datatypes
di_radtypes = dotdict({\
    'azav': dotdict({'reading_func': AZ_Avgs, 'dataname': 'AZ_Avgs'}),\
    'shav': dotdict({'reading_func': Shell_Avgs, 'dataname': 'Shell_Avgs'}),\
    'gav': dotdict({'reading_func': G_Avgs, 'dataname': 'G_Avgs'}),\
    'spec': dotdict({'reading_func': Shell_Spectra, 'dataname': 'Shell_Spectra'}),\
    'sslice': dotdict({'reading_func': Shell_Slices, 'dataname': 'Shell_Slices'}),\
    'mer': dotdict({'reading_func': Meridional_Slices, 'dataname': 'Meridional_Slices'}),\
    'eq': dotdict({'reading_func': Equatorial_Slices, 'dataname': 'Equatorial_Slices'}) })


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
sun.tmp_bcz = 2111256.4
sun.r_bcz = sun.rbcz = 5.0e10
sun.rbcz_nond = sun.rbcz/sun.r
# radius for the third density scale height (according to model S)
sun.r_nrho3 = 6.5860209e10 
sun.om_rz = 2.70e-6

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

# get a dictionary routine
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

# might as well get model S now
dir_modelS = os.environ['notes'] + '/Model_S/'
di_modelS = dotdict(get_dict(dir_modelS + 'Model_S.pkl'))

######################################
# FORMATTING ROUTINES FOR PRINT OUTPUT
######################################

def get_exp(num):
    if num != 0.:
        return int(np.floor(np.log10(np.abs(num))))
    else:
        return 1

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

# various formats for numbers
flt_fmt = '%1.3e'

lon_fmt = '%05.1f'
lat_fmt = '%+05.1f'

lon_fmt_tex = '%05.1f' + r'$^\circ$'
lat_fmt_tex = '%+05.1f' + r'$^\circ$'

#def lat_format(latval):
#    if latval < 0:
#        hemisphere = 'S'
#    else:
#        hemisphere = 'n'
#    return hemisphere + '%03.1f' %np.abs(latval)

def array_of_strings(arr):
    li = []
    for ele in arr:
        li.append(str(ele))
    return np.array(li)

def arr_to_str(a, fmt, nobra=False):
    st = ''
    for ele in a:
        st += (fmt + ' ') %ele
    if nobra:
        return st[:-1]
    else:
        return '[' + st[:-1] + ']'

def float_or_sci(num, SF=3):
    # determine the format (float or sci) for minimal spacing of a number
    # in scientific notation
    
    # get the exponent in scientific notation
    st_e = '%e' %num
    E = int(st_e.split('e')[1])
    if E >= 0 and not num == 0.: # number is >= 1
        w_E = SF + 3 # width for scientific notation 
        if SF < E + 2:
            w_F = E + 1
        else:
            w_F = SF + 1
    else:
        w_E = SF + 4
        w_F = SF + np.abs(E) + 1

    if w_F < w_E: # 'sci' wins in border cases
        return 'float'
    else:
        return 'sci'

def compactify_float(num, fmt_type, SF=3):
    ''' returns a string with minimal width of the chosen type
    fmt_type in ['float', 'sci'] and number of sig. figs (SF)'''
    if not fmt_type in ['float', 'sci']:
        print ("error, must choose 'float' or 'sci' for fmt_type")
        return 1

    # round the number first
    fmt_e = ('%1.') + ('%i' %(SF-1)) + 'e'
    st_e = fmt_e %num
    num = float(st_e)
    st_e = fmt_e %num

    # get mantissa, exp (as strings)
    mant, exp = st_e.split('e')
    # exponent (as int)
    E = int(exp)

    if fmt_type == 'sci':
        st = fmt_e %num
        # mantissa is ok. Compactify the exponent
        if exp[0] == '+':
            exp = exp[1:]

        # then remove leading zeros
        while exp[0] == '0':
            exp = exp[1:]
        # I feel like I might be left with a "-" for things close to zero. If so, remove the exponent completely
        if exp == '-':
            exp = ''

        # now rebuild the element
        if exp == '':
            st = mant + 'e0'
        else:
            st = mant + 'e' + exp

    if fmt_type == 'float':
        if E >= 0 and not num == 0.: # number is >= 1
            if SF >= E+2: # there is a floating point
                fmt_f = '%.' + ('%i' %(SF-E-1)) + 'f'
            else: # there is no floating point
                fmt_f = '%.0f'
        else: # number is < 1
            fmt_f = '%.' + ('%i' %(np.abs(E) + SF - 1)) + 'f'
        st = fmt_f %num

    return st
   
def arr_to_str_tab(a, SF=3, fmt=None, header='Row'): # ideal for latex tables
    starr = []
    buffst = ' & '
    st = header + buffst

    for ele in a:
        if fmt is None:
            elest = compactify_float(ele, float_or_sci(ele, SF=SF), SF=SF)
        else:
            elest = fmt %ele

        # add the element to the string 
        st += elest + buffst

    # now remove the final buffst
    st = st[:-len(buffst)]

    # and replace it with \\
    st = st + '\\\\'
    return st

# basic array utilities
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
        return np.array(out)

def isall(arg):
    if np.isscalar(arg):
        if arg == 'all':
            return True
        else:
            return False
    else:
        return False

# dictionary utils
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

def reverse_dict(di):
    return dict((v, u) for u, v in di.items())

########################
# BASIC COMPUTE ROUTINES
########################

# Here is the erf-style function I always use and its derivative
def psifunc(r, r0, delta):
    arr = np.zeros_like(r)
    for i in range(len(r)):
        rloc = r[i]
        if rloc <= r0:
            arr[i] = 0.
        elif rloc < r0 + delta and rloc > r0:
            x = (rloc - r0)/delta
            arr[i] = 1.0 - (1.0 - x**2.0)**2.0
        else:
            arr[i] = 1.0
    return arr

def dpsifunc(r, r0, delta): # might need the derivative
    arr = np.zeros_like(r)
    for i in range(len(r)):
        rloc = r[i]
        if rloc <= r0:
            arr[i] = 0.
        elif rloc < r0 + delta and rloc > r0:
            x = (rloc - r0)/delta
            arr[i] = 4.*(1. - x**2.)*x
        else:
            arr[i] = 1.0
    return arr

def inds_from_vals(arr, arrvals):
    arrvals = make_array(arrvals)
    nind = len(arrvals)
    indarr = np.zeros(nind, 'int')
    for i in range(nind):
        indarr[i] = np.argmin(np.abs(arr - arrvals[i]))
    return indarr

def is_an_int(string):
    # obviously, first check if it's actually an int
    if isinstance(string, np.int32) or isinstance(string, np.int64):
        return True
    # otherwise see if it's an int in disguise
    else:
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

def rms(array, axis=None):
    if np.size(array) == 0:
        return 0
    else:
        return np.sqrt(np.mean(array**2, axis=axis))

def minabs(array):
    return np.min(np.abs(array))

def maxabs(array):
    return np.max(np.abs(array))

def close_to_zero(array):
    # see if an array has values too close to zero to divide by
    tol = 1.0e-12
    the_rms = rms(array)
    if the_rms > 1.0e-100: # we can divide by the rms
        if minabs(array)/the_rms < tol:
            return True
        else:
            return False
    else: # rms must be zero...it's an array of zeros!
        return True

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

def curlphi(arr_r, arr_t, rr, tt):
    nt, nr = np.shape(arr_r)
    rr_2d = rr.reshape((1, nr))
    return ( drad(rr_2d*arr_t, rr) - dth(arr_r, tt) )/rr_2d

def indefinite_integral(integrand, x, x0):
    # check on the order of x
    # if x is descending, reverse both arrays for the integration,
    # then reverse them back at the end
    integrand, x = np.copy(integrand), np.copy(x)
    reverse = False
    if x[-1] < x[0]:
        reverse = True
        x = x[::-1]
        integrand = integrand[::-1]

    # basic grid info
    nx = len(x)
    ix0 = np.argmin(np.abs(x - x0))
   
    # compute indefinite integral
    integral = np.zeros(nx)
    for ix in range(nx):
        if ix >= ix0:
            integral[ix] = simps(integrand[ix0:ix+1], x[ix0:ix+1])
        else:
            integral[ix] = -simps(integrand[ix:ix0+1], x[ix:ix0+1])
    if reverse:
        integral = integral[::-1]

    return integral

def definite_integral(integrand, x, x1, x2):
    # check on the order of x
    # if x is descending, reverse both arrays for the integration,
    # then reverse them back at the end
    integrand, x = np.copy(integrand), np.copy(x)
    reverse = False
    if x[-1] < x[0]:
        reverse = True
        x = x[::-1]
        integrand = integrand[::-1]

    # basic grid info
    nx = len(x)
    ix1 = np.argmin(np.abs(x - x1))
    ix2 = np.argmin(np.abs(x - x2))
   
    # compute definite integral
    return simps(integrand[ix1:ix2+1], x[ix1:ix2+1])

def opt_workload(n, nproc):
    # optimally distributes workload (n tasks) over processes (n workers)
    n_per_proc_min = int(np.floor(n/nproc)) # min workload
    n_per_proc_max = int(np.ceil(n/nproc)) # max workload
    # min/max workloads differ by 1
    r = n/nproc - n_per_proc_min # remainder: r sets optimal number of processes
    # to perform max workload
    nproc_max = int(np.floor(nproc*r))
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

def sliding_average(vals, times, delta_t):
    # simple sliding average
    nt = len(vals)
    slider = np.zeros_like(vals)
    intervals = np.zeros(nt)
    total_time = times[-1] - times[0]
    navg = int(nt*delta_t/total_time)
    over2 = navg//2
    navg = over2*2 + 1

    slider[0,...] = np.mean(vals[:over2, ...], axis=0)
    for it in range(1, nt):
        # begin with prior time's average
        slider[it,...] = slider[it-1,...]

        # calculate times to average over
        it1 = it - over2
        it2 = it + over2

        if it1 >= 0: # subtract out the part of the mean we already have
            slider[it, ...] -= vals[it1, ...]/navg
        else:
            it1 = 0 # need this to calculate interval
        if it2 < nt: # add in the part of the mean we don't have
            slider[it, ...] += vals[it2, ...]/navg
        else:
            it2 = nt - 1 # need this to calculate interval

        intervals[it] = times[it2] - times[it1]
    return slider, intervals

def sliding_average_bak(vals, times, delta_t):
    # simple sliding average
    nt = len(vals)
    slider = np.zeros_like(vals)
    intervals = np.zeros(nt)
    for it in range(nt):
        t0 = times[it]
        t1 = t0 - delta_t/2.
        t2 = t0 + delta_t/2.
        it1, it2 = inds_from_vals(times, [t1, t2])

        slider[it] = np.mean(vals[it1:it2+1], axis=0)
        intervals[it] = times[it2] - times[it1]
    return slider, intervals

def sliding_average_inds(vals, nsmooth):
    # simple sliding average
    nt = len(vals)
    slider = np.zeros(nt)
    for it in range(nt):
        it1 = it - nsmooth//2
        it2 = it1 + nsmooth
        if it1 < 0:
            it1 = 0
        if it2 > nt:
            it2 = nt
        slider[it] = np.mean(vals[it1:it2])
    return slider

# Nonlinear Fourier transforms
def my_nfft(times, arr, axis=0, window=False, renorm=False):
    # shift the times to lie in range -1/2, 1/2
    total_time = times[-1] - times[0]
    times_shift = (times - times[0])/total_time - 1/2
    # get equally spaced times
    times_eq = np.linspace(-1/2, 1/2, len(times))
    interpolant = interp1d(times_shift, arr, axis=axis)
    arr_interp = interpolant(times_eq)
    arr_fft = np.fft.fft(arr_interp, axis=axis)
    arr_fft = np.fft.fftshift(arr_fft, axes=axis)

    # apply a Hann window possibly
    if window:
        the_window = np.hanning(len(times))
        the_shape = np.ones(arr_interp.ndim, dtype='int')
        the_shape[axis] = len(times)
        the_window = the_window.reshape(the_shape)
        arr_fft_window = np.fft.fft(arr_interp*the_window, axis=axis)
        arr_fft_window = np.fft.fftshift(arr_fft_window, axes=axis)

        # normalize the windowed FFT to match the regular FFT
        if renorm:
            ratio = np.sqrt(np.sum(np.abs(arr_fft)**2, axis=axis)/\
                    np.sum(np.abs(arr_fft_window)**2, axis=axis) )
            the_shape = list(ratio.shape)
            the_shape.insert(axis, 1)
            ratio = ratio.reshape(the_shape)
        else:
            ratio = 1.
        arr_fft = arr_fft_window*ratio

    # may as well get frequencies here too
    delta_t = np.mean(np.diff(times))
    freq = np.fft.fftfreq(len(times), delta_t)
    freq = np.fft.fftshift(freq)
    # return everything
    return arr_fft, freq

def interp_nd(arr, x_old, x_new, axis=0):
    '''interpolates ND-array arr (of abitrary number of dimensions)
    along axis from old grid (x_old) to new grid (x_new). To ensure
    this works optimally, ensure the values in x_new lie in the range 
    of the values of x_old. Don't speculate!'''

    nx_old, nx_new = len(x_old), len(x_new)
    old_shape = np.shape(arr)
    new_shapeish = (nx_new,) + old_shape[:axis] + old_shape[axis+1:] 
    rest_of_shape = old_shape[:axis] + old_shape[axis+1:] 
    nrest = np.prod(rest_of_shape)
    arr_flattish = np.zeros((nx_old, nrest))
    arr_interp_flattish = np.zeros((nx_new, nrest))
    for ix in range(nx_old):
        arr_flattish[ix, :] = np.take(arr, indices=ix, axis=axis).flatten()
    for irest in range(nrest):
        f = interp1d(x_old, arr_flattish[:, irest])
        arr_interp_flattish[:, irest] = f(x_new)
    arr_interp = np.reshape(arr_interp_flattish, new_shapeish)

    # finally transpose the zeroth axis if necessary
    if axis > 0:
        arr_interp = np.swapaxes(arr_interp, 0, axis)
    return arr_interp

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
range_options = ['range', 'centerrange', 'leftrange', 'rightrange', 'iters', 'n', 'f', 'all']

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

def get_desired_range(int_file_list, clas):
    # Get first and last index (within the int_file_list) associated with the desired range

    # By default, the range will always be the last 100 files:
    nfiles = len(int_file_list)
    index_first, index_last = nfiles - 100, nfiles - 1

    # user can modify this default in a number of ways
    for key, val in clas.items():
        val = make_array(val)
        if key in ['range', 'centerrange', 'leftrange',\
                'rightrange']: # first arg will be iter no.
            # 'first' means first available file.
            # 'last' means last available file.
            desired_iter = val[0]
            if desired_iter == 'first':
                desired_iter = int_file_list[0]
            elif desired_iter == 'last':
                desired_iter = int_file_list[-1]
            else:
                desired_iter = int(float(desired_iter))
            index = np.argmin(np.abs(int_file_list - desired_iter))
        if key in ['centerrange', 'rightrange', 'leftrange']:
            # many options include an "ndatafiles" argument
            ndatafiles = int(val[1])
        if key in ['n', 'f']:
            ndatafiles = int(val[0])
        if key == 'range': # consider range between two specific files
            index_first = index # first arg is first desired iter
            # also need last desired iter
            desired_iter = val[1]
            if desired_iter == 'first':
                desired_iter = int_file_list[0]
            elif desired_iter == 'last':
                desired_iter = int_file_list[-1]
            else:
                desired_iter = int(float(desired_iter))
            index_last = np.argmin(np.abs(int_file_list - desired_iter))
        if key == 'centerrange': #range centered around specific file
            if ndatafiles % 2 == 0: #ndatafiles is even
                index_first = index - ndatafiles//2 + 1
                index_last = index + ndatafiles//2
            else:  #ndatafiles is odd
                index_first = index - ndatafiles//2
                index_last = index + ndatafiles//2
        if key == 'leftrange': # range with specific file first
            index_first = index
            index_last = index + ndatafiles - 1
        if key == 'rightrange': # range with specific file last
            index_last = index
            index_first = index - ndatafiles + 1
        if key == 'n': 
            # range with certain no. files ending with the last
            index_last = nfiles - 1
            index_first = nfiles - ndatafiles
        if key == 'f': 
            # range with certain no. files starting with the first
            index_first = 0
            index_last = ndatafiles - 1
        if key == 'all': # all files
            index_first = 0
            index_last = nfiles - 1

    # Check to see if either of the indices fall "out of bounds"
    # and if they do replace them with the first or last index
    if index_first < 0: 
        index_first = 0
    if index_last > nfiles - 1: 
        index_last = nfiles - 1

    # Return the desired indices
    return index_first, index_last

def get_file_lists(radatadir, clas):
    # Get file names in datadir and their integer counterparts
    # (only the ones in the desired range determined by args)
    # all the "action" occurs in get_desired_range() function below

    # get all files
    file_list, int_file_list, nfiles = get_file_lists_all(radatadir)

    # see if the last key is "iters" (specific iters)
    # if so, we need to do things slightly differently
    range_key = None
    for key, val in clas.items():
        val = make_array(val)
        if key in range_options:
            range_key = key
            range_val = val

    if range_key == 'iters':
        # first convert the range_val to all ints (they are strings rn)
        tmp = np.ones_like(range_val, dtype='int')
        for i in range(len(range_val)):
            if range_val[i] == 'first':
                tmp[i] = int_file_list[0]
            elif range_val[i] == 'last':
                tmp[i] = int_file_list[-1]
            else:
                tmp[i] = int(float(range_val[i]))
        range_val = tmp

        inds = inds_from_vals(int_file_list, range_val)
        file_list = file_list[inds]
        int_file_list = int_file_list[inds]
    else:
        # get the desired range from the range indices
        index_first, index_last = get_desired_range(int_file_list, clas)
        # Remove parts of file lists we don't need
        file_list = file_list[index_first:index_last + 1]
        int_file_list = int_file_list[index_first:index_last + 1]
    # in all cases, check the number of files we end up with 
    nfiles = len(file_list)

    # see if user wants to skip any files or get specific number 
    # (nfiles) in the range
    for key, val in clas.items():
        val = make_array(val)
        if key == 'skip':
            nskip = int(val[0])
            file_list = file_list[::skip]
            int_file_list = int_file_list[::skip]
            nfiles = len(int_file_list)
        if key == 'nfiles':
            ndesiredfiles = int(val[0])
            skip = nfiles//ndesiredfiles
            file_list = file_list[::skip]
            int_file_list = int_file_list[::skip]
            nfiles = len(int_file_list)

    return file_list, int_file_list, nfiles


########################################################################
# ROUTINES FOR CHARACTERIZING/READING RAYLEIGH POST-PROCESSED DATA FILES
########################################################################

def get_widest_range_file(datadir, dataname, stringent=True):
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
            if stringent:
                if dataname == get_dataname_from_file(datafile):
                    specific_files.append(datafile)
            else:
                if dataname in get_dataname_from_file(datafile):
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

def strip_dirname(dirname, dirdepth=None, wrap=False, ncut=None):
    full_dirname = os.path.abspath(dirname)
    splits = full_dirname.split('/')
    nsplit = len(splits)
    dirname_stripped = ''
    if dirdepth is None: # default dirdepth 1
        dirdepth = 1
    loopmax = min(dirdepth, nsplit)
    for i in range(loopmax):
        dirname_stripped = '/' + splits[nsplit-i-1] + dirname_stripped
    dirname_stripped = dirname_stripped[1:] # remove prepended /

    if ncut is None:
        ncut = 20
    if wrap and len(dirname_stripped) > ncut:
        # Split dirname_stripped into multiple lines if it is very long
        tmp = ''
        nlines = int(np.ceil(len(dirname_stripped)/ncut))
        for iline in range(nlines):
            tmp += dirname_stripped[iline*ncut:(iline+1)*ncut] 
            if iline < nlines - 1:
                tmp += '\n'
        dirname_stripped = tmp
    return dirname_stripped

def get_num_lines(st):
    # count number of lines in a string
    nlines = 1
    for char in st:
        if char == '\n':
            nlines += 1
    return nlines

def get_parameter(dirname, parameter):
    # read a parameter from main_input

    # deal only with lower-case
    parameter = parameter.lower()

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

        # these parameters still have a value
        # if they weren't specified
        if parameter in ['magnetism', 'rotation', 'advect_reference_state']:
            return False # if these weren't specified, they are false
        elif parameter in ['nu_type', 'kappa_type', 'eta_type']:
            return 1
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

def get_minmax(arr):
    newarr = make_array(arr)
    return (np.min(newarr), np.max(newarr))

def get_rminmax(dirname):
    ncheby, domain_bounds = get_domain_bounds(dirname)
    return domain_bounds[0], domain_bounds[-1]

def get_latminmax(dirname):
    tt_lat = get_grid_info(dirname).tt_lat
    return np.min(tt_lat), np.max(tt_lat)

#################################################################################
# Get basic radial coefficients (grid info, reference state) associated with sim.
#################################################################################

####################################
# Routines associated with grid info
####################################
def get_grid_info(dirname, verbose=False, fname=None, ntheta=None):
    # get basic grid info; try to read from grid_info file
    # or directly from main_input if grid_info doesn't exist
    di = dotdict()
    if fname is None: # default
        fname = 'grid_info'

    # get basic grid (colocation points and weights)
    if os.path.exists(dirname + '/' + fname):
        if verbose:
            print ("get_grid_info(): reading grid from grid_info")
            if not ntheta is None:
                print ("except setting nt = %i manually" %ntheta)
        gi = GridInfo(dirname + '/' + fname, '')
        # 1D arrays
        di.rr = rr = gi.radius
        di.rw = rr = gi.rweights
        if ntheta is None:
            di.tt = tt = gi.theta
            di.tw = tw = gi.tweights
        else:
            tt, tw = compute_theta_grid(ntheta)
            di.tt, di.tw = tt, tw
    else:
        if verbose:
            print ("get_grid_info(): inferring grid from main_input")
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

def compute_axial_H(rr, sint, rmin=None, rmax=None):
    # compute the axial distance H between the spheres r_min and r_max
    # assumes rr and sint are two arrays of compatible shape
    xx = rr*sint
    if rmin is None:
        rmin = np.min(rr)
    if rmax is None:
        rmax = np.max(rr)

    d = rmax - rmin
    xx_flat = xx.flatten()
    H_flat = np.zeros_like(xx_flat)
    for ix in range(len(xx_flat)):
        xx_loc = xx_flat[ix]
        if xx_loc > rmin: # outside tangent cylinder
            H_flat[ix] = 2*np.sqrt(rmax**2 - xx_loc**2)
        else:
            H_flat[ix] = np.sqrt(rmax**2 - xx_loc**2) - np.sqrt(rmin**2 - xx_loc**2)
    H = H_flat.reshape(np.shape(xx))
    return H

def interpret_rvals(dirname, rvals):
    # interpret array of rvals (array of strings), some could have the special keywords, rmin, rmid, rmax
    # but otherwise assumed to be float

    # check if None first
    if rvals is None:
        return None

    # get grid
    gi = get_grid_info(dirname)
    rr = gi.rr

    # get the actual min/max/mid of the full domain:
    ncheby, domain_bounds = get_domain_bounds(dirname)
    rmin, rmax = np.min(domain_bounds), np.max(domain_bounds)
    rmid = 0.5*(rmin + rmax)

    # replace these with their closest actual values in rr
    rmin, rmid, rmax = rr[inds_from_vals(rr, [rmin, rmid, rmax])]

    # compute desired rvals
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

def get_sliceinfo(dirname, dataname='Shell_Slices', fname=None):
    radatadir = dirname + '/' + dataname + '/'
    file_list, int_file_list, nfiles = get_file_lists_all(radatadir)
    if fname is None:
        fname = file_list[0]
    
    di = dotdict()
    if dataname in ['Shell_Slices', 'Shell_Spectra']:
        a = sliceinfo(fname, path=radatadir)
        di.nsamplevals = a.nr
        di.samplevals = a.radius
        di.isamplevals = a.inds
    else:
        if dataname == 'Meridional_Slices':
            a = Meridional_Slices(radatadir + fname, '')
            di.nsamplevals = a.nphi
            di.samplevals = 180*a.phi/np.pi
            di.isamplevals = a.phi_inds
        elif dataname == 'Equatorial_Slices':
            a = Equatorial_Slices(radatadir + fname, '')
            # no sampling locations
    di.qv = a.qv
    return di

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

def indefinite_radial_integral(dirname, arr, r0='rmin'):
    gi = get_grid_info(dirname)
    rr, rw = gi.rr, gi.rw
    
    r0 = interpret_rvals(dirname, r0)[0]
    ir0 = inds_from_vals(gi.rr, r0)[0]

    rmin, rmax = get_rminmax(dirname)
    prefactor = (rmax**3. - rmin**3.)/3.
    
    integral = np.zeros_like(rr)
    for ir in range(len(rr)):
        rval = rr[ir]
        if rval >= r0:
            ir1 = ir
            ir2 = ir0
            the_sign = 1.
        else:
            ir1 = ir0
            ir2 = ir
            the_sign = -1.
        integral[ir] = the_sign*prefactor*np.sum((arr*rw/rr**2.)[ir1:ir2+1])
    return integral

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

def get_eq(dirname, fname=None, verbose=False): 
    # return a human readable version of equation_coefficients
    # for [dirname], using either (in order of priority)

    # 1. binary file specified by fname
    # 2. equation_coefficients file
    # 3. custom_reference_binary file
    # 4. polytrope + transport coefs defined by main_input

    # human readable equation coefficients object to store reference state
    eq_hr = dotdict()

    # see if magnetism is on
    eq_hr.magnetism = get_parameter(dirname, 'magnetism')

    # figure out type of reference
    eq_hr.reference_type = get_parameter(dirname, 'reference_type')
    eq_hr.heating_type = get_parameter(dirname, 'heating_type')
    eq_hr.nu_type = get_parameter(dirname, 'nu_type')
    eq_hr.kappa_type = get_parameter(dirname, 'kappa_type')
    if eq_hr.magnetism:
        eq_hr.eta_type = get_parameter(dirname, 'eta_type')

    # read reference state from binary file or main_input
    if fname is None:
        if os.path.exists(dirname + '/' + 'equation_coefficients'):
            fname = 'equation_coefficients'
        elif os.path.exists(dirname + '/' + 'custom_reference_binary'):
            fname = 'custom_reference_binary'

    # by default, get info from equation_coefficients (if file exists)
    if not fname is None:
        if verbose:
            print ("get_eq(): reading reference state from file:", fname)
        # read in equation_coefficients
        eq = equation_coefficients()
        eq.read(dirname + '/' + fname)

        # keep the original functions/constants as "backup"
        eq_hr.constants = eq.constants
        eq_hr.functions = eq.functions

        # radius
        eq_hr.rr = eq.radius

        # do "universal" thermal variables (reference_type independent) first
        eq_hr.rho = eq.functions[0]
        eq_hr.dlnrho = eq.functions[7]
        eq_hr.d2lnrho = eq.functions[8]
        eq_hr.tmp = eq.functions[3]
        eq_hr.dlntmp = eq.functions[9]

        # N^2, gravity, pressure, and transport coefficients are semi-universal
        # the transport coefficients should be fully, but the documentation is not totally transparent with the diffusivities being scaled by val_top
        # this should change post-05/30/2024
        eq_hr.grav = eq.functions[1]/eq_hr.rho 
        eq_hr.nsq = eq_hr.grav*eq.functions[13] # actually this definition is universal 
                # for dim. anelastic, grav now = g/c_p and f_13 = dS/dr
        eq_hr.prs = eq_hr.rho*eq_hr.tmp
        eq_hr.nu = eq.functions[2]
        eq_hr.dlnu = eq.functions[10]
        eq_hr.kappa = eq.functions[4]
        eq_hr.dlnkappa = eq.functions[11]
        if eq_hr.magnetism:
            eq_hr.eta = eq.functions[6] # these are built-in to
            eq_hr.dlneta = eq.functions[12] # equation_coefficients as "zero"

        # rotation rate and period 
        # (in appropriate units based on chosen timescale)
        # this definition is universal
        eq_hr.om0 = eq.constants[0]/2.
        eq_hr.prot = 2.*np.pi/eq_hr.om0

        # volume-averaged angular momentum of shell 
        # (in appropriate units based on chosen timescale)
        gi = get_grid_info(dirname)
        amom_dens = (8*np.pi*eq_hr.om0/3)*eq_hr.rho*eq_hr.rr**4 # do the latitudinal integral analytically
        if not fname == 'equation_coefficients': # no proper weights to work with, just do simpson
            rmin, rmax = get_rminmax(dirname)
            vol = 4*np.pi/3*(rmax**3 - rmin**3)
            eq_hr.amom = definite_integral(amom_dens, eq_hr.rr, rmin, rmax)/vol
        else:
            eq_hr.amom = np.sum(amom_dens/eq_hr.rr**2*gi.rw)

        # some quantities interpretation depends on the reference_type
        # deal with that here
        if eq_hr.reference_type == 2:
            # assume gas is ideal and scale pressure and gravity
            eq_hr.c_p = get_parameter(dirname, 'pressure_specific_heat')
            eq_hr.gas_constant = (gamma_ideal-1)*eq_hr.c_p/gamma_ideal
            eq_hr.prs *= eq_hr.gas_constant
            eq_hr.grav *= eq_hr.c_p

            # scale transport coefficients
            eq_hr.nu *= eq.constants[4]
            eq_hr.kappa *= eq.constants[5]
            if eq_hr.magnetism:
                eq_hr.eta *= eq.constants[6]
            
            # heating
            eq_hr.heat = eq.constants[9]*eq.functions[5] 
            eq_hr.lum = eq.constants[9]

            # this quantity only exists for reference_type = 2 (along with c_p and gas_const
            eq_hr.dsdr = eq.constants[10]*eq.functions[13] 

            # thermal and viscous diffusion time 
            kappa_volav = volav_in_radius(dirname, eq_hr.kappa)
            eq_hr.tdt = (eq_hr.rr[0] - eq_hr.rr[-1])**2/kappa_volav

            nu_volav = volav_in_radius(dirname, eq_hr.nu)
            eq_hr.vdt = (eq_hr.rr[0] - eq_hr.rr[-1])**2/nu_volav

        else: # this assumes a general nondimensionalization
            # for ref type = 1, 3, 4, 5
            rmin, rmax = get_rminmax(dirname)
            vol = 4*np.pi/3*(rmax**3 - rmin**3)
            eq_hr.heat = eq.functions[5]
            integrand = eq_hr.heat*eq_hr.constants[9]/eq_hr.constants[7]
            if not fname == 'equation_coefficients': # no proper weights to work with, just do simpson
                eq_hr.lum = definite_integral(integrand, eq_hr.rr, rmin, rmax)
            else:
                eq_hr.lum = vol * np.sum(integrand/eq_hr.rr**2*gi.rw)

            # thermal diffusion time 
            eq_hr.tdt = 1./eq_hr.constants[5]

            # viscous diffusion time 
            eq_hr.vdt = 1./eq_hr.constants[4]

            # not sure about dsdr yet

    else: # no binary file; get everything from main_input
        # currently only have figured this out for dimensional anelastic
        if eq_hr.reference_type == 2:
            if verbose:
                print ("get_eq(): inferring reference state from main_input")

            # grid info
            gi = get_grid_info(dirname)
            eq_hr.rr = gi.rr 
            zero = np.zeros_like(eq_hr.rr)

            # polytrope parameters
            poly_nrho = get_parameter(dirname, 'poly_nrho')
            poly_n = get_parameter(dirname, 'poly_n')
            poly_rho_i = get_parameter(dirname, 'poly_rho_i')
            poly_mass = get_parameter(dirname, 'poly_mass')
            eq_hr.c_p = get_parameter(dirname, 'pressure_specific_heat')
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
            #eq_hr.nsq = poly.dsdr*eq_hr.grav/eq_hr.c_p
            #eq_hr.gamma = (poly_n+1)/poly_n
            eq_hr.gas_constant = eq_hr.c_p/(poly_n + 1)

            # get heating
            eq_hr.lum = get_parameter(dirname, 'luminosity')
            if eq_hr.heating_type is None: # default heating_type = 0
                eq_hr.heat = zero
            elif eq_hr.heating_type == 1: # heating proportional to pressure
                                    # "consant entropy heating"
                eq_hr.heat = eq_hr.prs - eq_hr.prs[0]
                integral = volav_in_radius(dirname, eq_hr.heat)
                integral *= get_vol(dirname)
                eq_hr.heat = eq_hr.heat*eq_hr.lum/integral
            elif eq_hr.heating_type == 4: # constant energy heating
                eq_hr.heat = zero + eq_hr.lum/get_vol(dirname)

            # get the transport coefficients
            nu_top = get_parameter(dirname, 'nu_top')
            kappa_top = get_parameter(dirname, 'kappa_top')

            if eq_hr.nu_type == 1:
                eq_hr.nu = zero + nu_top
                eq_hr.dlnu = zero
            elif eq_hr.nu_type == 2:
                nu_power = get_parameter(dirname, 'nu_power')
                eq_hr.nu = nu_top*(eq_hr.rho/eq_hr.rho[0])**nu_power
                eq_hr.dlnu = nu_power*eq_hr.dlnrho

            if eq_hr.kappa_type == 1:
                eq_hr.kappa = zero + kappa_top
                eq_hr.dlnkappa = zero
            elif eq_hr.kappa_type == 2:
                kappa_power = get_parameter(dirname, 'kappa_power')
                eq_hr.kappa = kappa_top*(eq_hr.rho/eq_hr.rho[0])**kappa_power
                eq_hr.dlnkappa = kappa_power*eq_hr.dlnrho
            
            if eq_hr.magnetism:
                eta_top = get_parameter(dirname, 'eta_top')
                if eq_hr.eta_type == 1:
                    eq_hr.eta = zero + eta_top
                    eq_hr.dlneta = zero
                elif eq_hr.eta_type == 2:
                    eta_power = get_parameter(dirname, 'eta_power')
                    eq_hr.eta = eta_top*(eq_hr.rho/eq_hr.rho[0])**eta_power
                    eq_hr.dlneta = eta_power*eq_hr.dlnrho

            # finally, get the rotation rate
            eq_hr.om0 = get_parameter(dirname, 'angular_velocity')
            eq_hr.prot = 2*np.pi/eq_hr.om0

    return eq_hr

def get_units(dirname, rvals=['rmin', 'rmax']):
    eq = get_eq(dirname)
    gi = get_grid_info(dirname)
    nsq = eq.dsdr*eq.grav

    # get background volavgs
    r1, r2 = interpret_rvals(dirname, rvals)
    nu_volav = volav_in_radius(dirname, eq.nu, r1, r2)
    rho_volav = volav_in_radius(dirname, eq.rho, r1, r2)
    tmp_volav = volav_in_radius(dirname, eq.tmp, r1, r2)
    grav_volav = volav_in_radius(dirname, eq.grav, r1, r2)
    kappa_volav  = volav_in_radius(dirname, eq.kappa, r1, r2)

    # output units
    di = dotdict()

    # basic parameters (length, time, rho, S)
    H = r2 - r1
    di.l = H
    di.t = 1/(2*eq.om0)
    vol = get_vol(dirname) # make sure to use the full volume
                # to calculate the non-radiative heat flux vs radius
    lstar = vol*np.sum(eq.heat*gi.rw)
    flux_rad = vol/(4*np.pi*eq.rr**2)*np.cumsum(eq.heat*gi.rw)
    flux_nonrad = lstar/(4*np.pi*eq.rr**2) - flux_rad
    flux_volav = volav_in_radius(dirname, flux_nonrad, 'rmid', 'rmax')
    di.s = flux_volav*H/(kappa_volav*rho_volav*tmp_volav)

    # derivative parameters (velocity, b field, etc.)
    di.rho = rho_volav
    di.v = di.l/di.t
    di.b = np.sqrt(4*np.pi*di.rho)*di.v
    di.e = di.rho*di.v**2

    # other background volav
    di.tmp = tmp_volav
    di.grav = grav_volav
    di.kappa = kappa_volav
    di.nu = nu_volav
    di.fluxnr = flux_volav
    di.lum = lstar
    di.vol = vol

    return di

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

def get_time_unit(dirname, tdt=False):
    # get basic time unit of simulation (rotation period or diffusion time)
    rotation = get_parameter(dirname, 'rotation')
    eq = get_eq(dirname)
    if rotation and not tdt:
        time_unit = eq.prot
        time_label = r'${\rm{P_{rot}}}$'
        simple_label = 'rotations'
    else:
        time_unit = eq.tdt
        time_label = simple_label = 'TDT'
    return time_unit, time_label, rotation, simple_label

def translate_times(time, dirname, translate_from='iter', verbose=False):
    # change between different time units (can translate from: 
    # iter, sim time (or measure sim time in TDT or P_rot)
    # TO USE MUST HAVE G_Avgs_trace file or equivalent 
    # (time-lat, time-rad, etc.)

    # Get a time trace data file for the translation
    datadir = dirname + '/data/'
    # first try the G_Avgs_trace...
    the_file = get_widest_range_file(datadir, 'G_Avgs_trace', stringent=False)
    if the_file is None: # next, try time-radius
        the_file = get_widest_range_file(datadir + '/timelat', '', stringent=False)
    if the_file is None: # finally, try time-radius
        the_file = get_widest_range_file(datadir + '/timerad', '', stringent=False)

    # finally, translate the times or else exit with error
    if the_file is None:
        print ("translate_times(): you need to have a trace file")
        print ("to use me! Exiting.")
        sys.exit()
    else:
        if verbose:
            print ("translate_times(): translating from " + the_file)
        di = get_dict(the_file)

    # Get times and iters from trace file
    times = di['times']
    iters = di['iters']

    # get the equation_coefficients file
    eq = get_eq(dirname)

    # translate the time
    if translate_from == 'iter':
        ind = np.argmin(np.abs(iters - time))
    elif translate_from == 'simt':
        ind = np.argmin(np.abs(times - time))
    elif translate_from == 'prot':
        ind = np.argmin(np.abs(times/eq.prot - time))
    elif translate_from == 'tdt':
        ind = np.argmin(np.abs(times/eq.tdt - time))

    # prepare the dictionary to return
    di = dotdict()
    di.val_simt = times[ind]
    di.val_iter = iters[ind]
    di.val_tdt = times[ind]/eq.tdt
    rotation = get_parameter(dirname, 'rotation')
    if rotation:
        di.val_prot = times[ind]/eq.prot
    return di

def get_time_string(dirname, iter1, iter2=None, oneline=False, threelines=False, iter0=None, floatwidth=None, floatprec=None):
    # see if user wants to subtract off base time
    if not iter0 is None:
        t0 = translate_times(iter0, dirname, translate_from='iter')['val_simt']
    else:
        t0 = None

    # Get the time range in sec
    t1 = translate_times(iter1, dirname, translate_from='iter')['val_simt']
    if not t0 is None:
        t1 -= t0
    if not iter2 == None:
        t2 = translate_times(iter2, dirname, translate_from='iter')['val_simt']
        if not t0 is None:
            t2 -= t0

    # Get the baseline time unit
    time_unit, time_label, rotation, simple_label = get_time_unit(dirname)

    # set the averaging-interval label
    if floatprec is None:
        if rotation:
            floatprec = 2
        else:
            floatprec = 4

    if floatwidth is None:
        fmt = '%.' + str(floatprec) + 'f' # measure rotations
    else:
        fmt = '%0' + str(floatwidth) + '.' + str(floatprec) + 'f'

    if not iter2 is None:
        if oneline:
            time_string = (('t = ' + fmt + ' to ' + fmt) %(t1/time_unit, t2/time_unit)) + ' ' + time_label + ' ' + r'$(\Delta t$' + ' = ' + (fmt %((t2 - t1)/time_unit)) + ' ' + time_label + ')'
        elif threelines:
            time_string =\
        't = ' + (fmt %(t1/time_unit)) + ' ' + time_label + ' to\n' + \
        't = ' + (fmt %(t2/time_unit)) + ' ' + time_label  + '\n' +\
        r'$(\Delta t$' + ' = ' + (fmt %((t2 - t1)/time_unit)) + ' ' + time_label + ')'
        else:
            time_string = (('t = ' + fmt + ' to ' + fmt) %(t1/time_unit, t2/time_unit)) + ' ' + time_label +\
                    '\n' + r'$\Delta t$' + ' = ' + (fmt %((t2 - t1)/time_unit)) + ' ' + time_label
    else:
        time_string = (('t = ' + fmt + ' ') %(t1/time_unit)) + time_label

    return time_string

##############################################################
# ROUTINES FOR FIELD AMPLITUDES, LENGTH SCALES AND NON-D NUMBERS
##############################################################

def field_amp(dirname, the_file=None, verbose=False):
    # Make empty dictionary for field-amplitude arrays
    di_out = dotdict()

    # See if run is magnetic
    magnetism = get_parameter(dirname, 'magnetism')

    # First get density
    eq = get_eq(dirname)
    rho = eq.rho
    rr = eq.rr
    nr = eq.nr
    di_out['rr'] = rr

    # Get data directory
    datadir = dirname + '/data/'

    # Read in the Shell_Avgs data
    if the_file is None: # default
        the_file = get_widest_range_file(datadir, 'Shell_Avgs')
    if verbose:
        print ('field_amp(): reading ' + the_file)
    di = get_dict(the_file)
    di_out['iter1'], di_out['iter2'] = get_iters_from_file(the_file)
    vals = di['vals']
    lut = di['lut']
    qv = di['qv']

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
    omsqr = vals[:, 1, lut[301]] 
    omsqt = vals[:, 1, lut[302]]
    omsqp = vals[:, 1, lut[303]]

    if 317 in qv:
        omsqr_fluc = vals[:, 0, lut[317]] 
        omsqt_fluc = vals[:, 0, lut[318]]
        omsqp_fluc = vals[:, 0, lut[319]]
    else: # just set it to zero for now
        omsqr_fluc = np.zeros(nr)
        omsqt_fluc = np.zeros(nr)
        omsqp_fluc = np.zeros(nr)

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

def length_scales(dirname, the_file=None, verbose=False):
    # Make empty dictionary for length_scale arrays
    di_out = dotdict(dict([]))

    # See if run is magnetic
    magnetism = get_parameter(dirname, 'magnetism')

    # First get mixing length scale + grid
    eq = get_eq(dirname)
    if not close_to_zero(eq.dlnrho):
        hrho = -1./eq.dlnrho
        di_out.hrho = hrho

    di_out['rr'] = rr = eq.rr

    # Get data directory
    datadir = dirname + '/data/'

    # Get field amplitudes
    if the_file is None: # default
        the_file = get_widest_range_file(datadir, 'Shell_Avgs')
    di_out['iter1'], di_out['iter2'] = get_iters_from_file(the_file)
    if verbose:
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

def get_term(dirname, vals, lut, quantity, verbose=False):
    if is_an_int(quantity):
        quantity = int(float(quantity))
        if lut[quantity] < 4000: # it's easy
            return vals[..., lut[quantity]]
        # for some quantities, we can do contingencies
        elif quantity == 1404:
            adv_tot = vals[..., lut[1401]]
            adv_fluc = vals[..., lut[1402]]
            if verbose:
                print ("get_term(): getting 1404 from 1401 - 1402")
            return adv_tot - adv_fluc
        elif quantity == 1479:
            vr = vals[..., lut[1]]
            eq = get_eq(dirname)
            c11 = eq.constants[10]
            f14 = eq.functions[13]
            if verbose:
                print ("get_term(): getting 1479 from rho * T * dsdr * vr")
            return eq.rho*eq.tmp*c11*f14*vr
