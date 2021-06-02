# Author: Loren Matilsky
# Created: 11/10/2018
# Common routines for routines for post-processing Rayleigh data 

import numpy as np
import sys, os, pickle
from string_to_num import string_to_number_or_array
sys.path.append(os.environ['rapp'])
from reference_tools import equation_coefficients
from rayleigh_diagnostics import G_Avgs, Shell_Slices, ReferenceState,\
    TransportCoeffs, GridInfo
from compute_grid_info import compute_grid_info, compute_theta_grid,\
        compute_r_grid

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

#Thermodyanmic variables
c_P = 3.5e8
thermo_gamma = 5./3.
thermo_R = c_P*(1. - 1./thermo_gamma)

# I am now calling r_m the base of the convection zone, 
# while r_i (the inner shell radius) can vary
rhom = 0.18053428
Tm = 2111256.4
rm = 5.0e10
ro = 6.5860209e10 # Radii consistent with the bottom 3 density scale 
        # heights in the Sun rho_i above corresponds to the density
        # at the base of the convection zone

# width of print messages in parallel routines
lent = 50
buff_frac = 0.05 # default buffer to set axes limits

# character to make text bold
bold_char_begin = "\033[1m"
bold_char_end = "\033[0m"

letters = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l',\
    'm', 'n', 'o', 'p', 'q', 'r', 's', 't', 'u', 'v', 'w', 'x', 'y', 'z']

def make_bold(st):
    return bold_char_begin + st + bold_char_end

def format_time(seconds):
    if seconds < 1: # output in milliseconds
        return "%5i ms" %(round(seconds*1000))
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

def inds_from_vals(arr, arrvals):
    nind = len(arrvals)
    indarr = np.zeros(nind, 'int')
    for i in range(nind):
        indarr[i] = np.argmin(np.abs(arr - arrvals[i]))
    return indarr

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

def get_desired_range(int_file_list, args):
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

def get_file_lists(radatadir, args):
    # Get file names in datadir and their integer counterparts
    # (only the ones in the desired range determined by args)

    # get all files
    file_list, int_file_list, nfiles = get_file_lists_all(radatadir)
    # get the desired range
    index_first, index_last = get_desired_range(int_file_list, args)
    # Remove parts of file lists we don't need
    file_list = file_list[index_first:index_last + 1]
    int_file_list = int_file_list[index_first:index_last + 1]
    nfiles = index_last - index_first + 1
    return file_list, int_file_list, nfiles

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

    if wrap and len(dirname_stripped) > 25:
        # Split dirname_stripped into two lines if it is very long
        dirname_stripped = dirname_stripped[:25] + '\n' +\
                dirname_stripped[25:]
    return dirname_stripped

def is_an_int(string):
    len_str = len(string)
    bool_val = True
    for i in range(len_str):
        char = string[i]
        bool_val *= (char >= '0' and char <= '9')
    return(bool(bool_val))

def get_dataname_from_file(filename):
    just_file = filename.split('/')[-1] #strip the possible full path info
    return just_file.split('-')[0]

def get_iters_from_file(filename):
    filename_end = filename.split('-')[-1][:-4] 
    # (gets the [iter1]_[iter2].pkl and removes the trailing .ext)
    iters_st = filename_end.split('_')
    iter1, iter2 = int(iters_st[0]), int(iters_st[1])
    return iter1, iter2
        
def get_widest_range_file(datadir, dataname):
    # Find the desired file(s) in the data directory. If there are 
    # multiple, by default choose the one with widest range in the
    # trace/average/distribution
    # If there is no matching file, return the empty string
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

def get_satvals(field, posdef=False, logscale=False, symlog=False,\
        fullrange=False):
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
    elif symlog or fullrange:
        maxabs = np.max(np.abs(field))
        minmax = -maxabs, maxabs       
    else:
        sig = np.std(field)
        minmax = -3.*sig, 3.*sig
    # Make sure minmax isn't 0, 0
    tinybit = 1.0e-100
    minmax = minmax[0] - tinybit, minmax[1] + tinybit
    return minmax

def buff_minmax(minval, maxval):
    buff = buff_frac*(maxval - minval)
    return minval - buff, maxval + buff

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

def sci_format(num, ndec=1):
    exponent = get_exp(num)
    mantissa = num/10.**exponent
    return ((r'$%1.' + (r'%i' %ndec) + r'f\times10^{%i}$')\
            %(mantissa, exponent))

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
    
def is_positive(arr):
    how_many_positive_points = np.sum(arr >= 0.)
    how_many_overall_points = np.size(arr)
    # If the array is positive everywhere, the number of positive points
    # will equal the total number of points
    return (how_many_positive_points == how_many_overall_points)

def reverse_dict(di):
    return dict((v, u) for u, v in di.items())

def read_log(fname):
    f = open(fname, 'r')
    lines = f.readlines()
    iters = []
    delta_t = []
    iters_per_sec = []
    for line in lines:
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

    if len(iters_per_sec) > 0:
        di = dict({'iters': np.array(iters), 'delta_t': np.array(delta_t),\
                'iters_per_sec': np.array(iters_per_sec), 'ncpu': ncpu})
    else:
        di = dict({'iters': np.array(iters), 'delta_t': np.array(delta_t),\
                'ncpu': ncpu})

    return di

def print_tuple(tup, format_str, prepend=''):
    whole_str = prepend + '('
    for i in range(len(tup)):
        whole_str += (format_str %tup[i])
        if i < len(tup) - 1:
            whole_str += ', '
    whole_str += ')'
    print(whole_str)

def fill_str(stri, lent, char):
    len_loc = len(stri)
    nfill = lent - len_loc
    return stri + char*nfill

def get_parameter(dirname, parameter):
    f = open(dirname + '/main_input')
    lines = f.readlines()
    n = len(lines)
    try:
        for i in range(n):
            if (parameter in lines[i].lower() and '=' in lines[i] and \
                    lines[i][0] != '!' and not (parameter == 'tacho_r' and\
                    lines[i][:8] == 'tacho_r2')):
                line = lines[i]
        line = line[:] # test if line was assigned
    except:
        if parameter == 'magnetism' or parameter == 'use_extrema':
            return False # if magnetism wasn't specified, it is "False"
        else:
            raise Exception('The parameter ' + parameter + ' was not\n' +\
                            'specified in run: ' + dirname + '. \n' +\
                            'exiting NOW\n')
    
    # Make line lowercase
    line = line.lower()

    # Remove spaces and newline character (at the end of each line)
    line = line.replace(' ', '')
    line = line.replace('\n', '')

    # Remove possible trailing comma from line
    if line[-1] == ',':
        line = line[:-1]
 
    equals_index = line.index('=') # find where the actual number
        # or array starts (should be after the equals sign)
    num_string = line[equals_index + 1:]
    if '!' in num_string: # there was a comment after the equals statement
                        # throw it away!
        excl_index = num_string.index('!')
        num_string = num_string[:excl_index]
    return (string_to_number_or_array(num_string))

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

class eq_human_readable:
    """Rayleigh Universal Equation Coefficients Structure
    ----------------------------------
    self.nr          : number of radial points
    self.radius      : radial coordinates
    self.density     : density
    self.rho         : density
    self.dlnrho      : logarithmic derivative of density
    self.d2lnrho     : d_by_dr of dlnrho
    self.temperature : temperature
    self.T           : temperature
    self.dlnT        : logarithmic derivative of temperature
    self.pressure    : pressure (rho*R*T)
    self.P           : pressure (rho*R*T)
    self.dsdr        : radial entropy gradient
    self.gravity     : gravity 
    self.heating     : volumetric heating (Q) 
    self.nu          : momentum diffusivity (kinematic viscosity)
    self.dlnu        : logarithmic derivative of the viscosity
    self.kappa       : temperature diffusivity (thermometric conductivity)
    self.dlkappa     : logarithmic derivative of the temp. diffusivity
    self.eta :       : magnetic diffusivity 
    self.dlneta      : logarithmic derivative of magnetic diffusivity
    self.lum         : (scalar) stellar luminosity
    """
    def __init__(self, nr):
        self.nr = nr
        self.density = np.zeros(nr)
        self.rho = np.zeros(nr) # same as density
        self.dlnrho = np.zeros(nr)
        self.d2lnrho = np.zeros(nr)
        self.temperature = np.zeros(nr)
        self.T = np.zeros(nr) # same as temperature
        self.dlnT = np.zeros(nr)
        self.pressure = np.zeros(nr)
        self.P = np.zeros(nr)
        self.gravity = np.zeros(nr) 
        self.g = np.zeros(nr) # same as gravity
        self.dsdr = np.zeros(nr)
        self.heating = np.zeros(nr)
        self.Q = np.zeros(nr) # same as heating
        self.nu = np.zeros(nr)
        self.dlnu = np.zeros(nr)
        self.kappa = np.zeros(nr)
        self.dlnkappa = np.zeros(nr)
        self.eta = np.zeros(nr) # these should stay zero 
        self.dlneta = np.zeros(nr) # if magnetism = False
        self.lum = 0.0 # luminosity

def get_eq(dirname, fname='equation_coefficients'): # return an eq_human_readable class associated with
    # [dirname], either using equation_coefficients or 
    # transport/reference files
    if os.path.exists(dirname + '/' + fname):
        # by default, get info from equation_coefficients (if file exists)
        eq = equation_coefficients()
        eq.read(dirname + '/' + fname)
        eq_hr = eq_human_readable(eq.nr)

        eq_hr.radius = eq.radius
        eq_hr.density = eq.functions[0]
        eq_hr.rho = eq_hr.density
        eq_hr.dlnrho = eq.functions[7]
        eq_hr.d2lnrho = eq.functions[8]
        eq_hr.temperature = eq.functions[3]
        eq_hr.T = eq_hr.temperature
        eq_hr.pressure = thermo_R*eq_hr.rho*eq_hr.T
        eq_hr.P = eq_hr.pressure
        eq_hr.dlnT = eq.functions[9]
        eq_hr.gravity = eq.functions[1]/eq_hr.rho*c_P
        eq_hr.g = eq_hr.gravity
        eq_hr.dsdr = eq.functions[13]
        eq_hr.heating = eq.constants[9]*eq.functions[5]
        eq_hr.Q = eq_hr.heating
        eq_hr.nu = eq.constants[4]*eq.functions[2]
        eq_hr.dlnu = eq.functions[10]
        eq_hr.kappa = eq.constants[5]*eq.functions[4]
        eq_hr.dlnkappa = eq.functions[11]
        eq_hr.eta = eq.constants[6]*eq.functions[6] # these are built-in to
        eq_hr.dlneta = eq.functions[12] # equation_coefficients as "zero"
        eq_hr.lum = eq.constants[9]
        # if magnetism = False
    else:
        ref = ReferenceState(dirname + '/reference')
        eq_hr = eq_human_readable(ref.nr)

        eq_hr.radius = ref.radius
        eq_hr.density = ref.density
        eq_hr.rho = eq_hr.density
        eq_hr.dlnrho = ref.dlnrho
        eq_hr.d2lnrho = ref.d2lnrho
        eq_hr.temperature = ref.temperature
        eq_hr.T = eq_hr.temperature
        eq_hr.dlnT = ref.dlnt
        eq_hr.pressure = thermo_R*eq_hr.rho*eq_hr.T
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
        eq_hr.lum = get_parameter(dirname, 'luminosity')
    return eq_hr

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

def compute_Prot(dirname):
    try:
        Om0 = get_parameter(dirname, 'angular_velocity')
    except:
        eq = equation_coefficients()
        eq.read(dirname + '/equation_coefficients')
        Om0 = eq.constants[0]/2.
    return 2*np.pi/Om0

def translate_times(time, dirname, translate_from='iter'):
    # Get the baseline time unit
    rotation = get_parameter(dirname, 'rotation')
    if rotation:
        time_unit = compute_Prot(dirname)
        time_label = 'prot'
    else:
        time_unit = compute_tdt(dirname)
        time_label = 'tdt'

    # First, if translating from an iter, just use an individual data file
    if translate_from == 'iter': # just read in individual G_Avgs file
        file_list, int_file_list, nfiles = get_file_lists_all(dirname + '/G_Avgs')
        if nfiles > 0: # first see if we can use the G_Avgs data
            #print ("translate_times(): translating using G_Avgs data")
            funct = G_Avgs
            radatadir = dirname + '/G_Avgs'
            iiter = np.argmin(np.abs(int_file_list - time))
            a = funct(radatadir + '/' + file_list[iiter], '')
            jiter = np.argmin(np.abs(a.iters - time))
            val_sec = a.time[jiter]
            val_iter = a.iters[jiter]
            val_day = val_sec/86400.
            val_unit = val_sec/time_unit
        else:
            file_list, int_file_list, nfiles = get_file_lists_all(dirname +\
                    '/Shell_Slices')
            if nfiles > 0: # next see if we can use Shell_Slices data
                #print ("translate_times(): translating using Shell_Slices data")
                funct = Shell_Slices
                radatadir = dirname + '/Shell_Slices'
                iiter = np.argmin(np.abs(int_file_list - time))
                a = funct(radatadir + '/' + file_list[iiter], '')
                jiter = np.argmin(np.abs(a.iters - time))
                val_sec = a.time[jiter]
                val_iter = a.iters[jiter]
                val_day = val_sec/86400.
                val_unit = val_sec/time_unit
            else: # Finally use G_Avgs_trace or time-latitude data
                datadir = dirname + '/data/'
                try:        
                    the_file = get_widest_range_file(datadir, 'G_Avgs_trace')
                    di = get_dict(the_file)
                    #print ("translate_times(): translating using G_Avgs_trace file")
                except:
                    the_file = get_widest_range_file(datadir, 'time-latitude')
                    di = get_dict(the_file)
                    print ("translate_times(): translating using time-latitude file")

                # Get times and iters from trace file
                times = di['times']
                iters = di['iters']
                ind = np.argmin(np.abs(iters - time))
                val_sec = times[ind]
                val_iter = iters[ind]
                val_day = times[ind]/86400.
                val_unit = times[ind]/time_unit

    else: # otherwise, hopefully you computed some time traces beforehand!
        # Get the data directory
        datadir = dirname + '/data/'
     
        try:        
            the_file = get_widest_range_file(datadir, 'G_Avgs_trace')
            di = get_dict(the_file)
            #print ("translate_times(): translating using G_Avgs_trace file")
        except:
            the_file = get_widest_range_file(datadir, 'time-latitude')
            di = get_dict(the_file)
            #print ("translate_times(): translating using time-latitude file")

        # Get times and iters from trace file
        times = di['times']
        iters = di['iters']

        if translate_from in ['prot', 'tdt', 'unit']:
            ind = np.argmin(np.abs(times/time_unit - time))
        elif translate_from == 'day':
            ind = np.argmin(np.abs(times/86400. - time))
        elif translate_from == 'sec':
            ind = np.argmin(np.abs(times - time))

        val_sec = times[ind]
        val_iter = iters[ind]
        val_day = times[ind]/86400.
        val_unit = times[ind]/time_unit

    return dict({'val_sec': val_sec,'val_iter': val_iter,\
            'val_day': val_day, 'val_unit': val_unit,\
            'time_unit': time_unit, 'time_label': time_label})

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

def get_domain_bounds(dirname):
    try:
        rmin, rmax = get_parameter(dirname, 'rmin'),\
                get_parameter(dirname, 'rmax')
        nr = get_parameter(dirname, 'n_r')
        domain_bounds = (rmin, rmax)
        ncheby = (nr,)
    except:
        domain_bounds = tuple(get_parameter(dirname, 'domain_bounds'))
        ncheby = tuple(get_parameter(dirname, 'ncheby'))
    return ncheby, domain_bounds

def field_amp(dirname):
    # Make empty dictionary for field-amplitude arrays
    di_out = dict([])

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

    # Read in the Shell_Avgs data
    the_file = get_widest_range_file(datadir, 'Shell_Avgs')
    if not the_file == '':
        print ('field_amp(): Getting velocity (maybe B field) amps from '\
                + the_file)
        di = get_dict(the_file)
        di_out['iter1'], di_out['iter2'] = get_iters_from_file(datadir +\
                the_file)
        vals = di['vals']
        lut = di['lut']
        try:
            # Read in velocity-squared of flows
            # get this from kinetic energy
            vsqr = 2.0*vals[:, 0, lut[402]]/eq.rho
            vsqt = 2.0*vals[:, 0, lut[403]]/eq.rho
            vsqp = 2.0*vals[:, 0, lut[404]]/eq.rho

            vsqr_fluc = 2.0*vals[:, 0, lut[410]]/eq.rho
            vsqt_fluc = 2.0*vals[:, 0, lut[411]]/eq.rho
            vsqp_fluc = 2.0*vals[:, 0, lut[412]]/eq.rho

            vsqr_mean = vsqr - vsqr_fluc
            vsqt_mean = vsqt - vsqt_fluc
            vsqp_mean = vsqp - vsqp_fluc

            di_out['vamp_r'] = np.sqrt(vsqr)
            di_out['vamp_t'] = np.sqrt(vsqt)
            di_out['vamp_p'] = np.sqrt(vsqp)
            di_out['vamp_pol'] = np.sqrt(vsqr + vsqt)
            di_out['vamp_hor'] = np.sqrt(vsqt + vsqp)
            di_out['vamp'] = np.sqrt(vsqr + vsqt + vsqp)

            di_out['vamp_r_fluc'] = np.sqrt(vsqr_fluc)
            di_out['vamp_t_fluc'] = np.sqrt(vsqt_fluc)
            di_out['vamp_p_fluc'] = np.sqrt(vsqp_fluc)
            di_out['vamp_pol_fluc'] = np.sqrt(vsqr_fluc + vsqt_fluc)
            di_out['vamp_hor_fluc'] = np.sqrt(vsqt_fluc + vsqp_fluc)
            di_out['vamp_fluc'] = np.sqrt(vsqr_fluc + vsqt_fluc + vsqp_fluc)

            di_out['vamp_r_mean'] = np.sqrt(vsqr_mean)
            di_out['vamp_t_mean'] = np.sqrt(vsqt_mean)
            di_out['vamp_p_mean'] = np.sqrt(vsqp_mean)
            di_out['vamp_pol_mean'] = np.sqrt(vsqr_mean + vsqt_mean)
            di_out['vamp_hor_mean'] = np.sqrt(vsqt_mean + vsqp_mean)
            di_out['vamp_mean'] = np.sqrt(vsqr_mean + vsqt_mean + vsqp_mean)

        except:
            print ("field_amplitudes(): one or more quantities needed for")
            print("velocity-squared were not output for Shell_Avgs data")
            print("failed to compute the velocity amplitudes (vamps)")

        if magnetism:
            try:
                # Read in B-squared fields
                # get this from magnetic energy
                eightpi = 8.*np.pi
                bsqr = eightpi*vals[:, 0, lut[1102]]
                bsqt = eightpi*vals[:, 0, lut[1103]]
                bsqp = eightpi*vals[:, 0, lut[1104]]

                bsqr_fluc = eightpi*vals[:, 0, lut[1110]]
                bsqt_fluc = eightpi*vals[:, 0, lut[1111]]
                bsqp_fluc = eightpi*vals[:, 0, lut[1112]]

                bsqr_mean = bsqr - bsqr_fluc
                bsqt_mean = bsqt - bsqt_fluc
                bsqp_mean = bsqp - bsqp_fluc

                di_out['bamp_r'] = np.sqrt(bsqr)
                di_out['bamp_t'] = np.sqrt(bsqt)
                di_out['bamp_p'] = np.sqrt(bsqp)
                di_out['bamp_pol'] = np.sqrt(bsqr + bsqt)
                di_out['bamp_hor'] = np.sqrt(bsqt + bsqp)
                di_out['bamp'] = np.sqrt(bsqr + bsqt + bsqp)

                di_out['bamp_r_fluc'] = np.sqrt(bsqr_fluc)
                di_out['bamp_t_fluc'] = np.sqrt(bsqt_fluc)
                di_out['bamp_p_fluc'] = np.sqrt(bsqp_fluc)
                di_out['bamp_pol_fluc'] = np.sqrt(bsqr_fluc + bsqt_fluc)
                di_out['bamp_hor_fluc'] = np.sqrt(bsqt_fluc + bsqp_fluc)
                di_out['bamp_fluc'] = np.sqrt(bsqr_fluc + bsqt_fluc + bsqp_fluc)

                di_out['bamp_r_mean'] = np.sqrt(bsqr_mean)
                di_out['bamp_t_mean'] = np.sqrt(bsqt_mean)
                di_out['bamp_p_mean'] = np.sqrt(bsqp_mean)
                di_out['bamp_pol_mean'] = np.sqrt(bsqr_mean + bsqt_mean)
                di_out['bamp_hor_mean'] = np.sqrt(bsqt_mean + bsqp_mean)
                di_out['bamp_mean'] = np.sqrt(bsqr_mean + bsqt_mean + bsqp_mean)

            except:
                print ("field_amplitudes(): one or more quantities needed for")
                print("B-squared (ME) were not output for Shell_Avgs data")
                print("failed to compute the B-field amplitudes (bamps)")

    # For consistency also compute the shell depth
    di_out['shell_depth'] = np.max(rr) - np.min(rr)

    # Return the dictionary 
    return di_out

def length_scales(dirname):
    # Make empty dictionary for length_scale arrays
    di_out = dict([])

    # See if run is magnetic
    magnetism = get_parameter(dirname, 'magnetism')

    # First get mixing length scale (the one we know will be there)
    eq = get_eq(dirname)
    L_rho = -1./eq.dlnrho
    rr = eq.radius
    nr = len(rr)
    di_out['rr'] = rr
    di_out['nr'] = nr
    di_out['L_rho'] = L_rho

    # Get data directory
    datadir = dirname + '/data/'

    # Get field amplitudes
    di_field_amp = field_amp(dirname)

    # Read in the Shell_Avgs data
    the_file = get_widest_range_file(datadir, 'Shell_Avgs')
    if not the_file == '':
        print ('length_scales(): Getting vorticity from ' + the_file)
        di = get_dict(the_file)
        di_out['iter1'], di_out['iter2'] = get_iters_from_file(the_file)
        vals = di['vals']
        lut = di['lut']
        try:
            # Read in enstrophy of convective flows
            vortsqr = vals[:, 0, lut[317]] 
            vortsqt = vals[:, 0, lut[318]]
            vortsqp = vals[:, 0, lut[319]]
            vortsqh = vortsqt + vortsqp
            enstr = vortsqr + vortsqt + vortsqp
            # Read in velocity-squared of convective flows
            vsq = di_field_amp['vamp_fluc']**2.0
            # Compute length scale and put it in dictionary
            L_omr = (vsq/vortsqr)**0.5
            L_omh = (vsq/vortsqh)**0.5
            L_om = (vsq/enstr)**0.5
            di_out['L_omr'] = L_omr
            di_out['L_omh'] = L_omh
            di_out['L_om'] = L_om
        except:
            print ("length_scales(): one or more quantities needed for enstrophy")
            print("were not output for Shell_Avgs data")
            print("failed to compute L_om")

        if magnetism:
            print ('Getting del x B currents from ' + the_file)
            try:
                # Read in current of convective fields
                del_crossB2 = vals[:, 0, lut[1015]] + vals[:, 0, lut[1018]] + vals[:, 0, lut[1021]]
                # Read in B-squared of convective flows
                B2 = di_field_amp['bamp_fluc']**2.0
                # Compute length scale and put it in dictionary
                L_J = (B2/del_crossB2)**0.5
                di_out['L_J'] = L_J
            except:
                print ("one or more quantities needed for current or")
                print("were not output for Shell_Avgs data")
                print("failed to compute L_J, assigning it L_om ")
                di_out['L_J'] = L_om

    # Read in the Shell_Spectra data
    the_file = get_widest_range_file(datadir, 'Shell_Spectra')
    if not the_file == '':
        print ('length_scales(): Reading Shell_Spectra data from ' +\
                the_file)
        di = get_dict(the_file)
        lpower = di['lpower']
        nell = np.shape(lpower)[0]
        lvals = np.arange(nell)
        lvals = lvals.reshape((nell, 1))
        lut = di['lut']
        rr_spec = di['rvals']
        # get the convective power    
        vrsq_power = lpower[:, :, lut[1], 2]
        vtsq_power = lpower[:, :, lut[2], 2] 
        vpsq_power = lpower[:, :, lut[3], 2] 
        vhsq_power = vtsq_power + vpsq_power
        vsq_power = vrsq_power + vtsq_power + vpsq_power
        # Compute rms l-values
        l_rms_vr = np.sum(vrsq_power*(lvals + 1.), axis=0)/np.sum(vrsq_power, axis=0)
        l_rms_vh = np.sum(vhsq_power*(lvals + 1.), axis=0)/np.sum(vhsq_power, axis=0)
        l_rms_v = np.sum(vsq_power*(lvals + 1.), axis=0)/np.sum(vsq_power, axis=0)
        # Compute lengthscales and add to dictionary
        pir = np.pi*rr_spec
        L_vr = pir/l_rms_vr
        L_vh = pir/l_rms_vh
        L_v = pir/l_rms_v
        di_out['rr_spec'] = rr_spec
        di_out['ir_spec'] = inds_from_vals(rr, di_out['rr_spec'])
        di_out['nr_spec'] = len(rr_spec)
        di_out['L_vr'] = L_vr
        di_out['L_vh'] = L_vh
        di_out['L_v'] = L_v

        if magnetism:
            Brsq_power = lpower[:, :, lut[801], 2]
            Btsq_power = lpower[:, :, lut[802], 2] 
            Bpsq_power = lpower[:, :, lut[803], 2] 
            Bmsq_power = Brsq_power + Btsq_power
            Bsq_power = Brsq_power + Btsq_power + Bpsq_power
            # Compute rms l-values
            l_rms_Bp = np.sum(Bpsq_power*(lvals + 1.), axis=0)/np.sum(Bpsq_power, axis=0)
            l_rms_Bm = np.sum(Bmsq_power*(lvals + 1.), axis=0)/np.sum(Bmsq_power, axis=0)
            l_rms_B = np.sum(Bsq_power*(lvals + 1.), axis=0)/np.sum(Bsq_power, axis=0)

            # Compute lengthscales and add to dictionary
            L_Bp = pir/l_rms_Bp
            L_Bm = pir/l_rms_Bm
            L_B = pir/l_rms_B
            di_out['L_Bp'] = L_Bp
            di_out['L_Bm'] = L_Bm
            di_out['L_B'] = L_v
    # For consistency also compute the shell depth
    di_out['shell_depth'] = np.max(rr) - np.min(rr)

    # Return the dictionary 
    return di_out

def get_grid_info(dirname):
    di_out = dict({})
    gi = GridInfo(dirname + '/grid_info', '')
    # 1D arrays
    di_out['rr'] = gi.radius
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
    return di_out

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
    di_out['Re_fluc'] = di_amp['vamp_fluc']*shell_depth/eq.nu
    di_out['Re_mean'] = di_amp['vamp_mean']*shell_depth/eq.nu

    di_out['Rehrho'] = di_amp['vamp']*hrho/eq.nu
    di_out['Rehrho_fluc'] = di_amp['vamp_fluc']*hrho/eq.nu
    di_out['Rehrho_mean'] = di_amp['vamp_fluc']*hrho/eq.nu

    L_om = di_len['L_om']
    di_out['Revort'] = di_amp['vamp']*L_om/eq.nu
    di_out['Revort_fluc'] = di_amp['vamp_fluc']*L_om/eq.nu
    di_out['Revort_mean'] = di_amp['vamp_mean']*L_om/eq.nu

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
        di_out['Respec_fluc'] = (di_amp['vamp_fluc']/eq.nu)[ir_spec]*L_v
        di_out['Respec_mean'] = (di_amp['vamp_mean']/eq.nu)[ir_spec]*L_v

    if magnetism: # magnetic Reynolds numbers Rm
        L_J = di_len['L_J']
        di_out['Rm'] = di_amp['vamp']*L_J/eq.eta
        di_out['Rm_fluc'] = di_amp['vamp_fluc']*L_J/eq.eta
        di_out['Rm_mean'] = di_amp['vamp_mean']*L_J/eq.eta

        if have_spec:
            L_B = di_len['L_B']
            di_out['Rmspec'] = (di_amp['vamp']/eq.eta)[ir_spec]*L_B
            di_out['Rmspec_fluc'] = (di_amp['vamp_fluc']/eq.eta)[ir_spec]*L_B
            di_out['Rmspec_mean'] = (di_amp['vamp_mean']/eq.eta)[ir_spec]*L_B

    if rotation: # Rossby numbers
        Om0 = 2*np.pi/compute_Prot(dirname)
        di_out['Ro'] = di_amp['vamp']/(2.0*Om0*shell_depth)
        di_out['Ro_fluc'] = di_amp['vamp_fluc']/(2.0*Om0*shell_depth)
        di_out['Ro_mean'] = di_amp['vamp_mean']/(2.0*Om0*shell_depth)

        di_out['Rohrho'] = di_amp['vamp']/(2.0*Om0*hrho)
        di_out['Rohrho_fluc'] = di_amp['vamp_fluc']/(2.0*Om0*hrho)
        di_out['Rohrho_mean'] = di_amp['vamp_mean']/(2.0*Om0*hrho)

        di_out['Rovort'] = di_amp['vamp']/(2.0*Om0*L_om)
        di_out['Rovort_fluc'] = di_amp['vamp_fluc']/(2.0*Om0*L_om)
        di_out['Rovort_mean'] = di_amp['vamp_mean']/(2.0*Om0*L_om)

        if have_spec:
            di_out['Rospec'] = (di_amp['vamp']/eq.eta)[ir_spec]/(2.0*Om0*L_v)
            di_out['Rospec_fluc'] = (di_amp['vamp_fluc']/eq.eta)[ir_spec]/(2.0*Om0*L_v)
            di_out['Rospec_mean'] = (di_amp['vamp_mean']/eq.eta)[ir_spec]/(2.0*Om0*L_v)

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
def thin_data(vals, ntot):
    nx = np.shape(vals)[0]
    nskip = nx//ntot
    if not nskip in [0, 1]: #for ntot < 2*nx, do nothing
        vals_new = vals[::nskip]
    else:
        vals_new = vals
    return vals_new

def array_of_strings(arr):
    li = []
    for ele in arr:
        li.append(str(ele))
    return np.array(li)

def get_time_unit(dirname):
    rotation = get_parameter(dirname, 'rotation')
    if rotation:
        time_unit = compute_Prot(dirname)
        time_label = r'${\rm{P_{rot}}}$'
    else:
        time_unit = compute_tdt(dirname)
        time_label = r'${\rm{TDT}}$'
    return time_unit, time_label, rotation

def get_time_info(dirname, iter1, iter2):
    # Get the time range in sec
    t1 = translate_times(iter1, dirname, translate_from='iter')['val_sec']
    t2 = translate_times(iter2, dirname, translate_from='iter')['val_sec']

    # Get the baseline time unit
    time_unit, time_label, rotation = get_time_unit(dirname)

    # set the averaging-interval label
    if rotation:
        time_string = ('t = %.1f to %.1f ' %(t1/time_unit, t2/time_unit))\
                + time_label + '\n' + (r'$\Delta t = %.1f\ $'\
                %((t2 - t1)/time_unit)) + time_label
    else:
        time_string = ('t = %.3f to %.3f ' %(t1/time_unit, t2/time_unit))\
                + time_label + '\n' + (r'$\Delta t = %.3f\ $'\
                %((t2 - t1)/time_unit)) + time_label
    return time_string

def my_mkdir(dirname):
    if not os.path.isdir(dirname):
        os.makedirs(dirname)
    return dirname

def arr_to_str(a, fmt):
    st = ''
    for ele in a:
        st += (fmt + ' ') %ele
    return '[' + st[:-1] + ']'

class dotdict(dict):
    """dot.notation access to dictionary attributes"""
    __getattr__ = dict.get
    __setattr__ = dict.__setitem__
    __delattr__ = dict.__delitem__

def update_kwargs(kwargs_supplied, kwargs_default):
    kwargs = {**kwargs_default} # start with default kwargs
    for key, val in kwargs_supplied.items():
        if key in kwargs_default: # only update arguments that are
            # specified in the default set
            kwargs[key] = val
        else:
            print ("you specified an invalid keyword arg: ", key)
    return kwargs
