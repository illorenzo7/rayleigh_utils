# Author: Loren Matilsky
# Created: 11/10/2018
# Common routines for routines for post-processing Rayleigh data 

import numpy as np
import sys, os, pickle
from string_to_num import string_to_number_or_array
sys.path.append(os.environ['rapp'])
from reference_tools import equation_coefficients
from rayleigh_diagnostics import G_Avgs, Shell_Slices
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
thermo_R = c_P*thermo_gamma*(1. - 1./thermo_gamma)

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

def get_file_lists(radatadir):
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

range_options = ['-range', '-centerrange', '-leftrange', '-rightrange',\
        '-n', '-f', '-all', '-iter']
n_options = len(range_options)

def get_desired_range(int_file_list, args):
    nargs = len(args)
    nfiles = len(int_file_list)
    # By default, the range will always be the last 100 files:
    index_first, index_last = nfiles - 100, nfiles - 1
    # Determine if user specified any sort of range
    # If not return None instead of (index_first, index_last) tuple
    user_specified_range = False
    for i in range(n_options):
        if range_options[i] in args:
            user_specified_range = True
    if user_specified_range:
        for i in range(nargs):
            arg = args[i]
            if arg in ['-range', '-centerrange', '-leftrange',\
                    '-rightrange', '-iter']: # first arg will be iter no.
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
            if arg in ['-centerrange', '-rightrange', '-leftrange']:
                # many options include an "ndatafiles" argument
                ndatafiles = int(args[i+2])
            elif arg in ['-n', '-f']:
                ndatafiles = int(args[i+1])
            if arg == '-range': # average between two specific files
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
            elif arg == '-centerrange': #range centered around specific file
                if ndatafiles % 2 == 0: #ndatafiles is even
                    index_first = index - ndatafiles//2 + 1
                    index_last = index + ndatafiles//2
                else:  #ndatafiles is odd
                    index_first = index - ndatafiles//2
                    index_last = index + ndatafiles//2
            elif arg == '-leftrange': # range with specific file first
                index_first = index
                index_last = index + ndatafiles - 1
            elif arg == '-rightrange': # range with specific file last
                index_last = index
                index_first = index - ndatafiles + 1
            elif arg == '-n': 
                # range with certain no. files ending with the last
                index_last = nfiles - 1
                index_first = nfiles - ndatafiles
            elif arg == '-f': 
                # range with certain no. files starting with the first
                index_first = 0
                index_last = ndatafiles - 1
            elif arg == '-all': # all files
                index_first = 0
                index_last = nfiles - 1
            elif arg == '-iter': # just get 1 iter
                index_first = index_last = index
        # Check to see if either of the indices fall "out of bounds"
        # and if they do replace them with the first or last index
        if index_first < 0: 
            index_first = 0
        if index_last > nfiles - 1: 
            index_last = nfiles - 1
        the_tuple = index_first, index_last
    else: # user screwed up specifying the range; return nothing
        the_tuple = None

    # Return the desired indices
    return the_tuple
            
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
                num = 16
            possible_iter = datafile[istart + len_name + num:istart + len_name + num + 8]
            if is_an_int(possible_iter):
                if data_name == 'G_Avgs' or data_name == 'Shell_Avgs':
                    # can't confuse "G_Avgs" or "Shell_Avgs"  with 
                    # "trace_G_Avgs"/"trace_Shell_Avgs"; 
                    # please NEVER make a run directory with "trace" in 
                    # the name!
                    if not 'trace' in datafile and \
                            not 'inte_from' in datafile:
                        specific_files.append(datafile)
                else:
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
    # Make sure minmax isn't 0, 0
    tinybit = 1.0e-100
    minmax = minmax[0] - tinybit, minmax[1] + tinybit
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
        if "teration" in line:
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
                'iters_per_sec': np.array(iters_per_sec)})
    else:
        di = dict({'iters': np.array(iters), 'delta_t': np.array(delta_t)})

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
                            'exiting ....\n')
    
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
        file_list, int_file_list, nfiles = get_file_lists(dirname + '/G_Avgs')
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
            file_list, int_file_list, nfiles = get_file_lists(dirname +\
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
            else: # Finally use trace_G_Avgs or time-latitude data
                datadir = dirname + '/data/'
                try:        
                    the_file = get_widest_range_file(datadir, 'trace_G_Avgs')
                    di = get_dict(datadir + the_file)
                    #print ("translate_times(): translating using trace_G_Avgs file")
                except:
                    the_file = get_widest_range_file(datadir, 'time-latitude')
                    di = get_dict(datadir + the_file)
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
            the_file = get_widest_range_file(datadir, 'trace_G_Avgs')
            di = get_dict(datadir + the_file)
            #print ("translate_times(): translating using trace_G_Avgs file")
        except:
            the_file = get_widest_range_file(datadir, 'time-latitude')
            di = get_dict(datadir + the_file)
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

def drad(arr, radius):
    nt, nr = np.shape(arr)
    two_dr = np.zeros((1,nr-2))
    two_dr[0, :] = radius[:nr-2] - radius[2:nr]     
    deriv = np.zeros_like(arr)
    deriv[:,1:nr-1] = (arr[:,:nr-2] - arr[:,2:nr])/two_dr
    deriv[:,0] = deriv[:,1]
    deriv[:,nr-1] = deriv[:,nr-2]
    return deriv

def dth(arr,theta):
    nt, nr = np.shape(arr)
    two_dt = np.zeros((nt-2 ,1))
    two_dt[:, 0] = theta[:nt-2] - theta[2:nt]     
    deriv = np.zeros_like(arr)
    deriv[1:nt-1, :] = (arr[:nt-2,:] - arr[2:nt,:])/two_dt
    deriv[0, :] = deriv[1, :]
    deriv[nt-1, :] = deriv[nt-2, :]
    return deriv

def drad_3d(arr, radius):
    nphi, nt, nr = np.shape(arr)
    two_dr = np.zeros((1, 1, nr-2))
    two_dr[0, 0, :] = radius[:nr-2] - radius[2:nr]     
    deriv = np.zeros_like(arr)
    deriv[:, :, 1:nr-1] = (arr[:, :, :nr-2] - arr[:, :, 2:nr])/two_dr
    deriv[:, :, 0] = deriv[:, :, 1]
    deriv[:, :, nr-1] = deriv[:, :, nr-2]
    return deriv

def dth_3d(arr, theta):
    nphi, nt, nr = np.shape(arr)
    two_dt = np.zeros((1, nt-2, 1))
    two_dt[0, :, 0] = theta[:nt-2] - theta[2:nt]     
    deriv = np.zeros_like(arr)
    deriv[:, 1:nt-1, :] = (arr[:, :nt-2, :] - arr[:, 2:nt, :])/two_dt
    deriv[:, 0, :] = deriv[:, 1, :]
    deriv[:, nt-1, :] = deriv[:, nt-2, :]
    return deriv

def dph_3d(arr):
    nphi, nt, nr = np.shape(arr)
    dphi = 2.*np.pi/nphi
    arr2 = np.roll(arr, -1, axis=0)
    arr1 = np.roll(arr, 1, axis=0)
    deriv = (arr2 - arr1)/2./dphi
    return deriv

def deriv_1d(x, y):
    n = len(x)
    deriv = np.zeros(n)
    
    two_dx = x[2:] - x[:n-2]
    two_dy = y[2:] - y[:n-2]
    deriv[1:n-1] = two_dy/two_dx
    deriv[0] = deriv[1]
    deriv[n-1] = deriv[n-2]
    return(deriv)
    
def int_1d(x, y):
    n = len(x)
    dx = x[1:] - x[:n-1]
    ymid = 0.5*(y[1:] + y[:n-1])
    return (np.sum(ymid*dx))

def int_r(arr, radius, use_extrema=False):
    int_arr = np.zeros_like(arr)
    nt, nr = np.shape(arr)
    rmin, rmax = np.min(radius), np.max(radius)
    r, rw = compute_r_grid(nr, rmin, rmax, use_extrema=use_extrema)
    rw_2d = rw.reshape((1, nr))

    int_arr[:, nr - 1] = 0.  # Integral from rmin to rmin: = 0
    for ir in range(nr - 1):
        r0 = r[ir]
        dr = (1./3.)*(r0**3 - rmin**3)*rw_2d/r**2/np.sum(rw[ir:])
        int_arr[:, ir] = np.sum((arr*dr)[:, ir:], axis=1)
    return int_arr

def int_th(arr):
    int_arr = np.zeros_like(arr)
    nt, nr = np.shape(arr)
    tt, tw = compute_theta_grid(nt)
    cost, sint = np.cos(tt), np.sin(tt)

    int_arr[nt - 1, :] = 0.  # Integral from theta=0 to theta=0: = 0

    for it in range(nt):
        theta0 = tt[it]
        dt = (1. - cost[it])*tw*np.sin(tt)/np.sum(tw[it:])
        dt_2d = dt.reshape((nt, 1))
        int_arr[it, :] = np.sum((arr*dt_2d)[it:, :], axis=0)

    return int_arr

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

def get_length_scales(dirname):
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

    # Read in the Shell_Avgs data
    the_file = get_widest_range_file(datadir, 'Shell_Avgs')
    if not the_file == '':
        print ('length_scales(): Getting velocity/vorticity from ' +\
                the_file)
        di = get_dict(datadir + the_file)
        di_out['iter1'], di_out['iter2'] = get_iters_from_file(datadir +\
                the_file)
        vals = di['vals']
        lut = di['lut']
        try:
            # Read in enstrophy of convective flows
            vortsqr = vals[:, lut[317]] 
            vortsqt = vals[:, lut[318]]
            vortsqp = vals[:, lut[319]]
            vortsqh = vortsqt + vortsqp
            enstr = vortsqr + vortsqt + vortsqp
            # Read in velocity-squared of convective flows
            vsq = vals[:, lut[422]] + vals[:, lut[423]] + vals[:, lut[424]]
            # Compute length scale and put it in dictionary
            L_omr = (vsq/vortsqr)**0.5
            L_omh = (vsq/vortsqh)**0.5
            L_om = (vsq/enstr)**0.5
            di_out['L_omr'] = L_omr
            di_out['L_omh'] = L_omh
            di_out['L_om'] = L_om
        except:
            print ("length_scales(): one or more quantities needed for enstrophy or")
            print("velocity-squared were not output for Shell_Avgs data")
            print("failed to compute L_om")

        if magnetism:
            print ('Getting B fields/currents from ' + the_file)
            try:
                # Read in current of convective fields
                del_crossB2 = vals[:, lut[1015]] + vals[:, lut[1018]] +\
                        vals[:, lut[1021]]
                # Read in B-squared of convective flows
                B2 = 8.*np.pi*(vals[:, lut[1110]] + vals[:, lut[1111]] +\
                        vals[:, lut[1112]])
                # Compute length scale and put it in dictionary
                L_J = (B2/del_crossB2)**0.5
                di_out['L_J'] = L_J
            except:
                print ("one or more quantities needed for current or")
                print("B-squared were not output for Shell_Avgs data")
                print("failed to compute L_J, assigning it L_om ")
                di_out['L_J'] = L_om

    # Read in the Shell_Spectra data
    the_file = get_widest_range_file(datadir, 'Shell_Spectra')
    if not the_file == '':
        print ('length_scales(): Reading Shell_Spectra data from ' +\
                the_file)
        di = get_dict(datadir + the_file)
        lpower = di['lpower']
        rinds = di['rinds']
        nell = di['nell']
        lvals = di['lvals'].reshape((nell, 1))
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
