# Routines to deal with command-line arguments (CLAs)
# Created: 04/17/2021

import sys, os
sys.path.append(os.environ['raco'] + '/quantities_util')
import numpy as np
from common import *
from varprops import *
from lut import *

def read_string_after_arg(args, i):
    # read string with CLA arg (i.e., the string after --arg)
    args_after = args[i+1:]
    nafter = len(args_after)
    iend_found = False
    for j in range(nafter):
        arg = args_after[j]
        if arg[:2] == '--' and not iend_found:
            iend = j
            iend_found = True
    if not iend_found:
        iend = nafter
    return args_after[:iend]

def read_cla_vals(args, i):
    # read values associated with CLA arg (from string after --arg)
    vals_string = read_string_after_arg(args, i)
    vals = []
    if len(vals_string) == 0: # must be a boolean set to True
        vals.append(True)
    for st in vals_string:
        vals.append(string_to_number_or_array(st))
    vals = np.array(vals)

    # if the array has only one value, make it not an array
    if len(vals) == 1:
        vals = vals[0]
    return vals

def read_clas_raw(args): 
    # get all command-line arguments and their values
    # avoid get_parameter stuff
    # and fancy argument designation stuff

    # these are the very basic (0th-order) command line arguments
    clas0 = dotdict()
    clas0['routinename'] = args[0].split('/')[-1][:-3]
    dirname = args[1]
    clas0['dirname'] = dirname
    clas0['tag'] = ''

    # get the other arguments
    clas = dotdict()
    args = args[2:]
    nargs = len(args)
    for i in range(nargs):
        arg = args[i]
        if '--' in arg:
            key = arg[2:]
            clas[key] = read_cla_vals(args, i)

    return dotdict(clas0), dotdict(clas)

def read_clas(args):
    # get all command-line arguments and their values
    # interpret certain special CLAs and values like 
    # --nosave and 
    # --rvals rmin rmax

    # first get basic info
    clas0 = dotdict()
    clas0['routinename'] = args[0].split('/')[-1][:-3]
    dirname = args[1]
    clas0['dirname'] = dirname
    clas0['datadir'] = dirname + '/data/'
    clas0['plotdir'] = dirname + '/plots/'
    clas0['saveplot'] = True
    clas0['showplot'] = True
    clas0['tag'] = ''

    # see if magnetism/rotation are on
    clas0['magnetism'] = get_parameter(dirname, 'magnetism')
    clas0['rotation'] = get_parameter(dirname, 'rotation')

    # get the other arguments
    clas = dotdict()
    args = args[2:]
    nargs = len(args)
    for i in range(nargs):
        # some arguments get treated differently 
        # (they have keyword shortcuts)
        # basic args first
        arg = args[i]
        if arg == '--noshow':
            clas0['showplot'] = False
        elif arg == '--nosave':
            clas0['saveplot'] = False
        elif arg == '--datadir':
            clas0['datadir'] = args[i+1] + '/'
        elif arg == '--plotdir':
            clas0['plotdir'] = args[i+1] + '/'
        elif arg == '--tag':
            clas0['tag'] = '_' + args[i+1]

        # then the rest of the args
        elif arg == '--width':
            clas['fig_width_inches'] = float(args[i+1])
        elif arg == '--subwidth':
            clas['sub_width_inches'] = float(args[i+1])
        elif arg == '--usefile':
            clas['the_file'] = args[i+1]
        elif arg == '--usefileaz':
            clas['the_file_az'] = args[i+1]
        elif arg == '--usefilesh':
            clas['the_file_sh'] = args[i+1]
        elif arg == '--nocontour':
            clas['plotcontours'] = False
        elif arg == '--nobound':
            clas['plotboundary'] = False
        elif arg == '--nolat':
            clas['plotlatlines'] = False

        elif arg == '--rvals':
            if args[i+1] == 'all': # just leave this one alone
                clas.rvals = 'all'
            else:
                string_after = read_string_after_arg(args, i)
                clas['rvals'] = interpret_rvals(clas0.dirname, string_after)
        elif arg == '--rzones':
            zone_heights = np.array([0.05] +\
                    np.linspace(0.125, 0.875, 7).tolist() + [0.95])
            string_after = read_string_after_arg(args, i)
            zone_bounds = np.sort(interpret_rvals(clas0.dirname, string_after))
            rvals = []
            ndom = len(zone_bounds) - 1
            for idom in range(ndom):
                rbot, rtop = zone_bounds[idom:idom+2]
                rvals += (rbot + (rtop - rbot)*zone_heights).tolist()
                if idom < ndom - 1:
                    rvals += [rtop]
            clas['rvals'] = np.array(rvals)
            clas['rzones'] = True # might need to no if this was specified
        elif arg == '--rrange': # specify range of rvals
            vals_after = read_string_after_arg(args, i)
            nrvals = int(vals_after[2])
            rbot, rtop = interpret_rvals(clas0.dirname, vals_after[:2])
            nrvals = int(nrvals)
            clas['rvals'] = np.linspace(rbot, rtop, nrvals)
        elif arg == '--nquadr':
            rmin, rmax = get_rminmax(clas0.dirname)
            nquadr = read_cla_vals(args, i)
            clas['rvals'] = np.linspace(rmin, rmax, nquadr + 1)

            # also want to save nquadr in this case
            key = arg[2:]
            clas[key] = read_cla_vals(args, i)

        elif arg == '--latrange':
            latmin, latmax, nlatvals = read_cla_vals(args, i)
            nlatvals = int(nlatvals)
            clas['latvals'] = np.linspace(latmin, latmax, nlatvals)
        elif arg == '--nquadlat':
            latmin, latmax = get_latminmax(clas0.dirname)
            nquadlat = read_cla_vals(args, i)
            clas['latvals'] = np.linspace(latmin, latmax, nquadlat + 1)

        elif arg == '--v':
            clas.verbose = True

        # everyone else is easy...
        elif '--' in arg:
            key = arg[2:]
            clas[key] = read_cla_vals(args, i)

    return clas0, clas
