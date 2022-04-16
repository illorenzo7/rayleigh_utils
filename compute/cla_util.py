# Routines to deal with command-line arguments (CLAs)
# It's a long one!
# Created: 04/17/2021

import numpy as np
#from common import get_parameter, get_domain_bounds, array_of_strings, is_an_int, rsun
from common import *
from varprops import *
from lut import *

def read_cla_vals(args, i):
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
    vals_string = args_after[:iend]

    vals = []
    if len(vals_string) == 0: # must be a boolean set to True
        vals.append(True)
    for st_full in vals_string:
        for st in st_full.split():
            vals.append(string_to_number_or_array(st))
    vals = np.array(vals)

    # if the array has only one value, make it not an array
    if len(vals) == 1:
        vals = vals[0]
    return vals

def read_clas_raw(args): # avoid get_parameter stuff
    # and fancy argument designation stuff
    clas0 = dict({})
    clas0['routinename'] = args[0].split('/')[-1][:-3]
    dirname = args[1]
    clas0['dirname'] = dirname
    clas0['tag'] = ''

    # get the other arguments
    clas = dict({})
    args = args[2:]
    nargs = len(args)
    for i in range(nargs):
        arg = args[i]
        if '--' in arg:
            key = arg[2:]
            clas[key] = read_cla_vals(args, i)

    return dotdict(clas0), dotdict(clas)

def read_clas(args):
    # first get basic info
    clas0 = dict({})
    clas0['routinename'] = args[0].split('/')[-1][:-3]
    dirname = args[1]
    clas0['dirname'] = dirname
    clas0['datadir'] = dirname + '/data/'
    clas0['plotdir'] = dirname + '/plots/'
    clas0['saveplot'] = True
    clas0['showplot'] = True
    clas0['tag'] = ''

    # see if magnetism/rotation are on
    magnetism = get_parameter(dirname, 'magnetism')
    clas0['magnetism'] = magnetism
    rotation = get_parameter(dirname, 'rotation')
    clas0['rotation'] = rotation

    # get the spherical domain
    ncheby, domain_bounds = get_domain_bounds(dirname)
    rMIN = domain_bounds[0]/rsun
    rMAX = domain_bounds[-1]/rsun
    shell_depth = rMAX - rMIN
    if len(ncheby) == 2:
        rBCZ = domain_bounds[1]/rsun
        shell_depth_cz = rMAX - rBCZ
        shell_depth_rz = rBCZ - rMIN

    # get the other arguments
    clas = dict({})
    args = args[2:]
    nargs = len(args)
    for i in range(nargs):
        # some arguments get treated differently
        # basic args first
        arg = args[i]
        if arg == '--noshow':
            clas0['showplot'] = False
        elif arg == '--nosave':
            clas0['saveplot'] = False
        # then the rest of the args
        elif arg == '--width':
            clas['fig_width_inches'] = float(args[i+1])
        elif arg == '--subwidth':
            clas['sub_width_inches'] = float(args[i+1])
        elif arg == '--usefile':
            clas['the_file'] = args[i+1]
        elif arg == '--nocontour':
            clas['plotcontours'] = False
        elif arg == '--nobound':
            clas['plotboundary'] = False
        elif arg == '--nolat':
            clas['plotlatlines'] = False
        elif arg == '--tag':
            clas0['tag'] = '_' + args[i+1]

        elif arg == '--depths':
            clas['rvals'] = rMAX - read_cla_vals(args, i)*d
        elif arg == '--depthscz':
            clas['rvals'] = rMAX - read_cla_vals(args, i)*shell_depth_cz
        elif arg == '--depthirz':
            clas['rvals'] = rBCZ - read_cla_vals(args, i)*shell_depth_rz
        elif arg == '--rvals':
            if args[i+1] == 'default':
                clas['rvals'] = get_default_rvals(clas0['dirname'])
            elif args[i+1] == 'all':
                clas['rvals'] = 'all'
            else:
                clas['rvals'] = read_cla_vals(args, i)
        elif arg == '--rvalscm':
            clas['rvals'] = read_cla_vals(args, i)/rsun
        elif arg == '--rrange':
            rbot, rtop, nrvals = read_cla_vals(args, i)
            nrvals = int(nrvals)
            clas['rvals'] = np.linspace(rtop, rbot, nrvals)

        # specify a desired time
        elif arg in ['--iter', '--prot', '--tdt', '--sec']:
            t_loc = float(args[i+1])
            di_trans = translate_times(t_loc, dirname, arg[2:])
            clas['val_iter'] = di_trans['val_iter']
        
        # desired quantity list (or group)
        elif arg == '--qvals': # able to specify either index or quantity name
            if not isall(args[i+1]): # if it's 'all', do nothing
                # qvals....make sure it's an integer array
                qvals = make_array(read_cla_vals(args, i))
                qvals = parse_quantities(qvals)[0]
                # leave qvals = 'all' alone; (calling script may want to 
                # use this for something specific
                # 04/13/22: at some point, find a cleaner way of doing this
            clas['qvals'] = qvals
        elif arg == '--latrange':
            latmin, latmax, nlatvals = read_cla_vals(args, i)
            nlatvals = int(nlatvals)
            clas['latvals'] = np.linspace(latmin, latmax, nlatvals)
        # everyone else is easy...
        elif '--' in arg:
            key = arg[2:]
            clas[key] = read_cla_vals(args, i)

    return dotdict(clas0), dotdict(clas)
