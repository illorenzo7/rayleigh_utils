# Routines to deal with command-line arguments (CLAs)
# It's a long one!
# Created: 04/17/2021

import numpy as np
#from common import get_parameter, get_domain_bounds, array_of_strings, is_an_int, rsun
from common import *
from varprops import get_quantity_group

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

def read_clas(args):
    # first get basic info
    clas0 = dict({})
    clas0['routinename'] = args[0].split('/')[-1][:-3]
    clas0['dirname'] = args[1]
    clas0['datadir'] = clas0['dirname'] + '/data/'
    clas0['plotdir'] = clas0['dirname'] + '/plots/'
    clas0['saveplot'] = True
    clas0['showplot'] = True
    clas0['tag'] = ''

    # see if magnetism is on
    magnetism = get_parameter(clas0['dirname'], 'magnetism')

    # get the other arguments
    # tag is default to make things easy
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
            clas['tag'] = '_' + args[i+1]
        elif arg == '--depths':
            clas['rvals'] = ro - read_cla_vals(args, i)*d
        elif arg == '--depthscz':
            dcz = ro - rm
            clas['rvals'] = ro - read_cla_vals(args, i)*dcz
        elif arg == '--depthirz':
            drz = rm - ri
            clas['rvals'] = rm - read_cla_vals(args, i)*drz
        elif arg == '--rvals':
            if args[i+1] == 'default':
                clas['rvals'] = get_default_rvals(clas0['dirname'])*rsun
            else:
                clas['rvals'] = read_cla_vals(args, i)*rsun
        elif arg == '--rvalscm':
            clas['rvals'] = read_cla_vals(args, i)
        elif arg == '--rrange':
            rbot, rtop, nrvals = read_cla_vals(args, i)
            rbot *= rsun
            rtop *= rsun
            nrvals = int(nrvals)
            clas['rvals'] = np.linspace(rtop, rbot, nrvals)
        elif arg == '--nrperzone':
            rvals = np.array([], dtype='float')
            nrperzone = int(args[i+1])
            basedepths = np.arange(1.0/(nperzone+1.0), 1.0,\
                    1.0/(nperzone+1.0))
            for idomain in range(ndomains):
                rbot = domain_bounds[ndomains - idomain - 1]
                rtop = domain_bounds[ndomains - idomain]
                if idomain == ndomains - 1:
                    rvals_to_add = rtop - (rtop - rbot)*basedepths[:-1]
                else:
                    rvals_to_add = rtop - (rtop - rbot)*basedepths
                rvals = np.hstack((rvals, rvals_to_add))
            clas['rvals'] = rvals
        elif arg == '--qvals':
            argvals = read_cla_vals(args, i)
            if np.isscalar(argvals): # either group or one int
                if is_an_int(argvals):
                    clas['qvals'] = np.array([int(argvals)])
                    clas['titles'] = array_of_strings(clas['qvals'])
                    clas['units'] = 'cgs'
                else:
                    the_qgroup = get_quantity_group(argvals, magnetism)
                    clas['qvals'] = the_qgroup['qvals']
                    clas['titles'] = the_qgroup['titles']
                    clas['units'] = the_qgroup['units']
                    clas['ncol'] = the_qgroup['ncol']
                    clas['tag'] = '_' + argvals
                    clas['totsig'] = the_qgroup['totsig']
            else:
                # this was a list of integers
                clas['qvals'] = read_cla_vals(args, i)
                clas['titles'] = array_of_strings(clas['qvals'])
                clas['units'] = 'cgs'
        elif arg == '--latrange':
            latmin, latmax, nlatvals = read_cla_vals(args, i)
            nlatvals = int(nlatvals)
            clas['latvals'] = np.linspace(latmin, latmax, nlatvals)
        # everyone else is easy...
        elif '--' in arg:
            key = arg[2:]
            clas[key] = read_cla_vals(args, i)

    if 'rvals' in clas:
        clas['rvals'] /= rsun # always keep rvals in solar radius units
        if np.isscalar(clas['rvals']):
            clas['rvals'] = np.array([clas['rvals']]) 
            # rvals should always be an array
    return clas0, clas
