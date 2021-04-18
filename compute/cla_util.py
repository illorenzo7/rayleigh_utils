# Routines to deal with command-line arguments (CLAs)
# It's a long one!
# Created: 04/17/2021

import numpy as np
from common import get_parameter, get_domain_bounds, array_of_strings
from quantity_groups import get_quantity_group

# set default CLAs

# data averaging stuff
clas_default = dict({})
clas_default['datadir'] = None
clas_default['radtype'] = 'azav'
clas_default['tag'] = ''
clas_default['qvals'] = None

# plotting stuff
clas_default['plotdir'] = None
clas_default['nlevs'] = 20
clas_default['rbcz'] = None
clas_default['rvals'] = None

clas_default['minmax'] = None
clas_default['minmaxrz'] = None
clas_default['ymin'] = None
clas_default['ymax'] = None
clas_default['xminmax'] = None
clas_default['xmin'] = None
clas_default['xmax'] = None

clas_default['the_file'] = None
clas_default['linthresh'] = None
clas_default['linscale'] = None
clas_default['linthreshrz'] = None
clas_default['linscalerz'] = None
clas_default['symlog'] = False
clas_default['plotcontours'] = True
clas_default['plotlatlines'] = True
clas_default['plotboundary'] = True
clas_default['saveplot'] = True
clas_default['showplot'] = True
clas_default['latvals'] = None

default_latvals = np.array([-85., -75., -60., -45., -30., -15., 0., 15., 30., 45., 60., 75., 85.])

def get_default_rvals(dirname):
    ncheby, domain_bounds = get_domain_bounds(dirname)
    ndomains = len(ncheby)
    ri, rm, ro = domain_bounds
    basedepths = np.array([0.05, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 0.95, 1.0])
    rvals = np.array([], dtype='float')
    for idomain in range(ndomains):
        rbot = domain_bounds[ndomains - idomain - 1]
        rtop = domain_bounds[ndomains - idomain]
        if idomain == ndomains - 1:
            rvals_to_add = rtop - (rtop - rbot)*basedepths[:-1]
        else:
            rvals_to_add = rtop - (rtop - rbot)*basedepths
        rvals = np.hstack((rvals, rvals_to_add))
    return rvals/rsun

def read_cla_vals(args, i, dtype='float'):
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
    val_args = args_after[:iend]
    if dtype == 'float':
        convert = float
    if dtype == 'int':
        convert = int
    if dtype == 'str':
        convert = str

    # store the cla value in the "vals" list
    vals = []
    for val_arg in val_args:
        for st in val_arg.split():
            vals.append(convert(st))

    # if the list has only one value, make it not a list
    if len(vals) == 1:
        vals = vals[0]
    return np.array(vals)

def read_cla_arbitrary(args, key, default=None, dtype='float'):
    nargs = len(args)
    out = default
    for i in range(nargs):
        arg = args[i]
        if arg == '--' + key:
            if i == nargs - 1: # must be boolean if appearing at end
                out = True
            else:
                if '--' in args[i+1]: # must be boolean if no following val
                    out = True
                else:
                    out = read_cla_vals(args, i, dtype=dtype)
    return out

def read_clas(dirname, args):
    # start with default CLAs, then change them
    clas = clas_default.copy()

    # see if magnetism is on
    magnetism = get_parameter(dirname, 'magnetism')

    nargs = len(args)
    for i in range(nargs):
        arg = args[i]
        if arg == '--datadir':
            clas['datadir'] = args[i+1]
        if arg == '--radtype':
            clas['radtype'] = args[i+1]
        if arg == '--plotdir':
            clas['plotdir'] = args[i+1]
        if arg == '--minmax':
            clas['minmax'] = read_cla_vals(args, i)
        if arg == '--ymin':
            clas['ymin'] = read_cla_vals(args, i)
        if arg == '--ymax':
            clas['ymax'] = read_cla_vals(args, i)
        if arg == '--xminmax':
            clas['xminmax'] = read_cla_vals(args, i)
        if arg == '--xmin':
            clas['xmin'] = read_cla_vals(args, i)
        if arg == '--xmax':
            clas['xmax'] = read_cla_vals(args, i)
        if arg == '-minmaxrz':
            minmaxrz = read_cla_vals(args, i)
        if arg == '--rbcz':
            clas['rbcz'] = float(args[i+1])
        if arg == '--nlevs':
            clas['nlevs'] = int(args[i+1])
        if arg == '--usefile':
            the_file = args[i+1]
            clas['the_file'] = the_file.split('/')[-1]
        if arg == '--nocontour':
            clas['plotcontours'] = False
        if arg == '--nobound':
            clas['plotboundary'] = False
        if arg == '--nolat':
            clas['plotlatlines'] = False
        if arg == '--symlog':
            clas['symlog'] = True
        if arg == '--linthresh':
            clas['linthresh'] = float(args[i+1])
        if arg == '--linscale':
            clas['linscale'] = float(args[i+1])
        if arg == '--linthreshrz':
            clas['linthreshrz'] = float(args[i+1])
        if arg == '--linscalerz':
            clas['linscalerz'] = float(args[i+1])
        if arg == '--tag':
            clas['tag'] = '_' + args[i+1]
        if arg == '--depths':
            clas['rvals'] = ro - read_cla_vals(args, i)*d
        if arg == '--depthscz':
            dcz = ro - rm
            clas['rvals'] = ro - read_cla_vals(args, i)*dcz
        if arg == '--depthirz':
            drz = rm - ri
            clas['rvals'] = rm - read_cla_vals(args, i)*drz
        if arg == '--rvals':
            if args[i+1] == 'default':
                clas['rvals'] = get_default_rvals(dirname)*rsun
            else:
                clas['rvals'] = read_cla_vals(args, i)*rsun
        if arg == '--rvalscm':
            clas['rvals'] = read_cla_vals(args, i)
        if arg == '--rrange':
            rbot, rtop, nrvals = read_cla_vals(args, i)
            nrvals = int(nrvals)
            clas['rvals'] = np.linspace(rtop, rbot, nrvals)
        if arg == '--nrperzone':
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
        if arg == '--qvals':
            if isinstance(args[i+1], str):
                the_qgroup = get_quantity_group(args[i+1], magnetism)
                clas['qvals'] = the_qgroup['qvals']
                clas['titles'] = the_qgroup['titles']
                clas['units'] = the_qgroup['units']
            else:
                clas['qvals'] = read_cla_vals(args, i, dtype='int')
                clas['titles'] = array_of_strings(qvals)
                clas['units'] = 'cgs'
        if arg == '--latvals':
            clas['latvals'] = read_cla_vals(args, i)
        if arg == '--latrange':
            latmin, latmax, nlatvals = read_cla_vals(args, i)
            nlatvals = int(nlatvals)
            clas['latvals'] = np.linspace(latmin, latmax, nlatvals)

    # deal with qvals if it is still None
    if clas['qvals'] is None:
        the_qgroup = get_quantity_group('default', magnetism)
        clas['qvals'] = the_qgroup['qvals']
        clas['titles'] = the_qgroup['titles']
        clas['units'] = the_qgroup['units']
    if not clas['rvals'] is None:
        clas['rvals'] /= rsun # always keep rvals in solar radius units
        if np.isscalar(clas['rvals']):
            clas['rvals'] = np.array([clas['rvals']]) 
            # rvals should always be an array
    return clas


