# Date created: 02/04/2019
import numpy as np
import sys, os
from common import get_file_lists, get_widest_range_file, strip_dirname, get_dict
from get_parameter import get_parameter

def translate_times(time, dirname, translate_from='prot'):
    datadir = dirname + '/data/'

    trace_G_Avgs_file = get_widest_range_file(datadir, 'trace_G_Avgs')
    di = get_dict(datadir + trace_G_Avgs_file)

    times = di['times']
    iters = di['iters']

    # Get global rotation rate, if present
    angular_velocity = get_parameter(dirname, 'angular_velocity')
    Prot = 2*np.pi/angular_velocity

    if translate_from == 'prot':
        ind = np.argmin(np.abs(times/Prot - time))
    elif translate_from == 'iter':
        ind = np.argmin(np.abs(iters - time))
    elif translate_from == 'day':
        ind = np.argmin(np.abs(times/86400. - time))
    elif translate_from == 'sec':
        ind = np.argmin(np.abs(times - time))

    val_sec = times[ind]
    val_iter = iters[ind]
    val_day = times[ind]/86400.
    val_prot = times[ind]/Prot

    return dict({'val_sec': val_sec,'val_iter': val_iter,\
            'val_day': val_day, 'val_prot': val_prot})
