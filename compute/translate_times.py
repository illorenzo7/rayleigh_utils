# Date created: 02/04/2019
import numpy as np
import sys, os
from common import get_file_lists, get_widest_range_file, strip_dirname, get_dict
from get_parameter import get_parameter
from rayleigh_diagnostics import G_Avgs, Shell_Slices
from time_scales import compute_Prot, compute_tdt

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
        if nfiles == 0: # probably this is a movie, for which there are no G_Avgs, but Shell_Slices
            file_list, int_file_list, nfiles = get_file_lists(dirname + '/Shell_Slices')
            print ("translate_times(): translating using Shell_Slices data")
            funct = Shell_Slices
            radatadir = dirname + '/Shell_Slices'
        else:
            print ("translate_times(): translating using G_Avgs data")
            funct = G_Avgs
            radatadir = dirname + '/G_Avgs'
        iiter = np.argmin(np.abs(int_file_list - time))
        a = funct(radatadir + '/' + file_list[iiter], '')
        jiter = np.argmin(np.abs(a.iters - time))
        val_sec = a.time[jiter]
        val_iter = a.iters[jiter]
        val_day = val_sec/86400.
        val_unit = val_sec/time_unit

    else: # otherwise, hopefully you computed some time traces beforehand!
        # Get the data directory
        datadir = dirname + '/data/'
     
        try:        
            the_file = get_widest_range_file(datadir, 'trace_G_Avgs')
            di = get_dict(datadir + the_file)
            print ("translate_times(): translating using trace_G_Avgs file")
        except:
            the_file = get_widest_range_file(datadir, 'time-latitude')
            di = get_dict(datadir + the_file)
            print ("translate_times(): translating using time-latitude file")

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
