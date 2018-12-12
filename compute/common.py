# Author: Loren Matilsky
# Last modified: 11/10/2018
# Common routines for routines for post-processing Rayleigh data 

import numpy as np
import os

# Read in all files from the Rayleigh data directory and sort them by name (number)
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

def get_desired_range(int_file_list, args):
    nargs = len(args)
    nfiles = len(int_file_list)
    for i in range(nargs):
        arg = args[i]
        if (arg == '-range'):    # Give user option to specify range of iterations to
                                 # average over. User can specify 'first' after '-range'
                                 # to begin averaging at the first available file;
                                 # User can specify 'last' (2 arguments after '-range')
                                 # to average until the last available file.
            desired_first_iter = int(args[i+1])
            if (desired_first_iter == 'first'):
                desired_first_iter = int_file_list[0]
            else:
                desired_first_iter = int(desired_first_iter) 
                
            desired_last_iter = args[i+2]
            if (desired_last_iter == 'last'):
                desired_last_iter = int_file_list[-1]
            else:
                desired_last_iter = int(desired_last_iter)
            # Find the index in the int_file_list closest to the user's preferences
            index_first = np.argmin(np.abs(desired_first_iter - int_file_list))
            index_last = np.argmin(np.abs(desired_last_iter - int_file_list))
        elif (arg == '-n'): # allow the user to specify a desired number of iterations
                            # to average over, ending with the last data file
            number_to_average = int(args[i+1])
            index_first = nfiles - number_to_average
            if index_first < 0:
                index_first = 0
            index_last = nfiles - 1
        elif (arg == '-centerrange'):
            desired_central_iter = int(args[i + 1])
            ndatafiles = int(args[i + 2])
            central_index = np.argmin(np.abs(desired_central_iter - int_file_list))
            if ndatafiles % 2 == 0: #ndatafiles is even
                index_first = central_index - ndatafiles//2 + 1
                index_last = central_index + ndatafiles//2
            else:  #ndatafiles is odd
                index_first = central_index - ndatafiles//2
                index_last = central_index + ndatafiles//2
        elif (arg == '-all'):
            index_first = 0
            index_last = nfiles - 1

    if (index_first < 0): index_first = 0
    if (index_last > nfiles - 1): index_last = nfiles - 1

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
    # Find the desired file(s) in the data directory. If there are multiple, by
    # default choose the one with widest range in the trace/average/distribution
    datafiles = os.listdir(datadir)
    len_name = len(data_name)
    specific_files = []
    for i in range(len(datafiles)):
        datafile = datafiles[i]
        if data_name in datafile:
            specific_files.append(datafile)
#            istart = datafile.find(data_name)
#            possible_iter = datafile[istart + len_name + 1:istart + len_name + 9]
#            if is_an_int(possible_iter):
#                specific_files.append(datafile)

    ranges = []
    iters1 = []
    iters2 = []
    for specific_file in specific_files:
        specific_file_stripped = specific_file[:-4] # get rid of '.npy'...
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
