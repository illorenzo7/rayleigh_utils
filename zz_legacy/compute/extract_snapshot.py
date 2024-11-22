import sys, os, shutil
import numpy as np

def is_an_int(string):
    # obviously, first check if it's actually an int
    if isinstance(string, int) or isinstance(string, np.int64):
        return True
    # otherwise, loop through each character and see if it's an int
    len_str = len(string)
    bool_val = True # start assuming it's an int, then see if we're wrong
    for i in range(len_str):
        char = string[i]
        bool_val *= (char >= '0' and char <= '9')
    return(bool(bool_val))

def get_dataname_from_file(filename):
    # route does what it's named to do...gets the dataname associated with post-processed 
    # data product name
    just_file = filename.split('/')[-1] #strip the possible full path info
    return just_file.split('-')[0]

def get_iters_from_file(filename):
    # route does what it's named to do...
    filename_end = filename.split('-')[-1][:-4] 
    # (gets the [iter1]_[iter2].pkl and removes the trailing .ext)
    iters_st = filename_end.split('_')
    iter1, iter2 = int(iters_st[0]), int(iters_st[1])
    return iter1, iter2

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

def my_copy(fname1, fname2, direct=False):
    if os.path.exists(fname1) and not os.path.exists(fname2):
        if direct:
            shutil.copytree(fname1, fname2)
        else:
            shutil.copy(fname1, fname2)

def my_mkdir(dirname):
    if not os.path.isdir(dirname):
        os.makedirs(dirname)
    return dirname

dir1 = sys.argv[1]
dir2 = sys.argv[2]

# copy over basic files
for fname in ['main_input', 'equation_coefficients', 'grid_info', 'jobinfo.txt', 'custom_reference_binary']:
    my_copy(dir1 + '/' + fname, dir2 + '/' + fname)

# copy over last checkpoint
checkpoint_strings = os.listdir(dir1 + '/Checkpoints')
checkpoint_numbers = np.zeros_like(checkpoint_strings, dtype='int')
for i in range(len(checkpoint_strings)):
    st = checkpoint_strings[i]
    if is_an_int(st):
        checkpoint_numbers[i] = int(st)
i_last_numbered_checkpoint = np.argmax(checkpoint_numbers)
last_numbered_checkpoint = checkpoint_strings[i_last_numbered_checkpoint]
# hopefully this will be the same as "last_checkpoint" but may not be

my_copy(dir1 + '/Checkpoints/' + last_numbered_checkpoint,\
            dir2 + '/Checkpoints/' + last_numbered_checkpoint, direct=True)

# copy over checkpoint log info
for fname in ['checkpoint_log', 'last_checkpoint']:
    my_copy(dir1 + '/Checkpoints/' + fname, dir2 + '/Checkpoints/' + fname)

# copy over last output files
for datadir in ['G_Avgs', 'Shell_Avgs', 'AZ_Avgs', 'Shell_Spectra',\
                'Shell_Slices', 'Equatorial_Slices', 'Meridional_Slices']:
    datafile_strings = os.listdir(dir1 + '/' + datadir)
    datafile_numbers = np.zeros_like(datafile_strings, dtype='int')
    for i in range(len(datafile_strings)):
        st = datafile_strings[i]
        if is_an_int(st):
            datafile_numbers[i] = int(st)
    
    if len(datafile_numbers) > 0:
        i_last = np.argmax(datafile_numbers)
        last_datafile = datafile_strings[i_last]
        my_mkdir(dir2 + '/' + datadir)
        my_copy(dir1 + '/' + datadir + '/' + last_datafile,\
                dir2 + '/' + datadir + '/' + last_datafile)

# copy over data files
datadir1 = dir1 + '/data/'
datadir2 = dir2 + '/data/'
my_mkdir(datadir2)
for dataname in ['G_Avgs', 'G_Avgs_trace', 'Shell_Avgs', 'AZ_Avgs', 'Shell_Spectra']:
    the_file1 = get_widest_range_file(datadir1, dataname)
    if not the_file1 is None:
        the_file2 = the_file1.replace(dir1, dir2)
        my_copy(the_file1, the_file2)
