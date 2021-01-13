# Routine to trace Rayleigh G_Avgs data in time
# Created by: Loren Matilsky
# On: 08/15/2019
############################################################################
# This routine computes the trace in time/longitude of quantities in the 
# Shell_Slices data for a particular simulation. 
#
# By default, the 8 variables are computed at each time, in a latitude strip
# defined by (clat, dlat) = ([central latitude for average],
# [range of latitudes to average over]) at the depths the shell slices were 
# sampled. 
#
# The strip range can be changed using the options -clat and -dlat, e.g., 
# -clat 60 -dlat 30, for a strip averaged between 45 and 75 degrees (North)
#
# By default, the routine traces over all Shell_Slices in the directory,
# though user can specify an alternate range, by, e.g.,
# -n 10 (last 10 files)
# -range iter1 iter2 (no.s for start/stop data files; iter2 can be "last")
# -centerrange iter0 nfiles (trace about central file iter0 over nfiles)
#
# The final datacube output ('vals') will have shape
# (nphi, niter, nr, nq), where nr and nq are the attributes of the shell slices

# Import relevant modules
import numpy as np
import pickle
import sys, os
sys.path.append(os.environ['rapp'])
sys.path.append(os.environ['raco'])
from rayleigh_diagnostics import Shell_Slices, AZ_Avgs
from common import *
from get_parameter import get_parameter

# Get the name of the run directory
dirname = sys.argv[1]
# Get the stripped name to use in file naming
dirname_stripped = strip_dirname(dirname)

# Find the relevant place to store the data, and create the directory if it
# doesn't already exist
datadir = dirname + '/data/'
if (not os.path.isdir(datadir)):
    os.makedirs(datadir)

# We are interested in latitude strips from the shell slices
radatadir = dirname + '/Shell_Slices/'

# Get all the file names in datadir and their integer counterparts
file_list, int_file_list, nfiles = get_file_lists(radatadir)

# Read in CLAs
args = sys.argv[2:]
nargs = len(args)

# By default, have the iter indices range over full file range
index_first, index_last = 0, nfiles - 1
for arg in args:
    if arg in ['-range', '-centerrange', '-leftrange', '-rightrange', '-n',\
            '-f', '-all', '-iter']:
        index_first, index_last = get_desired_range(int_file_list, args)

# Set other defaults
tag = ''
clat = 10.
dlat = 0. # by default do NOT average over latitude
remove_diffrot = True
for i in range(nargs):
    arg = args[i]
    if arg == '-tag':
        tag = args[i+1] + '_'
    elif arg == '-clat':
        clat = float(args[i+1])
    elif arg == '-dlat':
        dlat = float(args[i+1])
        
# Set the timetrace savename by the directory, what we are saving, 
# and first and last iteration files for the trace (and optional tag)
if clat >= 0.:
    hemisphere = 'N'
else:
    hemisphere = 'S'
    
savename = dirname_stripped + '_time-longitude_' + tag +\
        ('clat%s%02.0f_dlat%03.0f' %(hemisphere, np.abs(clat), dlat)) +\
        '_' + file_list[index_first] + '_' + file_list[index_last] + '.pkl'
savefile = datadir + savename    
print('Your data will be saved in the file %s.' %savename)

# Read in first Shell_Slices/AZ_Avgs file 
a0 = Shell_Slices(radatadir + file_list[index_first], '')
az0 = AZ_Avgs(dirname + '/AZ_Avgs/' + file_list[index_first], '')

# Read in grid info from AZ_Avgs slice
rr =az0.radius
ri, ro = np.min(rr), np.max(rr)
sint, cost = az0.sintheta, az0.costheta
tt = np.arccos(cost)
tt_lat = (np.pi/2 - tt)*180./np.pi
nt, nr = len(sint), len(rr)
nphi = 2*nt
lons = np.arange(0., 360., 360/nphi)

# Desired latitude range
ith1 = np.argmin(np.abs(tt_lat - (clat - dlat/2.)))
ith2 = np.argmin(np.abs(tt_lat - (clat + dlat/2.)))
lats_strip = tt_lat[ith1:ith2+1]

# Start building the time-longitude traces
print ('Considering Shell_Slices files %s through %s for the trace ...'\
        %(file_list[index_first], file_list[index_last]))

iter1, iter2 = int_file_list[index_first], int_file_list[index_last]
                      
vals = []
times = []
iters = []

count = 0 # don't know a priori how many times there will be to sample, 
            # so "count" as we go
    
for i in range(index_first, index_last + 1):
    print ('Adding Shell_Slices/%s to the trace ...' %file_list[i])
    if i == index_first:
        a = a0
    else:   
        a = Shell_Slices(radatadir + file_list[i], '')
                     
    ntimes_loc = a.niter
    for j in range(ntimes_loc):
        vals_strip = a.vals[:, ith1:ith2+1, :, :, j]
        vals_av = np.mean(vals_strip, axis=1)
        
        times.append(a.time[j])
        iters.append(a.iters[j])
        vals.append(vals_av.tolist())                      
        count += 1        

# Convert lists into arrays
vals = np.array(vals)
times = np.array(times)
iters = np.array(iters)

# Miscellaneous metadata 
niter = len(iters)
nq = a.nq
rinds = a0.inds
rvals = a0.radius
nrvals = a0.nr

# The computer congratulates itself on a job well done!
print ('Traced over %i Shell_Slices ...' %niter)

# Save the avarage
print ('Saving file at ' + savefile + ' ...')
f = open(savefile, 'wb')
pickle.dump({'vals': vals, 'times': times, 'iters': iters, 'qvals': a0.qv,\
        'rinds': rinds, 'rvals': rvals, 'nrvals': nrvals, 'lut': a0.lut,\
        'niter': niter, 'nq': nq, 'iter1': iter1, 'iter2': iter2, 'rr': rr,\
    'nr': nr, 'ri': ri, 'ro': ro, 'tt': tt, 'tt_lat': tt_lat, 'sint': sint,\
    'cost': cost,'nt': nt, 'lats_strip': lats_strip, 'clat': clat,\
    'dlat': dlat, 'lons':lons, 'nphi': nphi}, f, protocol=4)
f.close()
