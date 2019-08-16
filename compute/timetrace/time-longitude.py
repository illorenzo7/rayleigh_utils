# Routine to trace Rayleigh G_Avgs data in time
# Created by: Loren Matilsky
# On: 08/15/2019
############################################################################
# This routine computes the trace in time/longitude of quantities in the 
# Shell_Slices data for a particular simulation. 
#
# By default, the 8 variables are computed at each time, in a latitude strip defined by
# (clat, dlat) = ([central latitude for average], [range of latitudes to average over])
# at the depths the shell slices were sampled at.
# The strip range can be changed using the options -clat and -dlat, e.g., 
# -clat 60 -dlat 30
#
# By default, the routine traces over all Shell_Slices, though user can specify an 
# alternate range, by, e.g.,
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
sys.path.append(os.environ['rasource'] + '/post_processing')
sys.path.append(os.environ['co'])
from rayleigh_diagnostics import AZ_Avgs, Shell_Slices
from common import get_file_lists, get_desired_range, strip_dirname,\
    get_widest_range_file, get_dict
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
newrange = False
for arg in args:
    if arg in ['-range', '-centerrange', '-leftrange', '-rightrange', '-n', '-f', '-all', '-iter']:
        newrange = True
if newrange:
    index_first, index_last = get_desired_range(int_file_list, args)

# Set other defaults
tag = ''
clat = 10.
dlat = 20. # by default average over the first 20 degrees of the Northern hemisphere
remove_diffrot = True
for i in range(nargs):
    arg = args[i]
    if arg == '-tag':
        tag = args[i+1] + '_'
    elif arg == '-clat':
        clat = float(args[i+1])
    elif arg == '-dlat':
        dlat = float(args[i+1])
    elif arg == '-keepdr':
        remove_diffrot = False
        
# Set the timetrace savename by the directory, what we are saving, 
# and first and last iteration files for the trace (and optional tag)
if clat > 0.:
    hemisphere = 'N'
else:
    hemisphere = 'S'
    
savename = dirname_stripped + ('_time-longitude_clat%s%02.0f_dlat%02.0f' %(hemisphere, np.abs(clat), dlat)) +\
        '_' + tag + file_list[index_first] + '_' + file_list[index_last] + '.pkl'
savefile = datadir + savename    
print('Your data will be saved in the file %s.' %savename)

# Read in first Shell_Slices file 
a0 = Shell_Slices(radatadir + file_list[index_first], '')

# Also read in time-latitude file (needed for differential rotation)
tl_file = get_widest_range_file(datadir, 'time-latitude')
print("Getting time-latitude trace data from %s..." %tl_file)
di_tl = get_dict(datadir + tl_file)
vals_tl = di_tl['vals']
qvals_tl = di_tl['qvals']
times_tl = di_tl['times']

# Read in grid info from time-latitude trace data
rr, ri, ro = di_tl['rr'], di_tl['ri'], di_tl['ro']
sint, cost, tt_lat = di_tl['sint'], di_tl['cost'], di_tl['tt_lat']
nt, nr = len(sint), len(rr)
nphi = 2*nt
lons = np.arange(0., 360., 360/nphi)

# Desired latitude range
ith1 = np.argmin(np.abs(tt_lat - (clat - dlat/2.)))
ith2 = np.argmin(np.abs(tt_lat - (clat + dlat/2.)))
lats_strip = tt_lat[ith1:ith2+1]


# Compute the timetrace of the D.R. at our particular latitude 
# and at the range of depths in for the shell slices
iq_vphi = np.argmin(np.abs(np.array(qvals_tl) - 3))
vphi_tl = vals_tl[:, :, :, iq_vphi]
rr_tl = rr[di_tl['rinds']]
nr_tl = len(rr_tl)
rr_tl_3d = rr_tl.reshape((1, 1, nr_tl))
sint_3d = sint.reshape((1, nt, 1))
Omega_tl = vphi_tl/(rr_tl_3d*sint_3d)
Omega_vs_time = np.mean(Omega_tl[:, ith1:ith2+1, :], axis=1)

print ('Considering Shell_Slices files %s through %s for the trace ...' %(file_list[index_first], file_list[index_last]))

iter1, iter2 = int_file_list[index_first], int_file_list[index_last]

vals = []
                     
times = []
iters = []

phi0 = np.zeros(nr_tl, dtype='float') # keep track of advection by the DR--remember 
                        # different advections for each depth
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
        # Please please PLEASE sample the same depths for both Shell_Slices and time-latitude traces!
        vals_strip = a.vals[:, ith1:ith2+1, :, :, j]
        vals_av = np.mean(vals_strip, axis=1)
        
        # Can update the times/iters lists now; [vals_av] may need to be "rolled"
        # (to deal with advection by the D.R.) before being added to [vals]
        times.append(a.time[j])
        iters.append(a.iters[j])
        
        # Closest time in times_tl to the current time
        it_tl = np.argmin(np.abs(times_tl - a.time[j]))
        
        # Loop over radius and subtract off the differential rotation by "rolling" the phi axis
        # Please please PLEASE sample the same depths for both Shell_Slices and time-latitude traces!
        for ir in range(nr_tl):                                        
            # Get the average rotation rate and figure out how far to advect phi0
            Omega_now = Omega_vs_time[it_tl, ir]
            if count > 0:
                phi0[ir] += Omega_now*(times[count] - times[count - 1])

                # roll the averaged strip along phi (backward by the deflection of phi0)
                nroll = int(phi0[ir]/(2*np.pi)*nphi)
                vals_av[:, ir, :] = np.roll(vals_av[:, ir, :], -nroll, axis=0)
        
        vals.append(vals_av.tolist())                      
        count += 1
        
# Convert lists into arrays
vals = np.array(vals)
times = np.array(times)
iters = np.array(iters)

# Miscellaneous metadata 
niter = len(iters)
nq = a.nq
depths = di_tl['depths']
ndepths = di_tl['ndepths']
rinds = di_tl['rinds']
rr_depth = di_tl['rr_depth']
rr_height = di_tl['rr_height']
d = ro - ri
tt = di_tl['tt']
rr_2d = di_tl['rr_2d']
tt_2d = di_tl['tt_2d']
sint_2d = di_tl['sint_2d']
cost_2d = di_tl['cost_2d']
xx = di_tl['xx']
zz = di_tl['zz']

# The computer congratulates itself on a job well done!
print ('Traced over %i Shell_Slices ...' %niter)

# Save the avarage
print ('Saving file at ' + savefile + ' ...')
f = open(savefile, 'wb')
pickle.dump({'vals': vals, 'times': times, 'iters': iters, 'qvals': a0.qv,\
    'rinds': rinds, 'lut': a0.lut, 'niter': niter, 'ndepths': ndepths, 'nq': nq,\
    'iter1': iter1, 'iter2': iter2, 'rr': rr, 'rr_depth': rr_depth,\
    'rr_height': rr_height, 'nr': nr, 'ri': ri, 'ro': ro, 'd': d, 'tt': tt,\
    'tt_lat': tt_lat, 'sint': sint, 'cost': cost,'nt': nt,\
    'rr_2d': rr_2d, 'tt_2d': tt_2d, 'sint_2d': sint_2d, 'cost_2d': cost_2d,\
    'xx': xx, 'zz': zz, 'lats_strip': lats_strip, 'clat': clat, 'dlat': dlat,\
    'lons':lons, 'nphi': nphi}, f, protocol=4)
f.close()