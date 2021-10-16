# Routine to trace Rayleigh Shell_Slices mag. field data in time
# Created by: Loren Matilsky
# On: 02/03/2020
############################################################################

# Import relevant modules
import numpy as np
import pickle
import sys, os
sys.path.append(os.environ['rapp'])
sys.path.append(os.environ['raco'])
from rayleigh_diagnostics import Shell_Slices, AZ_Avgs
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
for arg in args:
    if arg in ['-range', '-centerrange', '-leftrange', '-rightrange', '-n',\
            '-f', '-all', '-iter']:
        index_first, index_last = get_desired_range(int_file_list, args)

savename = dirname_stripped + '_partial_wreath_trace_' +\
        file_list[index_first] + '_' + file_list[index_last] + '.pkl'
savefile = datadir + savename    
print('Your data will be saved in the file %s' %savename)

# Read in first Shell_Slices/AZ_Avgs file 
a0 = Shell_Slices(radatadir + file_list[index_first], '')
az0 = AZ_Avgs(dirname + '/AZ_Avgs/' + file_list[index_first], '')

# Read in grid info from AZ_Avgs slice
rr = az0.radius
ri, ro = np.min(rr), np.max(rr)
sint, cost = az0.sintheta, az0.costheta
tt = np.arccos(cost)
tt_lat = (np.pi/2 - tt)*180./np.pi
nt, nr = len(sint), len(rr)
nphi = 2*nt
lons = np.arange(0., 360., 360/nphi)

# Desired latitude range
lat1 = -25.
lat2 = 25.
ith1 = np.argmin(np.abs(tt_lat - lat1))
ith2 = np.argmin(np.abs(tt_lat - lat2))
ith0 = np.argmin(np.abs(tt_lat))

# Start building the time-longitude traces
print ('Considering Shell_Slices files %s through %s for the trace ...'\
        %(file_list[index_first], file_list[index_last]))

iter1, iter2 = int_file_list[index_first], int_file_list[index_last]

# Split field by hemisphere (sticking to low-latitude), sign of Bphi field,
# and average of field/rms of field
Bphi_av_S_plus = []
Bphi_av_N_plus = []
Bphi_rms_S_plus = []
Bphi_rms_N_plus = []
Bphi_av_S_minus = []
Bphi_av_N_minus = []
Bphi_rms_S_minus = []
Bphi_rms_N_minus = []

Br_av_S_plus = []
Br_av_N_plus = []
Br_rms_S_plus = []
Br_rms_N_plus = []
Br_av_S_minus = []
Br_av_N_minus = []
Br_rms_S_minus = []
Br_rms_N_minus = []

Bt_av_S_plus = []
Bt_av_N_plus = []
Bt_rms_S_plus = []
Bt_rms_N_plus = []
Bt_av_S_minus = []
Bt_av_N_minus = []
Bt_rms_S_minus = []
Bt_rms_N_minus = []

# also filling factors
ff_S_plus = []
ff_N_plus = []
ff_S_minus = []
ff_N_minus = []

times = []
iters = []

count = 0 # don't know a priori how many times there will be to sample, 
            # so "count" as we go

def av(arr):
    # Take average of array along phi/theta dimensions, 
    # leaving 1d array with length nr
    return np.mean(np.mean(arr, axis=0), axis=0)

def rms(arr):
    # Take rms of array along phi/theta dimensions, 
    # leaving 1d array with length nr
    nphi_loc, ntheta_loc, nr_loc = np.shape(arr)
    square = np.sum(np.sum(arr**2, axis=0), axis=0)
    return np.sqrt(square/(nphi_loc*ntheta_loc))

for i in range(index_first, index_last + 1):
    print ('Adding Shell_Slices/%s to the trace ...' %file_list[i])
    if i == index_first:
        a = a0
    else:   
        a = Shell_Slices(radatadir + file_list[i], '')
                     
    ntimes_loc = a.niter # I think this is always 1...but may as well be
                        # general
    for j in range(ntimes_loc):
        Bphi_S = a.vals[:, ith1:ith0+1, :, a.lut[803], j]
        Bphi_N = a.vals[:, ith0:ith2+1, :, a.lut[803], j]
        where_plus_S = np.where(Bphi_S >= 0.)
        where_minus_S = np.where(Bphi_S < 0.)
        where_plus_N = np.where(Bphi_N >= 0.)
        where_minus_N = np.where(Bphi_N < 0.)

        Br_S = a.vals[:, ith1:ith0+1, :, a.lut[801], j]
        Br_N = a.vals[:, ith0:ith2+1, :, a.lut[801], j]
        Bt_S = a.vals[:, ith1:ith0+1, :, a.lut[802], j]
        Bt_N = a.vals[:, ith0:ith2+1, :, a.lut[802], j]

        Bphi_S_plus = np.zeros_like(Bphi_S)
        Bphi_S_minus = np.zeros_like(Bphi_S)
        Bphi_N_plus = np.zeros_like(Bphi_N)
        Bphi_N_minus = np.zeros_like(Bphi_N)

        Br_S_plus = np.zeros_like(Bphi_S)
        Br_S_minus = np.zeros_like(Bphi_S)
        Br_N_plus = np.zeros_like(Bphi_N)
        Br_N_minus = np.zeros_like(Bphi_N)

        Bt_S_plus = np.zeros_like(Bphi_S)
        Bt_S_minus = np.zeros_like(Bphi_S)
        Bt_N_plus = np.zeros_like(Bphi_N)
        Bt_N_minus = np.zeros_like(Bphi_N)

        ones_S_plus = np.zeros_like(Bphi_S)
        ones_S_minus = np.zeros_like(Bphi_S)
        ones_N_plus = np.zeros_like(Bphi_N)
        ones_N_minus = np.zeros_like(Bphi_N)

        Bphi_S_plus[where_plus_S] = Bphi_S[where_plus_S]
        Bphi_S_minus[where_minus_S] = Bphi_S[where_minus_S]
        Bphi_N_plus[where_plus_N] = Bphi_N[where_plus_N]
        Bphi_N_minus[where_minus_N] = Bphi_N[where_minus_N]

        Br_S_plus[where_plus_S] = Br_S[where_plus_S]
        Br_S_minus[where_minus_S] = Br_S[where_minus_S]
        Br_N_plus[where_plus_N] = Br_N[where_plus_N]
        Br_N_minus[where_minus_N] = Br_N[where_minus_N]

        Bt_S_plus[where_plus_S] = Bt_S[where_plus_S]
        Bt_S_minus[where_minus_S] = Bt_S[where_minus_S]
        Bt_N_plus[where_plus_N] = Bt_N[where_plus_N]
        Bt_N_minus[where_minus_N] = Bt_N[where_minus_N]

        ones_S_plus[where_plus_S] = 1.
        ones_S_minus[where_minus_S] = 1.
        ones_N_plus[where_plus_N] = 1.
        ones_N_minus[where_minus_N] = 1.

        Bphi_av_S_plus_loc = av(Bphi_S_plus)
        Bphi_av_N_plus_loc = av(Bphi_N_plus)
        Bphi_rms_S_plus_loc = rms(Bphi_S_plus)
        Bphi_rms_N_plus_loc = rms(Bphi_N_plus)
        Bphi_av_S_minus_loc = av(Bphi_S_minus)
        Bphi_av_N_minus_loc = av(Bphi_N_minus)
        Bphi_rms_S_minus_loc = rms(Bphi_S_minus)
        Bphi_rms_N_minus_loc = rms(Bphi_N_minus)

        Br_av_S_plus_loc = av(Br_S_plus)
        Br_av_N_plus_loc = av(Br_N_plus)
        Br_rms_S_plus_loc = rms(Br_S_plus)
        Br_rms_N_plus_loc = rms(Br_N_plus)
        Br_av_S_minus_loc = av(Br_S_minus)
        Br_av_N_minus_loc = av(Br_N_minus)
        Br_rms_S_minus_loc = rms(Br_S_minus)
        Br_rms_N_minus_loc = rms(Br_N_minus)

        Bt_av_S_plus_loc = av(Bt_S_plus)
        Bt_av_N_plus_loc = av(Bt_N_plus)
        Bt_rms_S_plus_loc = rms(Bt_S_plus)
        Bt_rms_N_plus_loc = rms(Bt_N_plus)
        Bt_av_S_minus_loc = av(Bt_S_minus)
        Bt_av_N_minus_loc = av(Bt_N_minus)
        Bt_rms_S_minus_loc = rms(Bt_S_minus)
        Bt_rms_N_minus_loc = rms(Bt_N_minus)
           
        times.append(a.time[j])
        iters.append(a.iters[j])

        Bphi_av_S_plus.append(Bphi_av_S_plus_loc.tolist())
        Bphi_av_N_plus.append(Bphi_av_N_plus_loc.tolist())
        Bphi_rms_S_plus.append(Bphi_rms_S_plus_loc.tolist())
        Bphi_rms_N_plus.append(Bphi_rms_N_plus_loc.tolist())
        Bphi_av_S_minus.append(Bphi_av_S_minus_loc.tolist())
        Bphi_av_N_minus.append(Bphi_av_N_minus_loc.tolist())
        Bphi_rms_S_minus.append(Bphi_rms_S_minus_loc.tolist())
        Bphi_rms_N_minus.append(Bphi_rms_N_minus_loc.tolist())

        Br_av_S_plus.append(Br_av_S_plus_loc.tolist())
        Br_av_N_plus.append(Br_av_N_plus_loc.tolist())
        Br_rms_S_plus.append(Br_rms_S_plus_loc.tolist())
        Br_rms_N_plus.append(Br_rms_N_plus_loc.tolist())
        Br_av_S_minus.append(Br_av_S_minus_loc.tolist())
        Br_av_N_minus.append(Br_av_N_minus_loc.tolist())
        Br_rms_S_minus.append(Br_rms_S_minus_loc.tolist())
        Br_rms_N_minus.append(Br_rms_N_minus_loc.tolist())

        Bt_av_S_plus.append(Bt_av_S_plus_loc.tolist())
        Bt_av_N_plus.append(Bt_av_N_plus_loc.tolist())
        Bt_rms_S_plus.append(Bt_rms_S_plus_loc.tolist())
        Bt_rms_N_plus.append(Bt_rms_N_plus_loc.tolist())
        Bt_av_S_minus.append(Bt_av_S_minus_loc.tolist())
        Bt_av_N_minus.append(Bt_av_N_minus_loc.tolist())
        Bt_rms_S_minus.append(Bt_rms_S_minus_loc.tolist())
        Bt_rms_N_minus.append(Bt_rms_N_minus_loc.tolist())

        ff_S_plus.append(av(ones_S_plus).tolist())
        ff_S_minus.append(av(ones_S_minus).tolist())
        ff_N_plus.append(av(ones_N_plus).tolist())
        ff_N_minus.append(av(ones_N_minus).tolist())

        count += 1        

# Convert lists into arrays
times = np.array(times)
iters = np.array(iters)

Bphi_av_S_plus = np.array(Bphi_av_S_plus)
Bphi_av_N_plus = np.array(Bphi_av_N_plus)
Bphi_rms_S_plus = np.array(Bphi_rms_S_plus)
Bphi_rms_N_plus = np.array(Bphi_rms_N_plus)
Bphi_av_S_minus = np.array(Bphi_av_S_minus)
Bphi_av_N_minus = np.array(Bphi_av_N_minus)
Bphi_rms_S_minus = np.array(Bphi_rms_S_minus)
Bphi_rms_N_minus = np.array(Bphi_rms_N_minus)

Br_av_S_plus = np.array(Br_av_S_plus)
Br_av_N_plus = np.array(Br_av_N_plus)
Br_rms_S_plus = np.array(Br_rms_S_plus)
Br_rms_N_plus = np.array(Br_rms_N_plus)
Br_av_S_minus = np.array(Br_av_S_minus)
Br_av_N_minus = np.array(Br_av_N_minus)
Br_rms_S_minus = np.array(Br_rms_S_minus)
Br_rms_N_minus = np.array(Br_rms_N_minus)

Bt_av_S_plus = np.array(Bt_av_S_plus)
Bt_av_N_plus = np.array(Bt_av_N_plus)
Bt_rms_S_plus = np.array(Bt_rms_S_plus)
Bt_rms_N_plus = np.array(Bt_rms_N_plus)
Bt_av_S_minus = np.array(Bt_av_S_minus)
Bt_av_N_minus = np.array(Bt_av_N_minus)
Bt_rms_S_minus = np.array(Bt_rms_S_minus)
Bt_rms_N_minus = np.array(Bt_rms_N_minus)

ff_S_plus = np.array(ff_S_plus)
ff_N_plus = np.array(ff_N_plus)
ff_S_minus = np.array(ff_S_minus)
ff_N_minus = np.array(ff_N_minus)

# Miscellaneous metadata 
niter = len(iters)
rinds = a0.inds
rvals = a0.radius
nrvals = a0.nr

# The computer congratulates itself on a job well done!
print ('Traced over %i Shell_Slices ...' %niter)

# Save the avarage
print ('Saving file at ' + savefile)
f = open(savefile, 'wb')
pickle.dump({'times': times, 'iters': iters,\
        'rinds': rinds, 'rvals': rvals, 'nrvals': nrvals,\
        'niter': niter,  'iter1': iter1, 'iter2': iter2, 'rr': rr,\
    'nr': nr, 'ri': ri, 'ro': ro, 'tt': tt, 'tt_lat': tt_lat,\
    'sint': sint, 'cost': cost,'nt': nt, 'lons':lons, 'nphi': nphi,\
    'lat1': tt_lat[ith1], 'lat2': tt_lat[ith2], 'lat0': tt_lat[ith0],\

    'Bphi_av_S_plus': Bphi_av_S_plus,\
    'Bphi_av_N_plus': Bphi_av_N_plus,\
    'Bphi_rms_S_plus': Bphi_rms_S_plus,\
    'Bphi_rms_N_plus': Bphi_rms_N_plus,\
    'Bphi_av_S_minus': Bphi_av_S_minus,\
    'Bphi_av_N_minus': Bphi_av_N_minus,\
    'Bphi_rms_S_minus': Bphi_rms_S_minus,\
    'Bphi_rms_N_minus': Bphi_rms_N_minus,\

    'Br_av_S_plus': Br_av_S_plus,\
    'Br_av_N_plus': Br_av_N_plus,\
    'Br_rms_S_plus': Br_rms_S_plus,\
    'Br_rms_N_plus': Br_rms_N_plus,\
    'Br_av_S_minus': Br_av_S_minus,\
    'Br_av_N_minus': Br_av_N_minus,\
    'Br_rms_S_minus': Br_rms_S_minus,\
    'Br_rms_N_minus': Br_rms_N_minus,\

    'Bt_av_S_plus': Bt_av_S_plus,\
    'Bt_av_N_plus': Bt_av_N_plus,\
    'Bt_rms_S_plus': Bt_rms_S_plus,\
    'Bt_rms_N_plus': Bt_rms_N_plus,\
    'Bt_av_S_minus': Bt_av_S_minus,\
    'Bt_av_N_minus': Bt_av_N_minus,\
    'Bt_rms_S_minus': Bt_rms_S_minus,\
    'Bt_rms_N_minus': Bt_rms_N_minus,\

    'ff_S_plus': ff_S_plus,\
    'ff_N_plus': ff_N_plus,\
    'ff_S_minus': ff_S_minus,\
    'ff_N_minus': ff_N_minus}, f, protocol=4)
f.close()
