# Author: Loren Matilsky
# Created: 03/23/2020
# This script takes the time/longitudinally averaged <v_phi> and outputs
# it in the file "eq_vp" for Rayleigh later to read
import numpy as np
import sys, os
sys.path.append(os.environ['rapp'])
sys.path.append(os.environ['raco'])
from compute_grid_info import compute_grid_info
from common import *
from scipy.interpolate import interp1d
dirname = sys.argv[1]

fname = 'eq_vp'
azav_file = None
args = sys.argv[2:]
nargs = len(args)
outdir = dirname
ncheby_new = None
nt_new = None
for i in range(nargs):
    arg = args[i]
    if arg == '-usefile':
        azav_file = args[i+1]
    elif arg == '-fname':
        fname = args[i+1]
    elif arg == '-outdir':
        outdir = args[i+1]
    elif arg == '-ncheby':
        ncheby_li = args[i+1].split()
        for i in range(len(ncheby_li)):
            ncheby_li[i] = int(ncheby_li[i])
        ncheby_new = tuple(ncheby_li)
    elif arg == '-nt':
        nt_new = int(args[i+1])

# Read in the AZ_Avgs data
datadir = dirname + '/data/'
if azav_file is None:
    azav_file = get_widest_range_file(datadir, 'AZ_Avgs')
print("getting data from data/%s" %azav_file)
di = get_dict(datadir + azav_file)
vals = di['vals']
# Get <v_phi> and the associated grid
lut = di['lut']
mean_vp = vals[:, :, lut[3]]
nr = di['nr']
ncheby, domain_bounds = get_domain_bounds(dirname)
rr = di['rr']
nt = di['nt']
tt = di['tt']

if ncheby_new is None:
    # No new r - grid; the "new" grid equals the old grid
    try:
        rmin, rmax = get_parameter(dirname, 'rmin'),\
                get_parameter(dirname, 'rmax')
        nr_new = nr
        ncheby_new = (nr_new,)
    except:
        nr_new = nr
        ncheby_new = tuple(get_parameter(dirname, 'ncheby'))
    interpr = False
else:
    interpr = True

if nt_new is None:
    nt_new = nt
    interpt = False
else:
    interpt = True

# Compute the new grid if necessary and interpolate
#    nr, nt, nphi, r, rw, tt, cost, sint, tw, phi, dphi =\
if interpr or interpt:
    nr_new, nt_new, dummy, rr_new, dummy, tt_new, dummy, dummy, dummy,\
            dummy, dummy =\
        compute_grid_info(domain_bounds, ncheby_new, nt_new)

    # Always interpolate in both directions if either direction requested
    # if new grid = old grid in a particular direction, interpolation
    # in that direction will leave the data unchanged
    print ("Interpolation requested")
    print ("Original ncheby = ", ncheby)
    print ("Requested ncheby = ", ncheby_new)
    print ("Original ntheta = %i" %nt)
    print ("Requested ntheta = %i" %nt_new)

    # First interpolate in radius
    mean_vp_interp = np.zeros((nt, nr_new))
    for it in range(nt):
        func = interp1d(rr, mean_vp[it, :], kind='linear')
        mean_vp_interp[it, :] = func(rr_new)

    # Then interpolate in theta
    mean_vp_interp2 = np.zeros((nt_new, nr_new))
    for ir in range(nr_new):
        func = interp1d(tt, mean_vp_interp[:, ir], kind='linear',\
                fill_value = 'extrapolate')
        mean_vp_interp2[:, ir] = func(tt_new)

    mean_vp = mean_vp_interp2

# Write the data
print ("Writing <v_phi> to", outdir + '/' + fname)

f = open(outdir + '/' + fname, 'wb')
sigpi = np.array(314, dtype=np.int32)

f.write(sigpi.tobytes())
for ir in range(nr_new):
    for it in range(nt_new):
        f.write(mean_vp[it,ir].tobytes())
f.close()
