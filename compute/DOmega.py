# Author: Loren Matilsky
# Created: 05/18/2019
# This script computes the differential rotation contrast at the surface
# of the simulation from the equator to ~60 degrees, expressed as a 
# fraction of the frame rotation rate.
# Computes this fraction differential rotation for the Rayleigh run 
# directory indicated by [dirname]. To use an AZ_Avgs file
# different than the one associated with the longest averaging range, use
# -usefile [complete name of desired AZ_Avgs file]
# Saves plot in
# [dirname]_diffrot_[first iter]_[last iter].npy

import numpy as np
import pickle
import sys, os
sys.path.append(os.environ['rapp'])
sys.path.append(os.environ['raco'])
from common import get_widest_range_file, strip_dirname, get_dict

# Get directory name
dirname = sys.argv[1]

# Directory with data and plots, make the plotting directory if it doesn't
# already exist    
datadir = dirname + '/data/'
plotdir = dirname + '/plots/'
if (not os.path.isdir(plotdir)):
    os.makedirs(plotdir)

# Set defaults
my_boundstype = 'manual'
user_specified_minmax = False 
the_file = get_widest_range_file(datadir, 'AZ_Avgs')

# Read in CLAs (if any) to change default variable ranges and other options
args = sys.argv[2:]
nargs = len(args)
for i in range(nargs):
    arg = args[i]
    if arg == '-minmax':
        my_min, my_max = float(args[i+1]), float(args[i+2])
        user_specified_minmax = True
    elif arg == '-usefile':
        the_file = args[i+1]
        the_file = the_file.split('/')[-1]
        
# Read in AZ_Avgs data
print ('Getting data from ' + datadir + the_file)
di = get_dict(datadir + the_file)

iter1, iter2 = di['iter1'], di['iter2']
vals = di['vals']
lut = di['lut']

vr_av, vt_av, vp_av = vals[:, :, lut[1]], vals[:, :, lut[2]],\
        vals[:, :, lut[3]]

# Get necessary grid info
rr = di['rr']
cost = di['cost']
sint = di['sint']
tt_lat = di['tt_lat']
xx = di['xx']

# Get differential rotation in the rotating frame. 
Om = vp_av/xx
diffrot = Om*1.0e9/2/np.pi # rad/s --> nHz

# Get the frame rotation rate
Om0 = get_parameter(dirname, 'angular_velocity')
Om0 = Om0/(2*np.pi)*1e9 # rad/s --> nHz

# Compute the rotation rate difference between equator and 60 degrees
# at the outer surface
it0, it60 = np.argmin(np.abs(tt_lat - 0)), np.argmin(np.abs(tt_lat - 60))
minrot, maxrot = diffrot[it60, 0], diffrot[it0, 0]
DOm = maxrot - minrot

# Compute the differential rotation fraction!
fraction = DOm/Om0

# And print it
print("The differential rotation is %.1f nHz" %DOm)
print("The frame rate is %.1f nHz" %Om0)
print("The differential rotation fraction is %.4f" %fraction)

# Also write it as an empty file in [dirname]
# Remove all other "DOmega_is_" that might be present
# (in case this definition of DOmega is replaced and/or another
# time interval is used for the average
names = os.listdir(dirname)
for name in names:
    if "DOmega_is_" in name:
        os.remove(dirname + '/' + name)
fname = dirname + ("/00_DOmega_is_%.4f" %fraction)
f = open(fname, "w")
f.close()
