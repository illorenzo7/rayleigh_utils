# Author: Loren Matilsky
# Created: 02/25/2018
# Computes the thermal boundary layer thickness as defined in 
# Featherstone & Hindman 2016, eq. (27). Essentially we compute an "equivalent
# width" by 

import numpy as np
import sys
import os

# Get run directory
dirname = sys.argv[1]
radatadir = dirname + '/Shell_Avgs/'

datadir = dirname + '/data/'
if (not os.path.isdir(datadir)):
    os.makedirs(datadir)

entropy_dr = np.load(datadir + 's_dr_spherical_mean.npy')
min_entropy_dr = np.min(entropy_dr)

threshold_fraction = 0.06 # Some small number--we expect by the bottom of
        # the TBL, the entropy gradient, should be "nearly flat"
threshold_gradient = threshold_fraction*min_entropy_dr

# Get the grid structure
rr,tt,cost,sint,rr_depth,ri,ro,d = np.load(datadir + 'grid_info.npy')

# compute the two cells with gradients closest to the threshold value
diffs = np.abs(entropy_dr - threshold_gradient)
ir1 = np.argmin(diffs)

diffs_prime = np.copy(diffs)
diffs_prime[ir1] = np.inf
ir2 = np.argmin(diffs_prime)

li = [ir1, ir2]
li.sort()

r1 = rr[ir1]; r2 = rr[ir2]
entropy_dr1 = entropy_dr[ir1]; entropy_dr2 = entropy_dr[ir2]

# Get slope
m = (entropy_dr1 - entropy_dr2)/(r1 - r2) # note slope is negative
# Recall r1 > r2, since rr decreases with index, and ir2 > ir1
# Also, since all the entropy profiles are monotonically decreasing with 
# radius (increasing with index), we have 
# entropy_dr1 < threshold < entropy_dr2

r_tbl = r2 - (entropy_dr2 - threshold_gradient)/m

dirname_stripped = dirname.split('/')[-1]

delta_T = (ro - r_tbl)/1.e8

print ('r1 (Mm): %.2f, r2 (Mm): %.2f' %(r1/1.e8, r2/1.e8)) 
print ('r_tbl (Mm): %.2f' %(r_tbl/1.e8))
print ('s_dr1 : %.2e, s_dr2 (erg/g/K/cm): %.2e'
        %(entropy_dr1, entropy_dr2)) 
print ('s_threshold : %.2e' %threshold_gradient)
print ('min(s_dr): %.2e' %min_entropy_dr)
print ('m: %.2e' %m)
print ('-----------------------------------')
print ('for the run %s, delta_T = %.3f Mm' %(dirname_stripped, delta_T))
np.save(datadir + 'delta_T.npy', delta_T)
