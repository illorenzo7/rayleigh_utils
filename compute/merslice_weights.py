import numpy as np
import sys
import os
from diagnostic_reading import Meridional_Slice

dirname = sys.argv[1]
radatadir = dirname + '/Meridional_Slices/'

datadir = dirname + '/data/'
if (not os.path.isdir(datadir)):
    os.makedirs(datadir)

files = os.listdir(radatadir)
nfiles = len(files)
files.sort()

mer0 = Meridional_Slice(radatadir + files[-1], '')

nt = mer0.ntheta
nr = mer0.nr

cost = mer0.costheta
sint = mer0.sintheta
tt = np.arccos(cost)

rr = mer0.radius
ro = rr[0]
ri = rr[nr-1]

tt_prime = (tt[:nt-1] + tt[1:])/2.
rr_prime = (rr[:nr-1] + rr[1:])/2.

# Compute coordinate ranges associated with each collocation point
rr_a = np.zeros(nr); rr_b = np.zeros(nr)
tt_a = np.zeros(nt); tt_b = np.zeros(nt)

rr_a[0], rr_b[0] = rr_prime[0], ro
rr_a[1:nr-1], rr_b[1:nr-1] = rr_prime[1:nr-1], rr_prime[:nr-2]
rr_a[nr-1], rr_b[nr-1] = ri, rr_prime[nr-2]

tt_a[0], tt_b[0] = tt_prime[0], np.pi
tt_a[1:nt-1], tt_b[1:nt-1] = tt_prime[1:nt-1], tt_prime[:nt-2]
tt_a[nt-1], tt_b[nt-1] = 0., tt_prime[nt-2]


# Now compute volume weights, integrating over azimuth
# (Delta_phi --> 2*pi)
tt_factor = np.cos(tt_a) - np.cos(tt_b)
rr_factor = (1./3.)*(rr_b**3 - rr_a**3)

# Volume element
Delta_V = 2*np.pi*tt_factor.reshape((nt,1))*rr_factor.reshape((1,nr))

# Volume weight
Total_V = (4.*np.pi/3.)*(ro**3 - ri**3)
w_v = Delta_V/Total_V

# Radial area element
Delta_S_r = (2*np.pi)*tt_factor.reshape((nt, 1))*(rr**2).reshape((1,nr))

# Radial area weight
Total_A_r = (4.*np.pi*rr**2)
w_ar = Delta_S_r/Total_A_r.reshape((1,nr))

# Latitudinal area element
Delta_S_t = (2.*np.pi)*sint.reshape((nt,1))*\
        0.5*(rr_b**2 - rr_a**2).reshape((1,nr))

# Latitudinal area weight
Total_A_t = np.pi*sint*(ro**2 - ri**2)
w_at = Delta_S_t/Total_A_t.reshape((nt,1))

# Azimuthal area element
Delta_S_p = (tt_b - tt_a).reshape((nt,1))*\
        0.5*(rr_b**2 - rr_a**2).reshape((1,nr))

# Azimuthal area weight
Total_A_p = (np.pi/2)*(ro**2 - ri**2)
w_ap = Delta_S_p/Total_A_p

# Save everything
np.save(datadir + 'merslice_weights.npy',\
        (Delta_V, Delta_S_r, Delta_S_t, Delta_S_p,\
        w_v, w_ar, w_at, w_ap))
np.save(datadir + 'total_A_and_V.npy',\
        (Total_A_r, Total_A_t, Total_A_p, Total_V))
