###############################################
# Author: Loren Matilsky
# Date created: 02/28/2020
#
# This script computes the "cone"-averages latitudinal energy flux in the 
# spherical domain as a function of radius. 
# Plots various average energy fluxes, integrated over cones of opening angle
# theta
# Since Rayleigh has no "cone"-averages, we assume symmetry in phi and use the
# AZ_Avgs data
# the "-2" means we break up fluxes into 
# "conv", "circ", "cond", "visc", "eq", "tot"
# instead of 
# "enth", "enth_pp", "enth_mm", "ke", "cond", "visc", "eq", "tot"

import numpy as np
import sys, os
sys.path.append(os.environ['rapp'])
sys.path.append(os.environ['raco'])
from common import get_widest_range_file, strip_dirname, get_dict
from rayleigh_diagnostics import GridInfo, ReferenceState

# Get the run directory on which to perform the analysis
dirname = sys.argv[1]
dirname_stripped = strip_dirname(dirname)

# Directory with data and plots, make the plotting directory if it doesn't
# already exist    
datadir = dirname + '/data/'
plotdir = dirname + '/plots/'
if (not os.path.isdir(plotdir)):
    os.makedirs(plotdir)

# Find the Shell_Avgs file(s) in the data directory. If there are multiple, by
# default choose the one with widest range in the average
AZ_Avgs_file = get_widest_range_file(datadir, 'AZ_Avgs')

# Get command-line arguments to adjust the interval of averaging files
args = sys.argv[2:]
nargs = len(args)
for i in range(nargs):
    arg = args[i]
    if arg == '-usefile':
        AZ_Avgs_file = args[i+1]
        AZ_Avgs_file = AZ_Avgs_file.split('/')[-1]
       
print ('Getting data from ', datadir + AZ_Avgs_file, ' ...')
di = get_dict(datadir + AZ_Avgs_file)

# Get grid info
gi = GridInfo(dirname + '/grid_info')
nr, nt = gi.nr, gi.ntheta
rr, tt, cost, sint = gi.radius, gi.theta, gi.costheta, gi.sintheta
ri, ro = np.min(rr), np.max(rr)
shell_depth = ro - ri
rweights = gi.rweights

#Create the plot
lw = 1.5 # Bit thicker lines

# Read in the flux data
#print ('Getting AZ_Avgs data from %s ...' %AZ_Avgs_file)
vals = di['vals']
lut = di['lut']
qindex_eflux = lut[1456]
qindex_eflux_pp = lut[1459]
qindex_kflux = lut[1924]

# Get the fluxes in the whole meridional plane
eflux = vals[:, :, qindex_eflux]
kflux = vals[:, :, qindex_kflux]

try:
    eflux_pp = vals[:, :, qindex_eflux_pp]
except:
    # the fluctuating enthalpy flux wasn't computed directly
    # so compute it by subtracting out the flux due to meridional 
    # circulation
    ref = ReferenceState(dirname + '/reference')
    rho_bar = ref.density
    T_bar = ref.temperature
    S_azav = vals[:, :, lut[501]]
    P_azav = vals[:, :, lut[502]]
    vt_azav = vals[:, :, lut[2]]
    eflux_mm = (rho_bar*T_bar*S_azav + P_azav)*vt_azav
    eflux_pp = eflux - eflux_mm

# Compute eflux_mm if the "try" succeeded
eflux_mm = eflux - eflux_pp

# Compute the integrated fluxes
# At each point in the meridional plane we associate a "ring" of width dr and circumference 2 pi r sin(theta)
dr = rweights/rr**2/np.sum(rweights/rr**2)*shell_depth
areas = 2*np.pi*sint.reshape((nt, 1))*rr.reshape((1, nr))*\
        dr.reshape((1, nr))

eflux_int = np.sum(eflux*areas, axis=1)
eflux_pp_int = np.sum(eflux_pp*areas, axis=1)
eflux_mm_int = np.sum(eflux_mm*areas, axis=1)

kflux_int = np.sum(kflux*areas, axis=1)

solar_lum = 3.846e33 # We'll normalize by the solar luminosity

conv_flux = (eflux_pp_int + kflux_int)/solar_lum

print ("Max conv flux / lsun: %.4f" %np.max(conv_flux))
