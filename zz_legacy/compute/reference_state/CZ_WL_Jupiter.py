# Author: Loren Matilsky
# Created: 05/01/2019
# Modified to work with new cref framework: 10/20/2019
# New CLA interface: 08/30/2021
#
# Purpose: generate a binary file (for Rayleigh to read) that contains
# a reference state, consisting of a stable region (stiffness k) underlying 
# a neutrally stable CZ
# assumes g(r) = -G * mstar / r^2

# Parameters: output_dir (first argument), 
# Command-line options:
#
# --rmin, --rmax, --rt
# Inner, outer, and transition radii (default rmin = 0.1*rstar
# (i.e. very deep shell; don't need to use all of it)
# rt = 5e10 (solar convection zone in mind)
# rmax = 0.99*rsun (as close to surface as my simple
# hydro model allows)

# --mstar
# central mass, default mstar = msun
# 
# --rstar 
# stellar radius, default rstar = rsun
#
# --rhot
# Density at transition boundary, default sun.rho_bcz ~ 0.18
#
# --tmpt
# Temperature at transition boundary default sun.tmp_bcz ~ 2.1e6
#
# --k
# Stiffness k, default 2.0 (define in units of dsdr in RZ / (cp/rstar) for consistency
# 
# --delta
# Transition width delta, default 0.05*rsun
#
# --gamma
# specific heat ratio gamma, default gamma_ideal = 5/3
#
# --cp
# pressure specific heat cp, default sun.c_p = 3.5e8
#
# --mag
# Whether magnetism is True or False, default False (hydro)
#
# --fname
# File to save reference state in (default custom_reference_binary)
#
# --nr
# Default number of radial (evenly spaced) grid points. Default 10,000 (very fine)

def integrate_from_r0(integrand, rr, r0):
    # basic grid info
    nr = len(rr)
    ir0 = np.argmin(np.abs(rr - r0))
   
    # compute indefinite integral
    integral = np.zeros(nr)
    for ir in range(nr):
        integral[ir] = simps(integrand[ir:ir0 + 1], rr[ir:ir0 + 1])
        if ir <= ir0:
            integral[ir] *= -1
    return integral

def rho_and_t(nr=5000, alpha=0.25, beta=0.9, gam=5./3., delta=0.1, nrho=3.):

    # compute radial locations and grid
    rin = beta/(1.-beta)
    r0 = 1./(1.-beta)
    rout = r0 + alpha
    rr = np.linspace(rin, rout, nr)

    # compute gravity
    grav = (1.-beta**3)/(3*(1.-beta)**3) / rr**2

    # compute masking function
    psi_WL = np.zeros(nr)
    for ir in range(nr):
        rr_loc = rr[ir]
        if rr_loc >= r0 + delta:
            psi_WL[ir] = 1.
        elif rr_loc > r0:
            psi_WL[ir] = 1. - (1. - ((rr_loc-r0)/delta)**2 )**2

    # get entropy gradient (in this case, just masking function)
    dsdr = psi_WL

    # compute shape of N^2(r)
    nsq = grav*dsdr

    # rescale so that N^2 has volume-average over WL of 1
    fourpi = 4.*np.pi
    vol_WL = (fourpi/3.)*(rout**3 - r0**3)
    nsq_volav = simps(nsq*fourpi*rr**2, rr)/vol_WL
    nsq /= nsq_volav

    # get entropy 
    s = integrate_from_r0(dsdr, rr, r0)

    # get dissipation number
    npoly = 1./(1.-gamma)
    expp = np.exp(nrho/npoly)
    expm = np.exp(-nrho/npoly)

    numer = 3.*beta*(1.-beta)**2*(1.-expm)
    denom = (3.*beta/2.)*(1.-beta**2)*(1.-expm) -\
            (1-beta**3)*(beta - expm)
    diss = numer/denom

    # get temp. at r0
    numer = (1.-beta**3)*(1.-beta)
    denom = (3.*beta/2.)*(1.-beta**2)*(expp-1.) -\
            (1-beta**3)*(beta*expp - 1.)
    tmp0 = numer/denom

    # now get temperature
    integrand = grav*np.exp(-s)
    integral = integrate_from_r0(integrand, rr, r0)
    tmp = np.exp(s)*(tmp0 - diss*integral)

    # and finally get density, normalize to unity over CZ
    rho = np.exp(-(gamma/(gamma-1.))*s)*tmp**npoly
    vol_CZ = (fourpi/3.)*(r0**3 - rin**3)
    rho_volav = simps((rho*fourpi*rr**2)[:ir0+1], rr[:ir0+1])/vol_CZ
    rho /= rho_volav

    # and finally get derivatives

    # first derivatives
    dtdr = dsdr*tmp - diss*grav
    dlnt = dtdr/tmp
    dlnrho = (1.0/(gamma - 1.0))*(dlnT - gamma*dsdr)

    # second derivatives
    d2sdr2 = np.gradient(dsdr, rr)
    dgdr = -2.*grav/rr
    
    d2lnT = d2sdr2 + (diss/tmp)*(grav*dlnt - dgdr)
    d2lnrho = (1.0/(gamma - 1.0))*(d2lnT - gamma*d2sdr2)
    
    return rho, tmp, grav, nsq, dlnt, dlnrho, d2lnrho

import numpy as np
import sys, os
from arbitrary_atmosphere import arbitrary_atmosphere_nd

sys.path.append(os.environ['rapp'])
sys.path.append(os.environ['raco'])

from reference_tools import equation_coefficients
from common import *
from cla_util import *

# Get CLAs
args = sys.argv
clas0, clas = read_clas_raw(args)
dirname = clas0['dirname']

# Set default kwargs
kw_default = dotdict(dict({'rmin': 0.1*sun.r, 'rmax': 0.99*sun.r, 'rt': sun.r_bcz, 'mstar': msun, 'rstar': rsun, 'rhot': sun.rho_bcz, 'tmpt': sun.tmp_bcz, 'k': 2.0, 'delta': 0.05*sun.r, 'gamma': gamma_ideal, 'cp': sun.c_p, 'mag': False, 'fname': 'custom_reference_binary', 'nr': 10000}))

# overwrite defaults
kw = update_dict(kw_default, clas)

# check for bad keys
find_bad_keys(kw_default, clas, clas0['routinename'], justwarn=True)

# compute gas constant R
gas_constant = kw.cp*(kw.gamma-1)/kw.gamma
prst = kw.rhot*gas_constant*kw.tmpt

# compute reference state on super-fine grid to interpolate onto later    
rr = np.linspace(kw.rmax, kw.rmin, kw.nr) # keep radius in decreasing order for consistency with Rayleigh convention

# Define an entropy profile that is +1 for r < rt - delta, 0 for r > rt,
# and continuously differentiable (a quartic) in between

s = np.zeros(kw.nr)
dsdr = np.zeros(kw.nr)
d2sdr2 = np.zeros(kw.nr)

for i in range(kw.nr):
    rloc = rr[i]
    if rloc <= kw.rt - kw.delta:
        s[i] = (8.0/15.0)*(kw.k*kw.cp)*(kw.delta/kw.rstar) +\
                kw.k*kw.cp*(rloc/kw.rstar - kw.rt/kw.rstar)
        dsdr[i] = kw.k*kw.cp/kw.rstar
        d2sdr2[i] = 0.0
    elif rloc > kw.rt - kw.delta and rloc < kw.rt:
        x = (rloc - kw.rt)/kw.delta
        s[i] = (kw.k*kw.cp)*(kw.delta/kw.rstar)*((2.0/3.0)*x**3.0 - (1.0/5.0)*x**5.0)
        dsdr[i] = (kw.k*kw.cp/kw.rstar)*(1.0 - (1.0 - x**2.0)**2.0)
        d2sdr2[i] = (4.0/kw.delta)*(kw.k*kw.cp/kw.rstar)*(1.0 - x**2.0)*x
    else:
        s[i] = 0.0
        dsdr[i] = 0.0
        d2sdr2[i] = 0.0

# Make gravity due to a point mass at the origin
g = g_univ*kw.mstar/rr**2.0
dgdr = -2.0*g/rr

T, rho, p, dlnT, dlnrho, dlnp, d2lnrho =\
    arbitrary_atmosphere(rr, s, dsdr, d2sdr2, g,\
                         dgdr, kw.rt, kw.tmpt, prst, kw.cp, kw.gamma)

print(buff_line)
print("Computed atmosphere for RZ-CZ, ds/dr joined with a quartic")
print("nr         : %i" %kw.nr) 
print("rstar      : %1.8e cm" %kw.rstar) 
print("mstar      : %1.8e  g" %kw.mstar) 
print("rmin/rstar : %.8f    " %(kw.rmin/kw.rstar))
print("rt/rstar   : %.8f    " %(kw.rt/kw.rstar))
print("rmax/rstar : %.8f    " %(kw.rmax/kw.rstar))
print("delta/rstar: %.8f"  %(kw.delta/kw.rstar))
print("stiffness k: %1.8e"  %kw.k)
print("stiffness  : %1.8e"  %(kw.k*kw.cp/kw.rstar))
print("...where stiffness := (dsdr in RZ)")
print("and k = (stiffness)/(cp/rstar)")
print(buff_line)

# Now write to file using the equation_coefficients framework
eq = equation_coefficients(rr)

# Set only the thermodynamic functions/constants in this routine
# In other routines, we can set the heating and transport coefficients
# Only set c_4 = 1/(4*pi) if mag = True

print("Setting f_1, f_2, f_4, f_8, f_9, f_10, and f_14")
eq.set_function(rho, 1)
buoy = rho*g/kw.cp
eq.set_function(buoy, 2)
eq.set_function(T, 4)
eq.set_function(dlnrho, 8)
eq.set_function(d2lnrho, 9)
eq.set_function(dlnT, 10)
eq.set_function(dsdr, 14)

print("Setting c_2, c_3, c_7, and c_8")
eq.set_constant(1.0, 2) # multiplies buoyancy
eq.set_constant(1.0, 3) # multiplies pressure grad.
eq.set_constant(1.0, 8) # multiplies viscous heating

if kw.mag:
    print("magnetism = True, so setting c_4 = c_9 =  1/(4*pi), c_7 = 1")
    eq.set_constant(1.0/4.0/np.pi, 4) # multiplies Lorentz force
    eq.set_constant(1.0, 7) # multiplies eta in induction-diffusion term
    eq.set_constant(1.0/4.0/np.pi, 9) # multiplies Joule heating
else:
    print("magnetism = False, so setting c_4, c_7, c_9 = 0.0")
    eq.set_constant(0.0, 4) # multiplies Lorentz force
    eq.set_constant(0.0, 7) # multiplies eta in induction-diffusion term
    eq.set_constant(0.0, 9) # multiplies Ohmic heating

# c_10 will be set in the "generate_heating" scripts

# The "generate_transport" scripts will set the transport
# "radial shapes", and the constants c_5, c_6, c_7

the_file = dirname + '/' + kw.fname

print("Writing the atmosphere to %s" %the_file)
print("---------------------------------")
eq.write(the_file)



