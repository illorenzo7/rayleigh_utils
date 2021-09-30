# module to get spin-down history of a star using the cluster-calibrated
# nonlinear relationship in Barnes et al. (2010)
from scipy.optimize import minimize
import numpy as np

tau_sun = 35. # convective overturning time of the Sun
age_sun = 4.6e3 # age of Sun in Myr
P_sun = 25.38 # rotation period of the Sun (sidereal Carrington)
# k_I = 452.
# k_C = 0.646
k_I = 2*tau_sun*age_sun/P_sun**2
k_C = tau_sun/54. # make efolding in "C" sequence 54 Myr for Sun

def age_from_period(P, P0, tau, t0=0.): # equation (19) Barnes et al. 2010
    # returns t (stellar age in Myr)
    # must specify:
    # P   =  current rotation rate in days
    # P0 = initial rotation rate in days
    # tau = convective overturning time in days
    rhs = np.log(P/P0) + k_I*k_C*(P**2 - P0**2)/2/tau**2
    return tau/k_C*rhs + t0

def P0_from_age(t, P, tau, n=1000):
    Pmin = 0.05 # no stars are observed to rotate faster than this!
    period_sample = np.linspace(Pmin, P, n)
    ages_sample = age_from_period(P, period_sample, tau)
    i0 = np.argmin(np.abs(ages_sample - t))
    i1 = i0
    if i1 == n - 1:
        P0 = P
    else:
        i2 = i0+1
        age1 = ages_sample[i1]
        age2 = ages_sample[i2]
        P1 = period_sample[i1]
        P2 = period_sample[i2]
        P0 = P1 + (t - age1)*(P2 - P1)/(age2 - age1)
    return P0

def skumanich(t, P0, tau, t0=0.):
    return np.sqrt(P0**2 + 2*(t-t0)*tau/k_I)

def age_period_relation(P0, tau, t0=0., max_age=1.37e4, nages=1000): 
    # calculate by default for max age = 13.7 B yr (age of universe)
    # P0 is the rotation rate at (by default) t = t0 = 0 Myr

    Pmin = P0_from_age(t0, P0, tau)
    Pmax = skumanich(max_age, P0, tau, t0=t0)
    periods = np.linspace(Pmin, Pmax, nages)
    ages = age_from_period(periods, P0, tau, t0=t0)
    return ages, periods

#def period_from_age(t, P0, tau)
