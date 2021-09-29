# module to get spin-down history of a star using the cluster-calibrated
# nonlinear relationship in Barnes et al. (2010)
from scipy.optimize import minimize
import numpy as np

k_I = 452.
k_C = 0.646

def age_from_period(P, P_0, tau): # equation (19) Barnes et al. 2010
    # returns t (stellar age in Myr)
    # must specify:
    # P   =  current rotation rate in days
    # P_0 = initial rotation rate in days
    # tau = convective overturning time in days
    rhs = np.log(P/P_0) + k_I*k_C*(P**2 - P_0**2)/2/tau**2
    return tau/k_C*rhs

def Skumanich_law(t, P_0, tau):
    return np.sqrt(P_0**2 + 2*t*tau/k_I)

def age_period_relation(main_sequence_age, P_0, tau, nages=10000):
    print ('ms age = ', main_sequence_age)
    Skuman = Skumanich_law(main_sequence_age, P_0, tau)
    periods = np.linspace(P_0, Skuman, nages)
    ages = age_from_period(periods, P_0, tau)
    return ages, periods
