# Author: Loren Matilsky
# Created: well before 05/06/2019
# Extremely basic routines to compute very basic finite-difference
# derivatives (central difference). Will work best on evenly-spaced grids.
# For the 2d derivatives (drad and dth), we assume the array is organized 
# as arr(nt, nr) (i.e., theta <--> rows, radius <--> columns)
import numpy as np
from compute_grid_info import compute_grid_info, compute_theta_grid,\
        compute_r_grid

def drad(arr, radius):
    nt, nr = np.shape(arr)
    two_dr = np.zeros((1,nr-2))
    two_dr[0, :] = radius[:nr-2] - radius[2:nr]     
    deriv = np.zeros_like(arr)
    deriv[:,1:nr-1] = (arr[:,:nr-2] - arr[:,2:nr])/two_dr
    deriv[:,0] = deriv[:,1]
    deriv[:,nr-1] = deriv[:,nr-2]
    return deriv

def dth(arr,theta):
    nt, nr = np.shape(arr)
    two_dt = np.zeros((nt-2 ,1))
    two_dt[:, 0] = theta[:nt-2] - theta[2:nt]     
    deriv = np.zeros_like(arr)
    deriv[1:nt-1, :] = (arr[:nt-2,:] - arr[2:nt,:])/two_dt
    deriv[0, :] = deriv[1, :]
    deriv[nt-1, :] = deriv[nt-2, :]
    return deriv

def deriv_1d(x, y):
    n = len(x)
    deriv = np.zeros(n)
    
    two_dx = x[2:] - x[:n-2]
    two_dy = y[2:] - y[:n-2]
    deriv[1:n-1] = two_dy/two_dx
    deriv[0] = deriv[1]
    deriv[n-1] = deriv[n-2]
    return(deriv)
    
def int_1d(x, y):
    n = len(x)
    dx = x[1:] - x[:n-1]
    ymid = 0.5*(y[1:] + y[:n-1])
    return (np.sum(ymid*dx))

def int_r(arr, radius, use_extrema=False):
    int_arr = np.zeros_like(arr)
    nt, nr = np.shape(arr)
    rmin, rmax = np.min(radius), np.max(radius)
    r, rw = compute_r_grid(nr, rmin, rmax, use_extrema=use_extrema)
    rw_2d = rw.reshape((1, nr))

    int_arr[:, nr - 1] = 0.  # Integral from rmin to rmin: = 0
    for ir in range(nr - 1):
        r0 = r[ir]
        dr = (1./3.)*(r0**3 - rmin**3)*rw_2d/r**2/np.sum(rw[ir:])
        int_arr[:, ir] = np.sum((arr*dr)[:, ir:], axis=1)
    return int_arr

def int_th(arr):
    int_arr = np.zeros_like(arr)
    nt, nr = np.shape(arr)
    tt, tw = compute_theta_grid(nt)
    cost, sint = np.cos(tt), np.sin(tt)

    int_arr[nt - 1, :] = 0.  # Integral from theta=0 to theta=0: = 0

    for it in range(nt):
        theta0 = tt[it]
        dt = (1. - cost[it])*tw*np.sin(tt)/np.sum(tw[it:])
        dt_2d = dt.reshape((nt, 1))
        int_arr[it, :] = np.sum((arr*dt_2d)[it:, :], axis=0)

    return int_arr
