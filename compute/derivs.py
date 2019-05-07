# Author: Loren Matilsky
# Created: well before 05/06/2019
# Extremely basic routines to compute very basic finite-difference
# derivatives (central difference). Will work best on evenly-spaced grids.
# For the 2d derivatives (drad and dth), we assume the array is organized 
# as arr(nt, nr) (i.e., theta <--> rows, radius <--> columns)
import numpy as np

def drad(arr,radius):
    nt,nr = np.shape(arr)
    r_order = np.argsort(radius)
    two_dr = np.zeros((1,nr-2))
    two_dr[0,:] = radius[:nr-2] - radius[2:nr]     
    deriv = np.zeros_like(arr)
    deriv[:,1:nr-1] = (arr[:,:nr-2] - arr[:,2:nr])/two_dr
    deriv[:,0] = deriv[:,1]
    deriv[:,nr-1] = deriv[:,nr-2]
    return(deriv)

def dth(arr,theta):
    nt,nr = np.shape(arr)
    t_order = np.argsort(theta)
    two_dt = np.zeros((nt-2,1))
    two_dt[:,0] = theta[:nt-2] - theta[2:nt]     
    deriv = np.zeros_like(arr)
    deriv[1:nt-1,:] = (arr[:nt-2,:] - arr[2:nt,:])/two_dt
    deriv[0,:] = deriv[1,:]
    deriv[nt-1,:] = deriv[nt-2,:]
    return(deriv)

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
