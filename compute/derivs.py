import numpy as np
from scipy.interpolate import UnivariateSpline
from scipy.signal import convolve, gaussian

# Use cubic/quintic splines to compute derivatives numerically!

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

    # gaussian for smoothing
#    len_conv = 4
#    g = gaussian(20,len_conv)
#    g /= np.sum(g)
    # interpolate in radius for each latitude and then take derivative:
#    deriv = np.zeros_like(arr)
#    for tindex in range(nt):
#        arr_rslice = arr[tindex,:]
#        arr_rslice_ext = np.zeros(nr + len_conv)
#        arr_rslice_ext[:len_conv//2] = arr_rslice[0]
#        arr_rslice_ext[len_conv//2:nr+len_conv//2] = arr_rslice
#        arr_rslice_ext[nr + len_conv//2:] = arr_rslice[nr-1]
        
#        slice_smooth_ext = convolve(arr_rslice_ext,gr,mode='same')
#        slice_smooth = slice_smooth_ext[len_conv//2:nr+len_conv//2]
        
#        spline = UnivariateSpline(radius[r_order],slice_smooth[r_order])
#        deriv[tindex,:] = spline.derivative()(radius)

#    return(deriv)

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
    # interpolate in latitude for each radius and then take derivative
#    deriv = np.zeros_like(arr)

    # Gaussian to convolve with wiggly data
#    g = gaussian(30,10) 
#    g /= np.sum(g) #make Gaussian normalized to numerical convolution

#    for rindex in range(nr):
#        arr_tslice = arr[:,rindex]
#        arr_tslice_ext
        # smooth out slice to get rid of spurious wiggles
#        slice_smooth = convolve(arr_tslice, g,mode='same')


#        spline = UnivariateSpline(theta[t_order],slice_smooth[t_order])
#        deriv[:,rindex] = spline.derivative()(theta)

#    return(deriv)

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