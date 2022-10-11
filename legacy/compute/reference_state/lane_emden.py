import numpy as np

# define functions used in integrate
def f(x, y, z):
    return z

def g(x, y, z, n):
    return -y**n-2/x*z

# integrate
def rk_integrate(x_0, y_0, z_0, n, t_step):
    '''Integrates one time step'''
    k_0 = t_step * f(x_0, y_0, z_0)
    l_0 = t_step * g(x_0, y_0, z_0, n)
    k_1 = t_step * f(x_0+1/2*t_step, y_0+1/2*k_0, z_0+1/2*l_0)
    l_1 = t_step * g(x_0+1/2*t_step, y_0+1/2*k_0, z_0+1/2*l_0, n)
    k_2 = t_step * f(x_0+1/2*t_step, y_0+1/2*k_1, z_0+1/2*l_1)
    l_2 = t_step * g(x_0+1/2*t_step, y_0+1/2*k_1, z_0+1/2*l_1, n)
    k_3 = t_step * f(x_0+t_step, y_0+k_2, z_0+l_2)
    l_3 = t_step * g(x_0+t_step, y_0+k_2, z_0+l_2, n)
    x_1 = x_0 + t_step
    y_1 = y_0 + 1/6 * (k_0+2*k_1+2*k_2+k_3)
    z_1 = z_0 + 1/6 * (l_0+2*l_1+2*l_2+l_3)
    return (x_1, y_1, z_1)

def lane_emden(n, npoints):
    # define 
    t_step = 1e-4 # integrate super fine to find where x ends
    '''Apply runge-kutta integration to the lane-emden equation'''

    # define initial conditions
    x_0 = 1e-50 # need some small value that is close to 0 but not exactly 0
    y_0 = 1.0
    z_0 = 0

    # do the integration and compile the results into two lists
    while y_0>0:
        x_0, y_0, z_0 = rk_integrate(x_0, y_0, z_0, n, t_step)
        if type(y_0) == complex: # complex values of y_0 sometimes occur when n is not an integer
            break

    # now x_0 = xmax
    t_step = x_0/(npoints - 1)

    # repeat to get the correct number of points
    x_0 = 1e-50 # need some small value that is close to 0 but not exactly 0
    y_0 = 1.0
    z_0 = 0

    # do the integration and compile the results into two lists
    xs = [x_0]
    ys = [y_0]
    while y_0>0:
        x_0, y_0, z_0 = rk_integrate(x_0, y_0, z_0, n, t_step)
        if type(y_0) == complex: # complex values of y_0 sometimes occur when n is not an integer
            break

        xs.append(x_0)
        ys.append(y_0)

    # make arrays
    xs = np.array(xs)
    ys = np.array(ys)

    # rescale the x array
    xs /= np.max(xs)
    return (xs, ys)
