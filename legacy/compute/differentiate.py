"""
Module for performing differentiation

Supported options, default values listed in []:
    --data_reversed, is grid monotonically increasing/decreasing [False]
    --method, what method to use ['fd']

    --order, what finite difference order to use [6]

    --zeros, use chebyshev zeros grid [True]
    --matrix, use chebyshev-matrix method [False]
    --endpoints, use chebyshev zeros grid scaled to include endpoints [False]

    --window, specify windowing function for fourier method [None]

Methods for Non-Uniform grids:
    --method=fd: finite difference; 2, 4, 6 order
    --method=cheb: chebyshev grid, zeros & extrema
    --method=cheb-end: chebyshev zeros grid scaled to include endpoints

Methods for Uniform grids:
    --method=cfd: compact finite difference, 6th order
    --method=fourier: fourier transform

"""

from __future__ import print_function
import numpy as np
from scipy import linalg
#from numpy.fft import fft as FFT
#from numpy.fft import ifft as iFFT
#from numpy.fft import fftfreq
from numpy.fft import fft
from numpy.polynomial import chebyshev
from sys import exit

def derivative(data_in, grid_in, data_reversed=False, method='fd', **kwargs):
    """
    user interface to perform 1D differentiation
    """
    supported = ["fd", "cfd", "fourier", "cheb", "cheb-end"]
    supported_fd_orders = [2, 4, 6]

    method = method.lower()
    if (method not in supported):
        print("\n---ERROR: unknown differentiation method = {}".format(method))
        print("\tsupported methods = {}".format(supported))
        exit()

    passed_order = kwargs.pop("order", 6)
    try:
        order = int(passed_order)
    except:
        order = passed_order
    if (method == "fd" and order not in supported_fd_orders):
        print("\n---ERROR: unknown finite difference order = {}".format(order))
        print("\tsupported orders = {}".format(supported_fd_orders))
        exit()

    if (data_reversed):
        data = data_in[::-1]
        grid = grid_in[::-1]
    else:
        data = data_in
        grid = grid_in

    # Finite Difference, non-uniform grids
    if (method == 'fd'):
        if (order == 2):
            dfdx = centered_2nd_order(grid, data)
        elif (order == 4):
            dfdx = centered_4th_order(grid, data)
        elif (order == 6):
            dfdx = centered_6th_order(grid, data)
        else:
            print("\n---ERROR: should never be here...\n")
            exit()

    # Compact Finite Difference, uniform grid only
    elif (method == 'cfd'):
        dfdx = compact_fd_6(grid, data)

    # Fourier, uniform grid only
    elif (method == 'fourier'):
        print("\n---WARNING: fourier differentiation method is untested\n")
        window = kwargs.pop("window", None)
        dfdx = fourier(grid, data, window=window)

    # Chebyshev
    elif ('cheb' in method):
        endpoints = kwargs.pop("endpoints", False)
        if (method == 'cheb-end'):
            endpoints = True
        zeros = kwargs.pop("zeros", True)
        matrix_method = kwargs.pop("matrix", False)

        # assume it's a chebyshev grid to get domain size, but check to be sure
        domain = chebyshev.get_limits(grid, zeros=zeros, endpoints=endpoints)
        cheby_results = chebyshev.is_grid_chebyshev(grid, domain)
        is_cheby = cheby_results['chebyshev']; cheby_type = cheby_results['type']
        e = [cheby_results['zero-error'], cheby_results['extrema-error'],
             cheby_results['zero-endpoints-error']]
        if (not is_cheby):
            print("\n---ERROR: grid doesn't appear to be chebyshev, error = {}\n".format(e))
            exit()

        if (not matrix_method):
            ddx_chebyshev = chebyshev.derivative
        else:
            ddx_chebyshev = chebyshev.derivative_matrix

        dfdx = ddx_chebyshev(data, grid, domain, zeros=zeros, endpoints=endpoints)

    else:
        print("\n---ERROR: should never be here...\n")
        exit()

    if (data_reversed):
        dfdx = dfdx[::-1]

    return dfdx

def fourier(x, f, window=None):
    """
    Calculate first derivative on a uniform mesh using FFT
    """
    x = np.array(x); f = np.array(f)
    N = len(x); dx = np.mean(x[1:] - x[:-1])

    fk, freq, power, xwin = fft.FFTc(f, x, angular=True, window=window)
    i_freq = complex(0,1)*freq
    dfdx = np.real(fft.iFFTc(i_freq*fk))

    return dfdx

def compact_fd_6(x, f):
    """
    Calculate first derivative on a uniform mesh. Uses 6th order
    compact finite differences for the interior points and 6th order
    one-sided differences for near the boundary. Does a tridiagonal
    solve, which is ~O(N)
    """
    x = np.asarray(x); f = np.asarray(f)
    nx = len(x); dfdx = np.zeros((nx))

    # left/right boundary point
    dfdx[0]    = one_sided_6th(x[:7], f[:7])
    dfdx[nx-1] = one_sided_6th(x[-7:], f[-7:], right_edge=True)

    # first/last interior point
    dfdx[1]    = one_sided_6th(x[1:8], f[1:8])
    dfdx[nx-2] = one_sided_6th(x[-8:-1], f[-8:-1], right_edge=True)

    dx = np.mean(x[1:] - x[:-1])
    a = 14./9.
    b = 1./9.
    alpha = 1./3.

    # alpha*f'(i-1) + f'(i) + alpha*f'(i+1) =
    #       b*(f(i+2) - f(i-2))/4/dx + a*(f(i+1) - f(i-1))/2/dx
    # already did boundaries, but include them in the
    # tridiagonal solve since we have that info available
    u = np.zeros(nx); d = np.zeros(nx); l = np.zeros(nx); rhs = np.zeros(nx)
    u[:] = alpha; u[:3] = 0.0; u[-1] = 0.0
    d[:] = 1.0
    l[:] = alpha; l[-3:] = 0.0; l[0] = 0.0
    A = np.matrix([u,d,l])

    # this is for a uniform grid
    # f(i) = f[2:-2]
    rhs[2:-2:] = 0.25*b/dx*(f[4:] - f[:-4]) + 0.5*a/dx*(f[3:-1] - f[1:-3])

    rhs[0] = dfdx[0]; rhs[1] = dfdx[1]; rhs[-2] = dfdx[nx-2]; rhs[-1] = dfdx[nx-1]
    dfdx = linalg.solve_banded((1,1), A, rhs)

    return dfdx

def centered_2nd_order(x, f):
    """
    Calculate first derivative on a non-uniform mesh. Uses 2nd order
    centered differences for the interior points and 2nd order
    one sided differences on the boundary points.
    """
    x = np.asarray(x); f = np.asarray(f)
    nx = len(x); dfdx = np.zeros((nx))

    # interior points
    for i in range(1,nx-1):
        dx01 = x[i-1] - x[i]; dx10 = x[i+1] - x[i]
        num = dx01**2*f[i+1] + (dx10**2-dx01**2)*f[i] - dx10**2*f[i-1]
        den = dx10*dx01*(dx01 - dx10)
        dfdx[i] = num/den

    # left/right boundary point
    dfdx[0]    = one_sided_2nd(x[:3], f[:3])
    dfdx[nx-1] = one_sided_2nd(x[-3:], f[-3:], right_edge=True)

    return dfdx

def centered_4th_order(x, f):
    """
    Calculate first derivative on a non-uniform mesh. Uses 4th order
    centered differences for the interior points and 4th order
    one sided differences on the boundary points. Uses 5 pt stencil.
    Solves a 4x4 system using linalg.solve of O(N**3) every iteration.
    There are ~N iterations, so the cost is ~N*O(4**3) ~64*O(N)
    """
    x = np.asarray(x); f = np.asarray(f)
    nx = len(x); dfdx = np.zeros((nx))
    # solve the following system for f' using determinants
    # f(-2)- f0 = dxm20*f'+dxm20**2*f''/2+dxm20**3*f'''/6+dxm20**4*f''''/24
    # f(-1)- f0 = dxm10*f'+dxm10**2*f''/2+dxm10**3*f'''/6+dxm10**4*f''''/24
    # f(1) - f0 = dxp10*f'+dxp10**2*f''/2+dxp10**3*f'''/6+dxp10**4*f''''/24
    # f(2) - f0 = dxp20*f'+dxp20**2*f''/2+dxp20**3*f'''/6+dxp20**4*f''''/24
    # dxpj0 = x(j) - x(0)
    N = np.zeros((4,4))
    for i in range(2,nx-2):
        dxm20 = x[i-2] - x[i]; dxm10 = x[i-1] - x[i]
        dxp10 = x[i+1] - x[i]; dxp20 = x[i+2] - x[i]
        delta = 0.25*(abs(dxm20) + abs(dxm10) + abs(dxp10) + abs(dxp20))
        dxm20 /= delta; dxm10 /= delta; dxp10 /= delta; dxp20 /= delta
        N[0,:] = np.array([dxm20, dxm20**2, dxm20**3, dxm20**4])
        N[1,:] = np.array([dxm10, dxm10**2, dxm10**3, dxm10**4])
        N[2,:] = np.array([dxp10, dxp10**2, dxp10**3, dxp10**4])
        N[3,:] = np.array([dxp20, dxp20**2, dxp20**3, dxp20**4])

        # get denominator = determinant of N
        den = linalg.det(N)

        # get numerator = det of A with first column replaced by rhs
        N[:,0] = np.array([f[i-2], f[i-1], f[i+1], f[i+2]]) - f[i]
        num = linalg.det(N)

        # construct 1st derivative, multiply (num/den) by n!/delta**n for nth derivative
        dfdx[i] = (num/den)/delta

    # left/right boundary point
    dfdx[0]    = one_sided_4th(x[:5], f[:5])
    dfdx[nx-1] = one_sided_4th(x[-5:], f[-5:], right_edge=True)

    # first/last interior point
    dfdx[1]    = one_sided_4th(x[1:6], f[1:6])
    dfdx[nx-2] = one_sided_4th(x[-6:-1], f[-6:-1], right_edge=True)

    return dfdx

def centered_6th_order(x, f):
    """
    Calculate first derivative on a non-uniform mesh. Uses 6th order
    centered differences for the interior points and 6th order
    one sided differences on the boundary points. Uses 7 pt stencils
    Solves a 6x6 system using linalg.solve of O(N**3) every iteration.
    There are ~N iterations, so the cost is ~N*O(6**3) ~216*O(N)
    """
    x = np.asarray(x); f = np.asarray(f)
    nx = len(x); dfdx = np.zeros((nx))

    # interior points:
    # solve the following system for f' using determinants
    # f(-3)- f0 = dxm30*f'+dxm30**2*f''/2+dxm30**3*f'''/6+dxm30**4*f''''/24
    # f(-2)- f0 = dxm20*f'+dxm20**2*f''/2+dxm20**3*f'''/6+dxm20**4*f''''/24
    # f(-1)- f0 = dxm10*f'+dxm10**2*f''/2+dxm10**3*f'''/6+dxm10**4*f''''/24
    # f(1) - f0 = dxp10*f'+dxp10**2*f''/2+dxp10**3*f'''/6+dxp10**4*f''''/24
    # f(2) - f0 = dxp20*f'+dxp20**2*f''/2+dxp20**3*f'''/6+dxp20**4*f''''/24
    # f(3) - f0 = dxp30*f'+dxp30**2*f''/2+dxp30**3*f'''/6+dxp30**4*f''''/24
    # dxpj0 = x(j) - x(0); dxmj0 = x(-j) - x(0)
    N = np.zeros((6,6))
    for i in range(3,nx-3):
        dxm30 = x[i-3] - x[i]; dxm20 = x[i-2] - x[i]; dxm10 = x[i-1] - x[i]
        dxp10 = x[i+1] - x[i]; dxp20 = x[i+2] - x[i]; dxp30 = x[i+3] - x[i]
        delta = (abs(dxm30)+abs(dxm20)+abs(dxm10)+abs(dxp10)+abs(dxp20)+abs(dxp30))/6.
        dxm30 /= delta; dxm20 /= delta; dxm10 /= delta
        dxp10 /= delta; dxp20 /= delta; dxp30 /= delta
        N[0,:] = np.array([dxm30, dxm30**2, dxm30**3, dxm30**4, dxm30**5, dxm30**6])
        N[1,:] = np.array([dxm20, dxm20**2, dxm20**3, dxm20**4, dxm20**5, dxm20**6])
        N[2,:] = np.array([dxm10, dxm10**2, dxm10**3, dxm10**4, dxm10**5, dxm10**6])
        N[3,:] = np.array([dxp10, dxp10**2, dxp10**3, dxp10**4, dxp10**5, dxp10**6])
        N[4,:] = np.array([dxp20, dxp20**2, dxp20**3, dxp20**4, dxp20**5, dxp20**6])
        N[5,:] = np.array([dxp30, dxp30**2, dxp30**3, dxp30**4, dxp30**5, dxp30**6])

        # get denominator = determinant of N
        den = linalg.det(N)

        # get numerator = det of A with first column replaced by rhs
        N[:,0] = np.array([f[i-3], f[i-2], f[i-1], f[i+1], f[i+2], f[i+3]]) - f[i]
        num = linalg.det(N)

        # construct 1st derivative, multiply (num/den) by n!/delta**n for nth derivative
        dfdx[i] = (num/den)/delta

    # left/right boundary point
    dfdx[0]    = one_sided_6th(x[:7], f[:7])
    dfdx[nx-1] = one_sided_6th(x[-7:], f[-7:], right_edge=True)

    # first/last interior point
    dfdx[1]    = one_sided_6th(x[1:8], f[1:8])
    dfdx[nx-2] = one_sided_6th(x[-8:-1], f[-8:-1], right_edge=True)

    # second/second to last interior point
    dfdx[2]    = one_sided_6th(x[2:9], f[2:9])
    dfdx[nx-3] = one_sided_6th(x[-9:-2], f[-9:-2], right_edge=True)

    return dfdx

def one_sided_2nd(xx, ff, right_edge=False):
    """
    get 2nd order one sided difference

    input is x = (x0, x1, x2)
    input is f = (f0, f1, f2)

    this assumes x0 is the left edge, if it is the
    right edge, we reverse the incoming (x,f)

    output is dfdx at x0 and is a scalar
    """
    xx = np.asarray(xx); ff = np.asarray(ff)
    if (right_edge):
        xx = xx[::-1]; ff = ff[::-1]

    dx20 = xx[2] - xx[0]; dx10 = xx[1] - xx[0]; dx21 = xx[2] - xx[1]
    num = dx20**2*ff[1] - dx10**2*ff[2] + (dx10**2 - dx20**2)*ff[0]
    den = dx10*dx20*dx21

    return num/den

def one_sided_4th(xx, ff, right_edge=False):
    """
    get 4th order one sided difference

    input is x = (x0, x1, x2, x3, x4)
    input is f = (f0, f1, f2, f3, f4)

    this assumes x0 is the left edge, if it is the
    right edge, we reverse the incoming (x,f)

    output is dfdx at x0 and is a scalar

    single matrix solve of 5x5 which is O(N**3) ~O(125)
    """
    xx = np.asarray(xx); ff = np.asarray(ff)
    if (right_edge):
        xx = xx[::-1]; ff = ff[::-1]
    # generate matrix such that A*f' = f
    # where the rows of A are taylor expansions:
    # f0 = f0
    # f1 = f0 + dx10*f' + dx10**2*f''/2 + dx10**3*f'''/6 + dx10**4*f''''/24
    # f2 = f0 + dx20*f' + dx20**2*f''/2 + dx20**3*f'''/6 + dx20**4*f''''/24
    # f3 = f0 + dx30*f' + dx30**2*f''/2 + dx30**3*f'''/6 + dx30**4*f''''/24
    # f4 = f0 + dx40*f' + dx40**2*f''/2 + dx40**3*f'''/6 + dx40**4*f''''/24
    # rewrite it so dx --> dx/delta for matrix conditioning
    # this makes the unknown vector f^(n)*delta**n/n!

    dx10 = xx[1] - xx[0]; dx20 = xx[2] - xx[0]; dx30 = xx[3] - xx[0]; dx40 = xx[4] - xx[0]
    delta = 0.25*(abs(dx10) + abs(dx20) + abs(dx30) + abs(dx40))
    dx10 /= delta; dx20 /= delta; dx30 /= delta; dx40 /= delta

    N = np.zeros((5,5))
    N[0,:] = np.array([1., 0.,   0.,      0.,      0.     ])
    N[1,:] = np.array([1., dx10, dx10**2, dx10**3, dx10**4])
    N[2,:] = np.array([1., dx20, dx20**2, dx20**3, dx20**4])
    N[3,:] = np.array([1., dx30, dx30**2, dx30**3, dx30**4])
    N[4,:] = np.array([1., dx40, dx40**2, dx40**3, dx40**4])

    dn_f_dxn = linalg.solve(N, ff)

    # fNprime = dn_f_dxn[N] * N! / delta**N where N specifies what derivative
    f1prime = dn_f_dxn[1]/delta  # 1st deriv to O(dx**4)

    return f1prime

def one_sided_6th(xx, ff, right_edge=False):
    """
    get 6th order one sided difference

    input is x = (x0, x1, x2, x3, x4, x5, x6)
    input is f = (f0, f1, f2, f3, f4, x5, x6)

    this assumes x0 is the left edge, if it is the
    right edge, we reverse the incoming (x,f)

    output is dfdx at x0 and is a scalar

    single matrix solve of 7x7 which is O(N**3) ~O(343)
    """
    xx = np.asarray(xx); ff = np.asarray(ff)
    if (right_edge):
        xx = xx[::-1]; ff = ff[::-1]
    # generate matrix such that A*f' = f
    # where the rows of A are taylor expansions:
    # f0 = f0
    # f1 = f0 + dx10*f' + dx10**2*f''/2 + dx10**3*f'''/6 + dx10**4*f''''/24 + ...
    # f2 = f0 + dx20*f' + dx20**2*f''/2 + dx20**3*f'''/6 + dx20**4*f''''/24 + ...
    # f3 = f0 + dx30*f' + dx30**2*f''/2 + dx30**3*f'''/6 + dx30**4*f''''/24 + ...
    # f4 = f0 + dx40*f' + dx40**2*f''/2 + dx40**3*f'''/6 + dx40**4*f''''/24 + ...
    # f5 = f0 + dx50*f' + dx50**2*f''/2 + dx50**3*f'''/6 + dx50**4*f''''/24 + ...
    # f6 = f0 + dx60*f' + dx60**2*f''/2 + dx60**3*f'''/6 + dx60**4*f''''/24 + ...
    # rewrite it so dx --> dx/delta for matrix conditioning
    # this makes the unknown vector f^(n)*delta**n/n!

    dx10 = xx[1] - xx[0]; dx20 = xx[2] - xx[0]; dx30 = xx[3] - xx[0]
    dx40 = xx[4] - xx[0]; dx50 = xx[5] - xx[0]; dx60 = xx[6] - xx[0]
    delta = (abs(dx10) + abs(dx20) + abs(dx30) + abs(dx40) + abs(dx50) + abs(dx60))/6.
    dx10 /= delta; dx20 /= delta; dx30 /= delta; dx40 /= delta; dx50 /= delta; dx60 /= delta

    N = np.zeros((7,7))
    N[0,:] = np.array([1., 0.,   0.,      0.,      0.,      0.,      0.     ])
    N[1,:] = np.array([1., dx10, dx10**2, dx10**3, dx10**4, dx10**5, dx10**6])
    N[2,:] = np.array([1., dx20, dx20**2, dx20**3, dx20**4, dx20**5, dx20**6])
    N[3,:] = np.array([1., dx30, dx30**2, dx30**3, dx30**4, dx30**5, dx30**6])
    N[4,:] = np.array([1., dx40, dx40**2, dx40**3, dx40**4, dx40**5, dx40**6])
    N[5,:] = np.array([1., dx50, dx50**2, dx50**3, dx50**4, dx50**5, dx50**6])
    N[6,:] = np.array([1., dx60, dx60**2, dx60**3, dx60**4, dx60**5, dx60**6])

    dn_f_dxn = linalg.solve(N, ff)

    # fNprime = dn_f_dxn[N] * N! / delta**N where N specifies what derivative
    f1prime = dn_f_dxn[1]/delta     # 1st deriv to O(dx**6)

    return f1prime

