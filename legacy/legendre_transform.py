# my attempt to reproduce Nick's Legendre evaluation
import numpy as np
from compute_grid_info import compute_theta_grid

def compute_factorial_ratio(mval):
    # build sqrt( (2m-1)!!(2m-1)!!/(2m)!! ) stably
    ratio = 1.0
    for i in range(1, mval+1):
        ratio = ratio*((i-0.5)/i)**0.5  # ratio = ratio*(2m-1)/(2m)
    return ratio

def compute_plm(nt):
    # calculates P_lm(x) at the Legendre collocation points 
    # (dealias so l_max = 3 * ntheta // 2)
    # all in all there will be ~ (1/2) (2/3*ntheta)**2 * ntheta
    # ~ 0.22 ntheta**3 nonzero values in the ~ 0.44 * ntheta**3 array
    # plm[0:lmax, 0:mmax, 0:ntheta]
    # with mmax = lmax 
    # and only nonzero for im <= il

    # set up theta grid
    tt, tw = compute_theta_grid(nt)
    cost = np.cos(tt)

    # compute nell, lmax
    nell = 2*nt//3
    lmax = nell-1

    # output array
    plm = np.zeros((nell, nell, nt))

    # loop over mvals
    for mval in range(nell):
        # First, fill in the l=mval, mval+1 pieces (closed form expression)
        factorial_ratio = compute_factorial_ratio(mval)
        amp = ((mval+0.5)/(2*np.pi))**0.5
        amp *= factorial_ratio

        for it in range(nt):
            x = cost[it]

            # get the l = mval polynomial (analytical)
            # P_(l, l)(x) = (2l - 1)!! sqrt(1 - x**2)^(l/2)
            lval = mval
            plm[lval, lval, it] = (-1)**(lval)*amp*(1 - x**2)**(lval/2)

            # and the l = mval + 1 (if applicable)
            # P_(l+1, l)(x) = sqrt(2*l + 1) * x * P_(l, l)(x)
            if mval < lmax:
                plm[lval+1, lval, it] =\
                    plm[lval, lval, it]*x*(2*lval+1)**0.5

        # now use the general recursion for l > mval+1
        for lval in range(mval + 2, lmax + 1):
            for it in range(nt):
                x = cost[it]
                amp = (lval-1)**2 - mval**2
                amp /= (4*(lval-1)**2 - 1)
                amp = amp**0.5
                tmp = plm[lval-1, mval, it]*x - amp*plm[lval-2, mval, it]
                amp = (4*lval**2 - 1)/(lval**2 - mval**2)
                plm[lval, mval, it] = tmp*amp**0.5

    return plm
