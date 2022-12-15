# Plot location of Chebyshev, Legendre, and Fourier collocation points
# Can also plot radial/theta/phi locations to see how many points 
# are in various intervals
import matplotlib.pyplot as plt
import numpy as np
import sys, os
sys.path.append(os.environ['rapp'])
sys.path.append(os.environ['raco'])
from grid_util import compute_grid_info
from common import *
from cla_util import *

# Get CLAs
args = sys.argv
clas0, clas = read_clas(args)
dirname = clas0['dirname']
dirname_stripped = strip_dirname(dirname)

# get the current grid info and plot that unless user tells me not to
ncheby, domain_bounds = get_domain_bounds(dirname)

# allowed args + defaults
kwargs_default = dict({'ncheby': ncheby, 'dombounds': domain_bounds, 'nt': get_parameter(dirname, 'n_theta'), 'r0': domain_bounds[0], 'rvals': None, 'latvals': None, 'lonvals': None, 'rminmax': None, 'latminmax': None, 'lonminmax': None, 'rnorm': None})

# overwrite defaults
kw = update_dict(kwargs_default, clas)

# calculate problemsize grid
print (buff_line)
print ("plotting grid for:")
print ("nt = nphi/2 =", kw.nt)
print ("ncheby =", kw.ncheby)
print ("domain_bounds =", kw.dombounds)
nr, nt, nphi, rr, rw, tt, cost, sint, tw, phi, dphi =\
        compute_grid_info(kw.dombounds, kw.ncheby, kw.nt)

# adjust which portion of grid to look at (minmax kwargs)

# radial
if not kw.rminmax is None:
    irmin = np.argmin(np.abs(rr - kw.rminmax[0]))
    irmax = np.argmin(np.abs(rr - kw.rminmax[1]))
    rr = rr[irmax:irmin+1]

# latitude
tt_lat = 180./np.pi*(np.pi/2. - tt)
if not kw.latminmax is None:
    ilatmin = np.argmin(np.abs(tt_lat - kw.latminmax[0]))
    ilatmax = np.argmin(np.abs(tt_lat - kw.latminmax[1]))
    tt_lat = tt_lat[ilatmin:ilatmax+1]

# longitude
phi_lon = 180./np.pi*phi
if not kw.lonminmax is None:
    ilonmin = np.argmin(np.abs(tt_lon - kw.lonminmax[0]))
    ilonmax = np.argmin(np.abs(tt_lon - kw.lonminmax[1]))
    phi_lon = phi_lon[ilonmin:ilonmax+1]

# calculate grid spacing

# radial
dr = np.zeros(nr)
dr[:-1] = rr[:-1] - rr[1:]
dr[-1] = dr[-2]

# theta
dt = np.zeros(nt)
dt[:-1] = tt[:-1] - tt[1:]
dt[-1] = dt[-2]

# generate the plot (1 row of 3 subplots)
fig, axs = plt.subplots(1, 3, sharey=True, figsize=(12,4))

# plot delta_r
size = 0.5
plt.sca(axs[0])
plt.scatter(rr, dr, color='k', s=size)

# make the "rval" location red and big
ir0 = np.argmin(np.abs(rr -r0))
plt.scatter(kw.r0, dr[ir0], color='r', s=3*size)

# set xy axes properties
rmin, rmax = np.min(rr), np.max(rr)
H = rmax - rmin
plt.xlim(rmin - 0.01*H, rmax + 0.01*H)
plt.xlabel('r (cm)')
plt.ylabel(r'$\delta r/(r_{\rm{max}}-r_{\rm{min}})$')

# plot r0 * delta_theta
plt.sca(axs[1])
tt_lat = 180./np.pi*(np.pi/2. - tt)
plt.scatter(tt_lat, dt*kw.r0, color='k', s=size)

# set xy axes properties
plt.xlim(-90., 90.)
plt.xlabel('latitude (deg)')
plt.ylabel(r'$r_0\delta\theta$')

# plot r0 * sin(theta) * dphi
dphi = 2*np.pi/(2*nt)
axs[2].scatter(tt_lat, kw.r0*sint*dphi, color='k', s=size)

# set xy axes properties
plt.sca(axs[2])
plt.xlim(-90., 90.)
plt.xlabel('latitude (deg)')
plt.ylabel(r'$r_0\sin\theta\delta\phi$')

# get ticks everywhere
for ax in axs.flatten():
    plt.sca(ax)
    plt.minorticks_on()
    plt.tick_params(top=True, right=True, direction='in', which='both')


# make title(s)
axs[0].set_title(dirname_stripped + '\n' + 'Grid resolution ')
axs[1].set_title((r'$r_0=%1.3e\ cm$' %kw.r0) + '\n' +\
        (r'$\delta r_0=%1.3e$' %dr[ir0]))
axs[2].set_title((r'$r_0{\rm{max}}(\delta\theta)/{\rm{rnorm}}=%1.3e$' %(kw.r0*np.max(dt))) + '\n' +\
        (r'$r_0{\rm{max}}(\sin\theta\delta\phi)/{\rm{rnorm}}=%1.3e$' %(kw.r0*np.max(sint)*dphi)))

# Display the plot
plt.show()
