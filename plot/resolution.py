# Plot location of Chebyshev, Legendre, and Fourier collocation points
# Can also plot radial/theta/phi locations to see how many points 
# are in various intervals
import matplotlib.pyplot as plt
import numpy as np
import sys, os
sys.path.append(os.environ['rapp'])
sys.path.append(os.environ['raco'])
from compute_grid_info import compute_grid_info
from common import *

# Get CLAs
args = sys.argv
clas0, clas = read_clas(args)
dirname = clas0['dirname']
dirname_stripped = strip_dirname(dirname)

# get the current grid info and plot that unless user tells me not to
ncheby, domain_bounds = get_domain_bounds(dirname)

# allowed args + defaults
kwargs_default = dict({'ncheby': ncheby, 'rbounds': domain_bounds, 'nt': get_parameter(dirname, 'n_theta'), 'rvals': None, 'latvals': None, 'lonvals': None, 'rminmax': None, 'latmin': None, 'lonminmax': None, 'rnorm': None})

# overwrite defaults
kw = update_dict(kwargs_default, clas)

# calculate problemsize grid
print (buffline)
print ("plotting grid for:")
print ("nt = nphi/2 =", kw.nt)
print ("ncheby =", kw.ncheby)
print ("rbounds (domain_bounds) =", kw.rbounds)
nr, nt, nphi, rr, rw, tt, cost, sint, tw, phi, dphi =\
        compute_grid_info(kw.domain_bounds, kw.ncheby, kw.nt)

# adjust which portion of grid to look at (minmax kwargs)

# radial
if not kw.rminmax is None:
    irmin = np.argmin(np.abs(rr - kw.rminmax[0]))
    irmax = np.argmin(np.abs(rr - kw.rminmax[1]))
    rr = rr[irmax:irmin+1]
rmin, rmax = rr[-1], rr[0]

# latitude
tt_lat = 180./np.pi*(np.pi/2. - tt)
if not kw.latminmax is None:
    ilatmin = np.argmin(np.abs(tt_lat - kw.latminmax[0]))
    ilatmax = np.argmin(np.abs(tt_lat - kw.latminmax[1]))
    tt_lat = tt_lat[ilatmin:ilatmax+1]
latmin, latmax = tt_lat[0], tt_lat[-1]

# longitude
phi_lon = 180./np.pi*phi
if not kw.lonminmax is None:
    ilonmin = np.argmin(np.abs(tt_lon - kw.lonminmax[0]))
    ilonmax = np.argmin(np.abs(tt_lon - kw.lonminmax[1]))
    phi_lon = phi_lon[ilonmin:ilonmax+1]
lonmin, lonmax = phi_lon[0], phi_lon[-1]

# User can specify what to normalize the radius by
# By default, normalize by the solar radius
if kw.rnorm is None:
    kw.rnorm = rsun
rr /= rnorm                                           
print (buffline)
print ("normalizing grid spacing by rnorm =", kw.rnorm)

# calculate grid spacing

# radial
dr = np.zeros(nr)
dr[:-1] = rr[:-1] - rr[1:]
dr[-1] = dr[-2]

# theta
dt = np.zeros(nt)
dt[:-1] = tt[:-1] - tt[1:]
dt[-1] = dt[-2]

# Generate the plot (1 row of 3 subplots)
fig, axs = plt.subplots(1, 3, sharey=True, figsize=(12,4))

# Plot delta_r
size = 0.5
axs[0].scatter(rr, dr, color='k', s=size)
# Make the "rval" location red and big
axs[0].scatter(r0, dr[ir0], color='r', s=3*size)

# Plot r0 * delta_theta
tt_lat = 180./np.pi*(np.pi/2. - tt)
axs[1].scatter(tt_lat, r0*dt, color='k', s=size)

# Plot r0 * sin(theta) * dphi
dphi = 2*np.pi/(2*nt)
axs[2].scatter(tt_lat, r0*sint*dphi, color='k', s=size)

# Set x-limits for axs[0]
plt.sca(axs[0])
ri, ro = np.min(rr), np.max(rr)
Dr = ro - ri
plt.xlim(ri - 0.01*Dr, ro + 0.01*Dr)
plt.xlabel(r'$r/({\rm{rnorm}} = %.1e\ \rm{cm})$' %rnorm, fontsize=12, **csfont)
plt.ylabel(r'$\delta r/\rm{rnorm}$')
# Get ticks everywhere
plt.minorticks_on()
plt.tick_params(top=True, right=True, direction='in', which='both')

# ... for axs[1]
plt.sca(axs[1])
plt.xlim(-90., 90.)
plt.xlabel('latitude', fontsize=12, **csfont)
plt.ylabel(r'$r_0\delta \theta/\rm{rnorm}$')
# Get ticks everywhere
plt.minorticks_on()
plt.tick_params(top=True, right=True, direction='in', which='both')

# ... for axs[2]
plt.sca(axs[2])
plt.xlim(-90., 90.)
plt.xlabel('latitude', fontsize=12, **csfont)
plt.ylabel(r'$r_0\sin\theta\delta \phi/\rm{rnorm}$')
# Get ticks everywhere
plt.minorticks_on()
plt.tick_params(top=True, right=True, direction='in', which='both')

# Make title(s)
axs[0].set_title(dirname_stripped + '\n' + 'Grid Resolution ', **csfont)
axs[1].set_title((r'$r_0/{\rm{rnorm}}=%.3f$' %r0) + '\n' +\
        (r'$\delta r_0/{\rm{rnorm}}=%1.3e$' %dr[ir0]))
axs[2].set_title((r'$r_0{\rm{max}}(\delta\theta)/{\rm{rnorm}}=%1.3e$' %(r0*np.max(dt))) + '\n' +\
        (r'$r_0{\rm{max}}(\sin\theta\delta\phi)/{\rm{rnorm}}=%1.3e$' %(r0*np.max(sint)*dphi)))

# Display the plot
plt.show()
