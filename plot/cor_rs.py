# I added something on LCD
import numpy as np
import matplotlib.pyplot as plt
from diagnostic_reading import *
import sys

# Get directory of run on which to do an,,sis
# and associated sub-directories
dirname = sys.argv[1]
datadir = dirname + '/data/'
plotdir = dirname + '/plots/'

# Read in command-line arguments (clas)
#clas = sys.argv[2:]
#nclas = len(clas)
#for ii in range(nclas):
#    cla = clas[ii]
#    if (cla == '-save'):
#        saveplot = True

# Create 'plotdir' if it doesn't exist already

if (not os.path.isdir(plotdir)):
        os.makedirs(plotdir)

# Get basic geometry info
(rr,tt,cost,sint,rr_depth,ri,ro,d) = np.load(datadir + 'grid_info.npy')
nr = len(rr); nt = len(tt)


# Read in area weights for each cell in the grid
(area_weights,areas,total_area,tt_weights,rr_weights) =\
            np.load(datadir + 'area_weights.npy')

# Read in raw Reynolds Stress correlations
(vr2_p, vt2_p, vp2_p, vrvp_p, vrvt_p, vtvp_p, vr2_m, vt2_m, vp2_m, 
        vrvp_m, vrvt_m, vtvp_m, fplus, fminus) =\
                np.load(datadir + 'rs_pm.npy')

# Add to get total RS components (filling factors already included)
vr2_t = vr2_p + vr2_m
vt2_t = vt2_p + vt2_m
vp2_t = vp2_p + vp2_m

vrvt_t = vrvt_p + vrvt_m
vrvp_t = vrvp_p + vrvp_m
vtvp_t = vtvp_p + vtvp_m

# Normalize so that the stresses are in m^2/s^2
t4 = 1.e4 # 'ten to the 4'

# diagonal parts
vr2_p /= t4; vr2_m /= t4; vr2_t /= t4
vt2_p /= t4; vt2_m /= t4; vt2_t /= t4
vp2_p /= t4; vp2_m /= t4; vp2_t /= t4

# off-diagonal parts
vrvt_p /= t4; vrvt_m /= t4; vrvt_t /= t4
vrvp_p /= t4; vrvp_m /= t4; vrvp_t /= t4
vtvp_p /= t4; vtvp_m /= t4; vtvp_t /= t4

# Plot the correlations in 2d on a square grid with depth on the x-axis
# and latitude on the y-axis

latvals = (np.pi/2. - tt)*180./np.pi
rr_2d, lat_2d = np.meshgrid(rr_depth, latvals)

# Use 'nn' for "non-normalized"
fignames = ['nn_vrvt', 'nn_vrvt_m', 'nn_vrvt_p',
        'nn_vrvp', 'nn_vrvp_m', 'nn_vrvp_p',
        'nn_vtvp', 'nn_vtvp_m', 'nn_vtvp_p']

titles = ["<vr'vt'> total", "<vr'vt'> downflow", 
                "<vr'vt'> upflow",
          "<vr'vp'> total", "<vr'vp'> downflow", 
                "<vr'vp'> upflow",      
          "<vt'vp'> total", "<vt'vp'> downflow", 
                "<vt'vp'> upflow"]

# max/min values
min_vrvt, max_vrvt = min(np.min(vrvt_t), np.min(vrvt_m), np.min(vrvt_p)),\
        max(np.max(vrvt_t), np.max(vrvt_m), np.max(vrvt_p))

min_vrvp, max_vrvp = min(np.min(vrvp_t), np.min(vrvp_m), np.min(vrvp_p)),\
        max(np.max(vrvp_t), np.max(vrvp_m), np.max(vrvp_p))

min_vtvp, max_vtvp = min(np.min(vtvp_t), np.min(vtvp_m), np.min(vtvp_p)),\
        max(np.max(vtvp_t), np.max(vtvp_m), np.max(vtvp_p))

minmax = [(min_vrvt, max_vrvt),(min_vrvt, max_vrvt),(min_vrvt, max_vrvt),\
        (min_vrvp, max_vrvp), (min_vrvp, max_vrvp), (min_vrvp, max_vrvp),\
        (min_vtvp, max_vtvp), (min_vtvp, max_vtvp), (min_vtvp, max_vtvp)]

count = 0
for arr in [vrvt_t, vrvt_m, vrvt_p,\
        vrvp_t, vrvp_m, vrvp_p,\
        vtvp_t, vtvp_m, vtvp_p]:
    
    (mymin, mymax) = minmax[count]
    plt.pcolormesh(rr_2d, lat_2d, arr, vmin=mymin, vmax=mymax)
    plt.colorbar()
    plt.xlabel('(ro - r)/(ro - ri)')
    plt.ylabel('latitude (degrees)')
    plt.title(titles[count] + ' ((m/s)^2)')
    savename = plotdir + fignames[count] + '.png'
    plt.savefig(savename, dpi = 300)
    plt.close()
    count += 1

# Plot the stress as a function of radius (plot rr_depth = 0.05 as
# the 'bottom of the TBL' for reference

# First average correlations over latitude
ttw_2d = tt_weights.reshape((nt,1))

# For the components involving vt', change sign when averaging over the 
# Southern hemisphere to avoid meaningless cancellation
ttw_2d_as = np.copy(ttw_2d) # 'as' for antisymmetric
ttw_2d_as[:nt/2] = -ttw_2d[:nt/2]

vrvt_tr = np.sum(vrvt_t*ttw_2d_as,axis=0)
vrvp_tr = np.sum(vrvp_t*ttw_2d,axis=0)
vtvp_tr = np.sum(vtvp_t*ttw_2d_as,axis=0)

vrvt_pr = np.sum(vrvt_p*ttw_2d_as,axis=0)
vrvp_pr = np.sum(vrvp_p*ttw_2d,axis=0)
vtvp_pr = np.sum(vtvp_p*ttw_2d_as,axis=0)

vrvt_mr = np.sum(vrvt_m*ttw_2d_as,axis=0)
vrvp_mr = np.sum(vrvp_m*ttw_2d,axis=0)
vtvp_mr = np.sum(vtvp_m*ttw_2d_as,axis=0)

# cor (vr', vt')
plt.plot(rr_depth[1:nr-1], vrvt_tr[1:nr-1], label = 'total')
plt.plot(rr_depth[1:nr-1], vrvt_mr[1:nr-1], label = 'downflow')
plt.plot(rr_depth[1:nr-1], vrvt_pr[1:nr-1], label = 'upflow')
# Mark depth = 5% line
ax = plt.gca()
mymin, mymax = ax.get_ylim()
plt.plot(0.05*np.ones(100), np.linspace(mymin, mymax, 100), 'k--', label='depth = 5%')
plt.legend()
plt.ylim((mymin, mymax))
plt.xlim((0,1))
plt.xlabel('(ro - r)/(ro - ri)')
plt.ylabel("<vr'vt'> ((m/s)^2)")
plt.title("<vr'vt'> vs. depth")
plt.savefig(plotdir + 'nn_vrvt_r.png', dpi=300)
plt.close()

# cor(vr',vp')
plt.plot(rr_depth[1:nr-1], vrvp_tr[1:nr-1], label = 'total')
plt.plot(rr_depth[1:nr-1], vrvp_mr[1:nr-1], label = 'downflow')
plt.plot(rr_depth[1:nr-1], vrvp_pr[1:nr-1], label = 'upflow')
# Mark depth = 5% line
ax = plt.gca()
mymin, mymax = ax.get_ylim()
plt.plot(0.05*np.ones(100), np.linspace(mymin, mymax,100), 'k--', label='depth = 5%')
plt.legend()
plt.xlim((0,1))
plt.ylim((mymin, mymax))
plt.xlabel('(ro - r)/(ro - ri)')
plt.ylabel("<vr'vp'> ((m/s)^2)")
plt.title("<vr'vp'> vs. depth")
plt.savefig(plotdir + 'nn_vrvp_r.png', dpi=300)
plt.close()

# cor(vt',vp')
plt.plot(rr_depth[1:nr-1], vtvp_tr[1:nr-1], label = 'total')
plt.plot(rr_depth[1:nr-1], vtvp_mr[1:nr-1], label = 'downflow')
plt.plot(rr_depth[1:nr-1], vtvp_pr[1:nr-1], label = 'upflow')
# Mark depth = 5% line
ax = plt.gca()
mymin, mymax = ax.get_ylim()
plt.plot(0.05*np.ones(100), np.linspace(mymin,mymax,100), 'k--', label='depth = 5%')
plt.legend()
plt.xlim((0,1))
plt.ylim((mymin, mymax))
plt.xlabel('(ro - r)/(ro - ri)')
plt.ylabel("<vt'vp'> ((m/s)^2)")
plt.title("<vt'vp'> vs. depth")
plt.savefig(plotdir + 'nn_vtvp_r.png', dpi=300)
plt.close()


#######################################################
# Plot correlation as a function of latitude 

# First average correlations over radius
rrw_2d = rr_weights.reshape((1,nr))

vrvt_tt = np.sum(vrvt_t*rrw_2d,axis=1)
vrvp_tt = np.sum(vrvp_t*rrw_2d,axis=1)
vtvp_tt = np.sum(vtvp_t*rrw_2d,axis=1)

vrvt_pt = np.sum(vrvt_p*rrw_2d,axis=1)
vrvp_pt = np.sum(vrvp_p*rrw_2d,axis=1)
vtvp_pt = np.sum(vtvp_p*rrw_2d,axis=1)

vrvt_mt = np.sum(vrvt_m*rrw_2d,axis=1)
vrvp_mt = np.sum(vrvp_m*rrw_2d,axis=1)
vtvp_mt = np.sum(vtvp_m*rrw_2d,axis=1)


# cor (vr', vt')
plt.plot(latvals, vrvt_tt, label = 'total')
plt.plot(latvals, vrvt_mt, label = 'downflow')
plt.plot(latvals, vrvt_pt, label = 'upflow')
plt.legend()
plt.xlim((-90,90))
plt.xlabel('latitude (degrees)')
plt.ylabel("<vr'vt'> ((m/s)^2)")
plt.title('correlation vs. latitude')
plt.savefig(plotdir + 'nn_vrvt_t.png', dpi=300)
plt.close()

# cor (vr', vp')
plt.plot(latvals, vrvp_tt, label = 'total')
plt.plot(latvals, vrvp_mt, label = 'downflow')
plt.plot(latvals, vrvp_pt, label = 'upflow')
plt.legend()
plt.xlim((-90,90))
plt.xlabel('latitude (degrees)')
plt.ylabel("<vr'vp'> ((m/s)^2)")
plt.title('correlation vs. latitude')
plt.savefig(plotdir + 'nn_vrvp_t.png', dpi=300)
plt.close()

# cor (vt', vp')
plt.plot(latvals, vtvp_tt, label = 'total')
plt.plot(latvals, vtvp_mt, label = 'downflow')
plt.plot(latvals, vtvp_pt, label = 'upflow')
plt.legend()
plt.xlim((-90,90))
plt.xlabel('latitude (degrees)')
plt.ylabel("<vt'vp'> ((m/s)^2)")
plt.title('correlation vs. latitude')
plt.savefig(plotdir + 'nn_vtvp_t.png', dpi=300)
plt.close()

