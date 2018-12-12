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
                np.load(datadir + 'rs_pm_raw.npy')

# Multiply by filling factors to get total RS components
vr2_t = fplus*vr2_p + fminus*vr2_m
vt2_t = fplus*vt2_p + fminus*vt2_m
vp2_t = fplus*vp2_p + fminus*vp2_m

vrvt_t = fplus*vrvt_p + fminus*vrvt_m
vrvp_t = fplus*vrvp_p + fminus*vrvp_m
vtvp_t = fplus*vtvp_p + fminus*vtvp_m

# Compute correlation coefficients
cor_vrvt_t = vrvt_t/np.sqrt(vr2_t*vt2_t)
cor_vrvp_t = vrvp_t/np.sqrt(vr2_t*vp2_t)
cor_vtvp_t = vtvp_t/np.sqrt(vt2_t*vp2_t)

cor_vrvt_p = vrvt_p/np.sqrt(vr2_p*vt2_p)
cor_vrvp_p = vrvp_p/np.sqrt(vr2_p*vp2_p)
cor_vtvp_p = vtvp_p/np.sqrt(vt2_p*vp2_p)

cor_vrvt_m = vrvt_m/np.sqrt(vr2_m*vt2_m)
cor_vrvp_m = vrvp_m/np.sqrt(vr2_m*vp2_m)
cor_vtvp_m = vtvp_m/np.sqrt(vt2_m*vp2_m)


# Plot the correlations in 2d on a square grid with depth on the x-axis
# and latitude on the y-axis

latvals = (np.pi/2. - tt)*180./np.pi
rr_2d, lat_2d = np.meshgrid(rr_depth, latvals)

cornames = ['cor_vrvt', 'cor_vrvt_m', 'cor_vrvt_p',
        'cor_vrvp', 'cor_vrvp_m', 'cor_vrvp_p',
        'cor_vtvp', 'cor_vtvp_m', 'cor_vtvp_p']

titles = ["cor(vr', vt') total", "cor(vr', vt') downflow", 
                "cor(vr', vt') upflow",
          "cor(vr', vp') total", "cor(vr', vp') downflow", 
                "cor(vr', vp') upflow",      
          "cor(vt', vp') total", "cor(vt', vp') downflow", 
                "cor(vt', vp') upflow"]

count = 0
for arr in [cor_vrvt_t, cor_vrvt_m, cor_vrvt_p,\
        cor_vrvp_t, cor_vrvp_m, cor_vrvp_p,\
        cor_vtvp_t, cor_vtvp_m, cor_vtvp_p]:

    plt.pcolormesh(rr_2d, lat_2d, arr, vmin = -1., vmax = 1.)
    plt.colorbar()
    plt.xlabel('(ro - r)/(ro - ri)')
    plt.ylabel('latitude (degrees)')
    plt.title(titles[count])
    savename = plotdir + cornames[count] + '.png'
    plt.savefig(savename, dpi = 300)
    plt.close()
    count += 1

# Plot correlation as a function of radius (plot rr_depth = 0.05 as
# the 'bottom of the TBL' for reference

# First average correlations over latitude
ttw_2d = tt_weights.reshape((nt,1))

# For the components involving vt', change sign when averaging over the 
# Southern hemisphere to avoid meaningless cancellation
ttw_2d_as = np.copy(ttw_2d) # 'as' for antisymmetric
ttw_2d_as[:nt/2] = -ttw_2d[:nt/2]

cor_vrvt_tr = np.sum(cor_vrvt_t*ttw_2d_as,axis=0)
cor_vrvp_tr = np.sum(cor_vrvp_t*ttw_2d,axis=0)
cor_vtvp_tr = np.sum(cor_vtvp_t*ttw_2d_as,axis=0)

cor_vrvt_pr = np.sum(cor_vrvt_p*ttw_2d_as,axis=0)
cor_vrvp_pr = np.sum(cor_vrvp_p*ttw_2d,axis=0)
cor_vtvp_pr = np.sum(cor_vtvp_p*ttw_2d_as,axis=0)

cor_vrvt_mr = np.sum(cor_vrvt_m*ttw_2d_as,axis=0)
cor_vrvp_mr = np.sum(cor_vrvp_m*ttw_2d,axis=0)
cor_vtvp_mr = np.sum(cor_vtvp_m*ttw_2d_as,axis=0)

# cor (vr', vt')
plt.plot(rr_depth[1:nr-1], cor_vrvt_tr[1:nr-1], label = 'total')
plt.plot(rr_depth[1:nr-1], cor_vrvt_mr[1:nr-1], label = 'downflow')
plt.plot(rr_depth[1:nr-1], cor_vrvt_pr[1:nr-1], label = 'upflow')
# Mark depth = 5% line
plt.plot(0.05*np.ones(100), np.linspace(-1,1,100), 'k--', label='depth = 5%')
plt.legend()
plt.xlim((0,1))
plt.ylim((-1,1))
plt.xlabel('(ro - r)/(ro - ri)')
plt.ylabel("cor(vr', vt')")
plt.title("cor(vr', vt') vs. depth")
plt.savefig(plotdir + 'cor_vrvt_r.png', dpi=300)
plt.close()

# cor(vr',vp')
plt.plot(rr_depth[1:nr-1], cor_vrvp_tr[1:nr-1], label = 'total')
plt.plot(rr_depth[1:nr-1], cor_vrvp_mr[1:nr-1], label = 'downflow')
plt.plot(rr_depth[1:nr-1], cor_vrvp_pr[1:nr-1], label = 'upflow')
# Mark depth = 5% line
plt.plot(0.05*np.ones(100), np.linspace(-1,1,100), 'k--', label='depth = 5%')
plt.legend()
plt.xlim((0,1))
plt.ylim((-1,1))
plt.xlabel('(ro - r)/(ro - ri)')
plt.ylabel("cor(vr', vp')")
plt.title("cor(vr', vp') vs. depth")
plt.savefig(plotdir + 'cor_vrvp_r.png', dpi=300)
plt.close()

# cor(vt',vp')
plt.plot(rr_depth[1:nr-1], cor_vtvp_tr[1:nr-1], label = 'total')
plt.plot(rr_depth[1:nr-1], cor_vtvp_mr[1:nr-1], label = 'downflow')
plt.plot(rr_depth[1:nr-1], cor_vtvp_pr[1:nr-1], label = 'upflow')
# Mark depth = 5% line
plt.plot(0.05*np.ones(100), np.linspace(-1,1,100), 'k--', label='depth = 5%')
plt.legend()
plt.xlim((0,1))
plt.ylim((-1,1))
plt.xlabel('(ro - r)/(ro - ri)')
plt.ylabel("cor(vt', vp')")
plt.title("cor(vt', vp') vs. depth")
plt.savefig(plotdir + 'cor_vtvp_r.png', dpi=300)
plt.close()


#######################################################
# Plot correlation as a function of latitude 

# First average correlations over radius
rrw_2d = rr_weights.reshape((1,nr))

cor_vrvt_tt = np.sum(cor_vrvt_t*rrw_2d,axis=1)
cor_vrvp_tt = np.sum(cor_vrvp_t*rrw_2d,axis=1)
cor_vtvp_tt = np.sum(cor_vtvp_t*rrw_2d,axis=1)

cor_vrvt_pt = np.sum(cor_vrvt_p*rrw_2d,axis=1)
cor_vrvp_pt = np.sum(cor_vrvp_p*rrw_2d,axis=1)
cor_vtvp_pt = np.sum(cor_vtvp_p*rrw_2d,axis=1)

cor_vrvt_mt = np.sum(cor_vrvt_m*rrw_2d,axis=1)
cor_vrvp_mt = np.sum(cor_vrvp_m*rrw_2d,axis=1)
cor_vtvp_mt = np.sum(cor_vtvp_m*rrw_2d,axis=1)


# cor (vr', vt')
plt.plot(latvals, cor_vrvt_tt, label = 'total')
plt.plot(latvals, cor_vrvt_mt, label = 'downflow')
plt.plot(latvals, cor_vrvt_pt, label = 'upflow')
plt.legend()
plt.xlim((-90,90))
plt.ylim((-1,1))
plt.xlabel('latitude (degrees)')
plt.ylabel("cor(vr', vt')")
plt.title('correlation vs. latitude')
plt.savefig(plotdir + 'cor_vrvt_t.png', dpi=300)
plt.close()

# cor (vr', vp')
plt.plot(latvals, cor_vrvp_tt, label = 'total')
plt.plot(latvals, cor_vrvp_mt, label = 'downflow')
plt.plot(latvals, cor_vrvp_pt, label = 'upflow')
plt.legend()
plt.xlim((-90,90))
plt.ylim((-1,1))
plt.xlabel('latitude (degrees)')
plt.ylabel("cor(vr', vp')")
plt.title('correlation vs. latitude')
plt.savefig(plotdir + 'cor_vrvp_t.png', dpi=300)
plt.close()

# cor (vt', vp')
plt.plot(latvals, cor_vtvp_tt, label = 'total')
plt.plot(latvals, cor_vtvp_mt, label = 'downflow')
plt.plot(latvals, cor_vtvp_pt, label = 'upflow')
plt.legend()
plt.xlim((-90,90))
plt.ylim((-1,1))
plt.xlabel('latitude (degrees)')
plt.ylabel("cor(vt', vp')")
plt.title('correlation vs. latitude')
plt.savefig(plotdir + 'cor_vtvp_t.png', dpi=300)
plt.close()

