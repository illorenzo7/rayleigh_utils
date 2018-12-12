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

# Read in correlation data for vr' and s'
cor = np.load(datadir + 'cor_vr_s.npy')

cor_df = np.load(datadir + 'cor_vr_s_df.npy')
cor_uf = np.load(datadir + 'cor_vr_s_uf.npy')

cor_r = np.load(datadir + 'cor_vr_s_r.npy') 
cor_t = np.load(datadir + 'cor_vr_s_t.npy') 

cor_df_r = np.load(datadir + 'cor_vr_s_df_r.npy')
cor_df_t = np.load(datadir + 'cor_vr_s_df_t.npy')

cor_uf_r = np.load(datadir + 'cor_vr_s_uf_r.npy')
cor_uf_t = np.load(datadir + 'cor_vr_s_uf_t.npy')

# Plot the correlations in 2d on a square grid with depth on the x-axis
# and latitude on the y-axis

latvals = (np.pi/2. - tt)*180./np.pi
rr_2d, lat_2d = np.meshgrid(rr_depth, latvals)

cornames = ['cor_vr_s', 'cor_vr_s_df', 'cor_vr_s_uf']
titles = ["cor(vr', s') total", "cor(vr', s') downflow", "cor(vr', s') upflow"]

count = 0
for arr in [cor, cor_df, cor_uf]:
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
# Exclude problematic endponts

plt.plot(rr_depth[1:nr-1], cor_r[1:nr-1], label = 'total')
plt.plot(rr_depth[1:nr-1], cor_df_r[1:nr-1], label = 'downflow')
plt.plot(rr_depth[1:nr-1], cor_uf_r[1:nr-1], label = 'upflow')
# Mark depth = 5% line
plt.plot(0.05*np.ones(100), np.linspace(0,1,100), 'k--', label='depth = 5%')
plt.legend()
plt.xlim((0,1))
plt.ylim((0,1))
plt.xlabel('(ro - r)/(ro - ri)')
plt.ylabel("cor(vr', s')")
plt.title('correlation vs. depth')
plt.savefig(plotdir + 'cor_vr_s_r.png', dpi=300)
plt.close()

plt.plot(latvals, cor_t, label = 'total')
plt.plot(latvals, cor_df_t, label = 'downflow')
plt.plot(latvals, cor_uf_t, label = 'upflow')
plt.legend()
plt.xlim((-90,90))
plt.ylim((0,1))
plt.xlabel('latitude (degrees)')
plt.ylabel("cor(vr', s')")
plt.title('correlation vs. latitude')
plt.savefig(plotdir + 'cor_vr_s_t.png', dpi=300)
plt.close()

