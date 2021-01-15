import numpy as np
import matplotlib.pyplot as plt
from diagnostic_reading import ShellSlice, ReferenceState
import sys, os
from common import *
from get_parameter import get_parameter

dirname = sys.argv[1]
dirname_stripped = strip_dirname(dirname)
datadir = dirname + '/data/'

# Create 'datadir' if it doesn't exist already
if (not os.path.isdir(datadir)):
    os.makedirs(datadir)

om = get_parameter(dirname, 'angular_velocity')
cp = get_parameter(dirname, 'pressure_specific_heat')

slicedir = dirname + '/Shell_Slices/'
files = os.listdir(slicedir)
files.sort()

rr,tt,cost,sint,rr_depth,ri,ro,d = np.load(datadir + 'grid_info.npy')
nr = len(rr)
nt = len(tt)
nph = 2*nt
sint_2d = sint.reshape((1, nt))

ref = ReferenceState(dirname + '/reference', '')
g = ref.gravity/100

print ('Reading ' + slicedir + files[-1] + ' ...')
a = ShellSlice(slicedir + files[-1], '')
g_rvals = g[a.inds]

vr = a.vals[:, :, :, a.lut[1], :]/100
vt = a.vals[:, :, :, a.lut[2], :]/100
vp = a.vals[:, :, :, a.lut[3], :]/100
s = a.vals[:, :, :, a.lut[64], :]

smean_file = get_widest_range_file(datadir, 's_spherical_mean')
smean = np.load(datadir + smean_file)
s_fluc = s - smean[a.inds].reshape((1, 1, a.nr, 1))
v_drift = -s_fluc*g_rvals/2/om/cp

var_dict = {'vr': vr, 'vt': vt, 'vp': vp, 's': s_fluc, 'v_drift': v_drift}


for i in range(a.nr):
    vr_loc = vr[:, :, i, :]
    where_p = np.where(vr_loc > 0)
    where_m = np.where(vr_loc < 0)
    for varname in ['vr', 'vt', 'vp', 's', 'v_drift']:
        plotdir = dirname + '/plots/sslice_dist/' + varname + '/'
        if (not os.path.isdir(plotdir)):
            os.makedirs(plotdir)
        var = var_dict[varname]
        var_loc = var[:, :, i, :]
        var_p = var_loc[where_p]
        var_m = var_loc[where_m]
        
        histvals, binvals, patches = plt.hist(var_p, bins=100, density=True)
        meanval = np.sum(binvals[1:]*histvals)/np.sum(histvals)
        ymin, ymax = plt.gca().get_ylim()
        yvals = np.linspace(ymin, ymax, 100)
        plt.plot(meanval*np.ones(100), yvals, label='mean')
        plt.legend()
        plt.title(dirname_stripped + ', ' + varname + ', depth = %f, upflow' %rr_depth[a.inds[i]])
        plt.minorticks_on()
        plt.tick_params(top='on', right='on', direction='in', which='both')
        savefile = plotdir + dirname_stripped + '_' + varname + '_depth' +\
            str(i).zfill(3) + '_upflow.png'
        print ('saving ' + savefile + ' ...')
        plt.savefig(savefile, dpi=300)
        plt.close()
        
        histvals, binvals, patches = plt.hist(var_m, bins=100, density=True)
        meanval = np.sum(binvals[1:]*histvals)/np.sum(histvals)
        ymin, ymax = plt.gca().get_ylim()
        yvals = np.linspace(ymin, ymax, 100)
        plt.plot(meanval*np.ones(100), yvals, label='mean')
        plt.legend()

        plt.title(dirname_stripped + ', ' + varname + ', depth = %f, downflow' %rr_depth[a.inds[i]])
        plt.minorticks_on()
        plt.tick_params(top='on', right='on', direction='in', which='both')

        savefile = plotdir + dirname_stripped + '_' + varname + '_depth' +\
            str(i).zfill(3) + '_downflow.png'
        print ('saving ' + savefile + ' ...')
        plt.savefig(savefile, dpi=300)
        plt.close()
        
        histvals, binvals, patches = plt.hist(var_loc.flatten(), bins=100, density=True)
        meanval = np.sum(binvals[1:]*histvals)/np.sum(histvals)
        ymin, ymax = plt.gca().get_ylim()
        yvals = np.linspace(ymin, ymax, 100)
        plt.plot(meanval*np.ones(100), yvals, label='mean')
        plt.legend()
        plt.title(dirname_stripped + ', ' + varname + ', depth = %f, full' %rr_depth[a.inds[i]])
        plt.minorticks_on()
        plt.tick_params(top='on', right='on', direction='in', which='both')
        savefile = plotdir + dirname_stripped + '_' + varname + '_depth' +\
            str(i).zfill(3) + '_full.png'
        print ('saving ' + savefile + ' ...')
        plt.savefig(savefile, dpi=300)
        plt.close()
