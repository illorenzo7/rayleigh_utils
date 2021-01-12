import sys, os
import numpy as np
from common import strip_dirname, get_widest_range_file, get_iters_from_file,\
    interpx, interpy

# Get the name of the run directory
dirname = sys.argv[1]
# Get the stripped name to use in file naming
dirname_stripped = strip_dirname(dirname)

# Find the relevant place to store the data, and create the directory if it
# doesn't already exist
datadir = dirname + '/data/'
plotdir = dirname + '/plots/vr_s_dist_df/'
if (not os.path.isdir(datadir)):
    os.makedirs(datadir)

# Get grid info (if it's not computed already using grid_info.py, this will fail)
rr,tt,cost,sint,rr_depth,ri,ro,d = np.load(datadir + 'grid_info.npy')
rr /= 100; ro /= 100; ri /= 100; d /= 100 # trajectories are in meter-land
nr = len(rr)
nt = len(tt)

# Get location of thermal BL 
ir_tbl = np.load(datadir + dirname_stripped + '_ir_tbl.npy')
rr_tbl = rr[ir_tbl]

# Get distribution of vr vs. vphi
vr_vp_dist_file = get_widest_range_file(datadir, 'vr_vp_dist_df_1_lowlat_tbl')
vr_vp_dist_tbl = np.load(datadir + vr_vp_dist_file)
vr_vp_dist_tbl_totalcount = np.sum(vr_vp_dist_tbl)

vr_bincenters, vp_bincenters = np.load(datadir + 'vrvp_bincenters.npy')
vrvals_orig, vpvals_orig = vr_bincenters[ir_tbl], vp_bincenters[ir_tbl]
        
# Load the trajectories
trajectories = np.load(datadir + dirname_stripped + '_eq_trajectories.npy')
ntraj = len(trajectories)
npoints = (len(trajectories[0,:]) - 9)//5

probs = trajectories[:, -6]
pen_depths = trajectories[:, -1]
where_in_domain = np.where(rr_tbl - pen_depths > ri)

min_rr = rr_tbl - np.max(pen_depths)

ir1 = np.argmin(np.abs(rr - rr_tbl))
if rr[ir1] > rr_tbl: # Want rr[ir1] to be just less than rr_tbl
    ir1 += 1

ir2 = np.argmin(np.abs(rr - min_rr))
if rr[ir2] < min_rr:
    ir2 -= 1

nr_synth = ir2 - ir1 + 1

vrvals_traj = trajectories[where_in_domain, npoints:2*npoints]
vpvals_traj = trajectories[where_in_domain, 2*npoints:3*npoints]

vrmin_traj = np.min(vrvals_traj)
vrmax_traj = 0
vpmin_traj = np.min(vpvals_traj)
vpmax_traj = np.max(vpvals_traj)

# Get bininfo for vr/vp (df)
nbins_vr = 100
nbins_vp = 100

vr_space = (0 - vrmin_traj)/nbins_vr
vp_space = (vpmax_traj - vpmin_traj)/nbins_vp
bin_area = vr_space*vp_space

vr_bincenters = np.linspace(vrmin_traj + 0.5*vr_space,\
                            vrmax_traj - 0.5*vr_space, nbins_vr)
vp_bincenters = np.linspace(vpmin_traj + 0.5*vp_space,\
                            vpmax_traj - 0.5*vp_space, nbins_vp)

# Array to hold the distribution
vr_vp_dist_lowlat_synthetic = np.zeros((nr_synth, nbins_vr, nbins_vp))

count = 0
for ir in range(ir1, ir2):
#for ir in range(ir1, ir1+1): # for debugging purposes
    print ('Calculating synthetic vr/vp distribution for depth %03i' %ir + ' ...')
    
    r_want = rr[ir]
    
    for itraj in range(ntraj):
    #for itraj in range(1):
        if (rr_tbl - pen_depths[itraj] < r_want):
            prob_traj = probs[itraj]
            t_traj = trajectories[itraj, :npoints] 
            vr_traj= trajectories[itraj, npoints:2*npoints]
            vp_traj = trajectories[itraj, 2*npoints:3*npoints]
            r_traj = trajectories[itraj, 3*npoints:4*npoints]
            
            it_close1 = np.argmin(np.abs(r_traj - r_want))
            t_close1 = t_traj[it_close1]
            r_close1 = r_traj[it_close1]
            if (r_traj[it_close1] > r_want):
                it_close2 = it_close1 + 1
                if (it_close2 >= npoints):
                    it_close2 = it_close1 - 1
            else:
                it_close2 = it_close1 - 1
            t_close2 = t_traj[it_close2]
            r_close2 = r_traj[it_close2]
            t_want = interpx(t_close1, r_close1, t_close2, r_close2, r_want)
            
            vr_loc = interpy(t_close1, vr_traj[it_close1],\
                             t_close2, vr_traj[it_close2], t_want)
     
            vp_loc = interpy(t_close1, vp_traj[it_close1],\
                             t_close2, vp_traj[it_close2], t_want)       
    
            # bin vr, putting all values outside [minvr,maxvr] 
            # in the
            # closest bin in range (the outermost bins)
            if (vr_loc < vrmin_traj + vr_space):
                bin_index_vr = 0
            elif (vr_loc > vrmax_traj - vr_space):
                bin_index_vr = nbins_vr - 1
            else:
                bin_index_vr = np.floor((vr_loc - vrmin_traj)\
                        /vr_space)
                bin_index_vr = int(bin_index_vr)
            
            # bin vp according to the same process
            if (vp_loc < vpmin_traj + vp_space):
                bin_index_vp = 0
            elif (vp_loc > vpmax_traj - vp_space):
                bin_index_vp = nbins_vp - 1
            else:
                bin_index_vp = np.floor((vp_loc - vpmin_traj)\
                        /vp_space)
                bin_index_vp = int(bin_index_vp)
    
            vr_vp_dist_lowlat_synthetic[count, bin_index_vr,bin_index_vp] +=\
                prob_traj 
    count += 1

vr_vp_dist_lowlat_synthetic *= vr_vp_dist_tbl_totalcount

savefile = datadir + dirname_stripped + '_vr_vp_dist_1_lowlat_synthetic.npy'   
print ('Saving synthetic distribution at ' + savefile + ' ...') 
np.save(savefile, (vr_vp_dist_lowlat_synthetic, vr_bincenters, vp_bincenters,\
                   ir1, ir2, vr_vp_dist_tbl_totalcount))