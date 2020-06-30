import numpy as np
import time, pickle
import sys, os
from scipy.interpolate import RegularGridInterpolator
sys.path.append(os.environ['raco'])
sys.path.append(os.environ['rapp'])
from common import get_file_lists, strip_dirname
from translate_times import translate_times
from rayleigh_diagnostics import Spherical_3D
from get_parameter import get_parameter
from time_scales import compute_Prot, compute_tdt
from get_eq import get_eq

# Get command line arguments
dirname = sys.argv[1]
dirname_stripped = strip_dirname(dirname)

# Data with Shell_Slices
radatadir = dirname + "/Spherical_3D/"
savedir = dirname + "/Cartesian_3D/"
if not os.path.isdir(savedir):
    os.makedirs(savedir)

file_list, int_file_list, nfiles = get_file_lists(radatadir)

iiter = nfiles - 1 # by default plot the last iteration in Spherical_3D

args = sys.argv[2:]
nargs = len(args)
ncart = 400 # By default make 3D data cube with 400 points per side

# By default, interpolate all variables possible
interpgrid = True
interpv = True
interpB = True
interpthermo = True
interpvort = True
for i in range(nargs):
    arg = args[i]
    if arg == '-iter':
        desired_iter = int(args[i+1])
        iiter = np.argmin(np.abs(int_file_list - desired_iter))
    elif arg == '-sec':
        time = float(args[i+1])
        di_trans = translate_times(time, dirname, translate_from='sec')
        desired_iter = di_trans['val_iter']
        iiter = np.argmin(np.abs(int_file_list - desired_iter))
    elif arg == '-day':
        time = float(args[i+1])
        di_trans = translate_times(time, dirname, translate_from='day')
        desired_iter = di_trans['val_iter']
        iiter = np.argmin(np.abs(int_file_list - desired_iter))
    elif arg == '-prot':
        time = float(args[i+1])
        di_trans = translate_times(time, dirname, translate_from='prot')
        desired_iter = di_trans['val_iter']
        iiter = np.argmin(np.abs(int_file_list - desired_iter))
    elif arg == '-clon':
        clon = float(args[i+1])
    elif arg == '-ncart':
        ncart = int(args[i+1])
    elif arg == '-nogrid':
        intergrid = False
    elif arg == '-nov':
        interpv = False
    elif arg == '-noB':
        interpB = False
    elif arg == '-nothermo':
        interpthermo = False
    elif arg == '-novort':
        intervport = False

iter_val = int_file_list[iiter]
fname = file_list[iiter]

print ("Reading fluid variables for Spherical_3D/%08i" %iter_val)
# Always read vr_3d, since we use it for grid info
vr_3d = Spherical_3D(fname + '_0001', radatadir)

# Get spherical grid info
nr = vr_3d.nr
nt = vr_3d.ntheta
nphi = vr_3d.nphi
print('Spherical dimensions (nphi,nth,nr)=({:d},{:d},{:d})'.format(nphi,nt,nr,))
r = (vr_3d.r[::-1]/np.max(vr_3d.r)).reshape((1, 1, nr))
theta = (vr_3d.theta[::-1]).reshape((1, nt, 1))
phi = (np.linspace(0, 2*np.pi, nphi+1)).reshape((nphi+1, 1, 1))
zero = np.zeros((nphi, nt, nr))

n_x = n_y = n_z = ncart
# Set up the new Cartesian grid to interpolate onto
print("Interpolating onto Cartesian grid (nx,ny,nz)=({:d},{:d},{:d})".format(n_x,n_y,n_z,))
zero_cart = np.zeros((n_x, n_y, n_z))
x_u = zero_cart + np.linspace(-1, 1, n_x)[:, None, None]
y_u = zero_cart + np.linspace(-1, 1, n_y)[None, :, None]
z_u = zero_cart + np.linspace(-1, 1, n_z)[None, None, :]

r_u = (x_u**2 + y_u**2 + z_u**2)**0.5
theta_u = np.arccos(z_u/r_u)
phi_u = np.arctan2(y_u,x_u) + np.pi

# Write the grid info to dictionary and save it
if interpgrid:
    savefile = savedir + fname + '_grid.pkl'
    print ('Saving grid data at ' + savefile)
    f = open(savefile, 'wb')
    pickle.dump({'x': x_u, 'y': y_u, 'z': z_u, 'r': r_u, 'theta': theta_u,\
            'phi': phi_u}, f, protocol=4)
    f.close()

# Make flattened arrays of original coordinates and interpolation coordinates
original_coords_flat = (phi.flatten(), theta.flatten(), r.flatten())
new_coords_flat = (phi_u.flatten(), theta_u.flatten(), r_u.flatten())
points = np.array(new_coords_flat).T

start_overall = time.time()
# Start interpolating
if interpv:
    vt_3d = Spherical_3D(fname + '_0002', radatadir)
    vp_3d = Spherical_3D(fname + '_0003', radatadir)
    vr = vr_3d.vals
    vt = vt_3d.vals
    vp = vp_3d.vals

    # vr
    print ("Interpolating vr")
    start = time.time()
    data = np.copy(vr[:, ::-1, ::-1])
    data = np.append(data[:,:,:],np.expand_dims(data[0,:,:], axis=0), axis=0)
    F_interp = RegularGridInterpolator(original_coords_flat, data,\
            bounds_error=False, fill_value=0.0)
    vr_cart = F_interp(points).reshape((n_x, n_y, n_z))
    end = time.time()
    print ("Interpolating for vr took %.1f sec" %(end - start))

    # vt
    print ("Interpolating vt")
    start = time.time()
    data = np.copy(vt[:, ::-1, ::-1])
    data = np.append(data[:,:,:],np.expand_dims(data[0,:,:], axis=0), axis=0)
    F_interp = RegularGridInterpolator(original_coords_flat, data,\
            bounds_error=False, fill_value=0.0)
    vt_cart = F_interp(points).reshape((n_x, n_y, n_z))
    end = time.time()
    print ("Interpolating for vt took %.1f sec" %(end - start))

    # vp
    print ("Interpolating vp")
    start = time.time()
    data = np.copy(vp[:, ::-1, ::-1])
    data = np.append(data[:,:,:],np.expand_dims(data[0,:,:], axis=0), axis=0)
    F_interp = RegularGridInterpolator(original_coords_flat, data,\
            bounds_error=False, fill_value=0.0)
    vp_cart = F_interp(points).reshape((n_x, n_y, n_z))
    end = time.time()
    print ("Interpolating for vp took %.1f sec" %(end - start))

    # vr_prime
    print ("Interpolating vr_prime")
    start = time.time()
    data = np.copy(vr[:, ::-1, ::-1])
    data = data - np.mean(data, axis=0).reshape((1, nt, nr))
    data = np.append(data[:,:,:],np.expand_dims(data[0,:,:], axis=0), axis=0)
    F_interp = RegularGridInterpolator(original_coords_flat, data,\
            bounds_error=False, fill_value=0.0)
    vr_prime_cart = F_interp(points).reshape((n_x, n_y, n_z))
    end = time.time()
    print ("Interpolating for vr_prime took %.1f sec" %(end - start))

    # vt
    print ("Interpolating vt_prime")
    start = time.time()
    data = np.copy(vt[:, ::-1, ::-1])
    data = data - np.mean(data, axis=0).reshape((1, nt, nr))
    data = np.append(data[:,:,:],np.expand_dims(data[0,:,:], axis=0), axis=0)
    F_interp = RegularGridInterpolator(original_coords_flat, data,\
            bounds_error=False, fill_value=0.0)
    vt_prime_cart = F_interp(points).reshape((n_x, n_y, n_z))
    end = time.time()
    print ("Interpolating for vt_prime took %.1f sec" %(end - start))

    # vp_prime
    print ("Interpolating vp_prime")
    start = time.time()
    data = np.copy(vp[:, ::-1, ::-1])
    data = data - np.mean(data, axis=0).reshape((1, nt, nr))
    data = np.append(data[:,:,:],np.expand_dims(data[0,:,:], axis=0), axis=0)
    F_interp = RegularGridInterpolator(original_coords_flat, data,\
            bounds_error=False, fill_value=0.0)
    vp_prime_cart = F_interp(points).reshape((n_x, n_y, n_z))
    end = time.time()
    print ("Interpolating for vp_prime took %.1f sec" %(end - start))

    # Write the grid info to dictionary and save it
    savefile = savedir + fname + '_v.pkl'
    print ('Saving velocity data at ' + savefile)
    f = open(savefile, 'wb')
    pickle.dump({'vr': vr_cart, 'vt': vt_cart, 'vp': vp_cart,\
            'vr_prime': vr_prime_cart, 'vt_prime': vt_prime_cart,\
            'vp_prime': vp_prime_cart}, f, protocol=4)
    f.close()

    # Compute the Cartesian velocities
    print ("Computing Cartesian interpolated velocities, vx, vy, vz")
    vx = -(vr_cart*np.sin(theta_u) + vt_cart*np.cos(theta_u))*np.cos(phi_u) +\
            vp_cart*np.sin(phi_u)
    vy = -(vr_cart*np.sin(theta_u) + vt_cart*np.cos(theta_u))*np.sin(phi_u) -\
            vp_cart*np.cos(phi_u)
    vz = vr_cart*np.cos(theta_u)  - vt_cart*np.sin(theta_u)

    vx_prime = -(vr_prime_cart*np.sin(theta_u) +\
            vt_prime_cart*np.cos(theta_u))*np.cos(phi_u) +\
            vp_prime_cart*np.sin(phi_u)
    vy_prime = -(vr_prime_cart*np.sin(theta_u) +\
            vt_prime_cart*np.cos(theta_u))*np.sin(phi_u) -\
            vp_prime_cart*np.cos(phi_u)
    vz_prime = vr_prime_cart*np.cos(theta_u)  - vt_prime_cart*np.sin(theta_u)

    # Write the Cartesian velocity data to dictionary and save it
    savefile = savedir + fname + '_vcart.pkl'
    print ('Saving velocity data at ' + savefile)
    f = open(savefile, 'wb')
    pickle.dump({'vx': vx, 'vy': vy, 'vz': vz,\
        'vx_prime': vx_prime, 'vy_prime': vy_prime, 'vz_prime': vz_prime},\
        f, protocol=4)
    f.close()

if interpthermo:
    # First interpolate reference state
    eq = get_eq(dirname)
    rhobar = eq.density.reshape((1, 1, nr)) + zero
    pbar = eq.pressure.reshape((1, 1, nr)) + zero
    tbar = eq.temperature.reshape((1, 1, nr)) + zero

    # p
    print ("Interpolating rhobar")
    start = time.time()
    data = np.copy(rhobar[:, ::-1, ::-1])
    data = np.append(data[:,:,:],np.expand_dims(data[0,:,:], axis=0), axis=0)
    F_interp = RegularGridInterpolator(original_coords_flat, data,\
            bounds_error=False, fill_value=0.0)
    rhobar_cart = F_interp(points).reshape((n_x, n_y, n_z))
    end = time.time()
    print ("Interpolating for rhobar took %.1f sec" %(end - start))

    # pbar
    print ("Interpolating pbar")
    start = time.time()
    data = np.copy(pbar[:, ::-1, ::-1])
    data = np.append(data[:,:,:],np.expand_dims(data[0,:,:], axis=0), axis=0)
    F_interp = RegularGridInterpolator(original_coords_flat, data,\
            bounds_error=False, fill_value=0.0)
    pbar_cart = F_interp(points).reshape((n_x, n_y, n_z))
    end = time.time()
    print ("Interpolating for rhobar took %.1f sec" %(end - start))

    # tbar
    print ("Interpolating tbar")
    start = time.time()
    data = np.copy(tbar[:, ::-1, ::-1])
    data = np.append(data[:,:,:],np.expand_dims(data[0,:,:], axis=0), axis=0)
    F_interp = RegularGridInterpolator(original_coords_flat, data,\
            bounds_error=False, fill_value=0.0)
    tbar_cart = F_interp(points).reshape((n_x, n_y, n_z))
    end = time.time()
    print ("Interpolating for tbar took %.1f sec" %(end - start))

    # Now the actual thermal perturbations
    s_3d = Spherical_3D(fname + '_0501', radatadir)
    p_3d = Spherical_3D(fname + '_0502', radatadir)

    s = s_3d.vals
    p = p_3d.vals

    # p
    print ("Interpolating p")
    start = time.time()
    data = np.copy(p[:, ::-1, ::-1])
    data = np.append(data[:,:,:],np.expand_dims(data[0,:,:], axis=0), axis=0)
    F_interp = RegularGridInterpolator(original_coords_flat, data,\
            bounds_error=False, fill_value=0.0)
    p_cart = F_interp(points).reshape((n_x, n_y, n_z))
    end = time.time()
    print ("Interpolating for p took %.1f sec" %(end - start))

    # s
    print ("Interpolating s")
    start = time.time()
    data = np.copy(s[:, ::-1, ::-1])
    data = np.append(data[:,:,:],np.expand_dims(data[0,:,:], axis=0), axis=0)
    F_interp = RegularGridInterpolator(original_coords_flat, data,\
            bounds_error=False, fill_value=0.0)
    s_cart = F_interp(points).reshape((n_x, n_y, n_z))
    end = time.time()
    print ("Interpolating for s took %.1f sec" %(end - start))

    # p_prime
    print ("Interpolating p_prime")
    start = time.time()
    data = np.copy(p[:, ::-1, ::-1])
    data = data - np.mean(data, axis=0).reshape((1, nt, nr))
    data = np.append(data[:,:,:],np.expand_dims(data[0,:,:], axis=0), axis=0)
    F_interp = RegularGridInterpolator(original_coords_flat, data,\
            bounds_error=False, fill_value=0.0)
    p_prime_cart = F_interp(points).reshape((n_x, n_y, n_z))
    end = time.time()
    print ("Interpolating for p_prime took %.1f sec" %(end - start))

    # s_prime
    print ("Interpolating s_prime")
    start = time.time()
    data = np.copy(s[:, ::-1, ::-1])
    data = data - np.mean(data, axis=0).reshape((1, nt, nr))
    data = np.append(data[:,:,:],np.expand_dims(data[0,:,:], axis=0), axis=0)
    F_interp = RegularGridInterpolator(original_coords_flat, data,\
            bounds_error=False, fill_value=0.0)
    s_prime_cart = F_interp(points).reshape((n_x, n_y, n_z))
    end = time.time()
    print ("Interpolating for s_prime took %.1f sec" %(end - start))

    # Write the grid info to dictionary and save it
    savefile = savedir + fname + '_thermo.pkl'
    print ('Saving thermo data at ' + savefile)
    f = open(savefile, 'wb')
    pickle.dump({'rhobar': rhobar_cart, 'pbar': pbar_cart, 'tbar': tbar_cart,\
            's': s_cart, 'p': p_cart, 's_prime': s_prime_cart,\
            'p_prime': p_prime_cart},\
        f, protocol=4)
    f.close()

if interpom:
    omr_3d = Spherical_3D(fname + '_0301', radatadir)
    omt_3d = Spherical_3D(fname + '_0302', radatadir)
    omp_3d = Spherical_3D(fname + '_0303', radatadir)
    omr = omr_3d.vals
    omt = omt_3d.vals
    omp = omp_3d.vals

    # omr
    print ("Interpolating omr")
    start = time.time()
    data = np.copy(omr[:, ::-1, ::-1])
    data = np.append(data[:,:,:],np.expand_dims(data[0,:,:], axis=0), axis=0)
    F_interp = RegularGridInterpolator(original_coords_flat, data,\
            bounds_error=False, fill_value=0.0)
    omr_cart = F_interp(points).reshape((n_x, n_y, n_z))
    end = time.time()
    print ("Interpolating for omr took %.1f sec" %(end - start))

    # omt
    print ("Interpolating omt")
    start = time.time()
    data = np.copy(omt[:, ::-1, ::-1])
    data = np.append(data[:,:,:],np.expand_dims(data[0,:,:], axis=0), axis=0)
    F_interp = RegularGridInterpolator(original_coords_flat, data,\
            bounds_error=False, fill_value=0.0)
    omt_cart = F_interp(points).reshape((n_x, n_y, n_z))
    end = time.time()
    print ("Interpolating for omt took %.1f sec" %(end - start))

    # omp
    print ("Interpolating omp")
    start = time.time()
    data = np.copy(omp[:, ::-1, ::-1])
    data = np.append(data[:,:,:],np.expand_dims(data[0,:,:], axis=0), axis=0)
    F_interp = RegularGridInterpolator(original_coords_flat, data,\
            bounds_error=False, fill_value=0.0)
    omp_cart = F_interp(points).reshape((n_x, n_y, n_z))
    end = time.time()
    print ("Interpolating for omp took %.1f sec" %(end - start))

    # omr_prime
    print ("Interpolating omr_prime")
    start = time.time()
    data = np.copy(omr[:, ::-1, ::-1])
    data = data - np.mean(data, axis=0).reshape((1, nt, nr))
    data = np.append(data[:,:,:],np.expand_dims(data[0,:,:], axis=0), axis=0)
    F_interp = RegularGridInterpolator(original_coords_flat, data,\
            bounds_error=False, fill_value=0.0)
    omr_prime_cart = F_interp(points).reshape((n_x, n_y, n_z))
    end = time.time()
    print ("Interpolating for vr_prime took %.1f sec" %(end - start))

    # vt
    print ("Interpolating vt_prime")
    start = time.time()
    data = np.copy(vt[:, ::-1, ::-1])
    data = data - np.mean(data, axis=0).reshape((1, nt, nr))
    data = np.append(data[:,:,:],np.expand_dims(data[0,:,:], axis=0), axis=0)
    F_interp = RegularGridInterpolator(original_coords_flat, data,\
            bounds_error=False, fill_value=0.0)
    vt_prime_cart = F_interp(points).reshape((n_x, n_y, n_z))
    end = time.time()
    print ("Interpolating for vt_prime took %.1f sec" %(end - start))

    # vp_prime
    print ("Interpolating vp_prime")
    start = time.time()
    data = np.copy(vp[:, ::-1, ::-1])
    data = data - np.mean(data, axis=0).reshape((1, nt, nr))
    data = np.append(data[:,:,:],np.expand_dims(data[0,:,:], axis=0), axis=0)
    F_interp = RegularGridInterpolator(original_coords_flat, data,\
            bounds_error=False, fill_value=0.0)
    vp_prime_cart = F_interp(points).reshape((n_x, n_y, n_z))
    end = time.time()
    print ("Interpolating for vp_prime took %.1f sec" %(end - start))

    # Write the grid info to dictionary and save it
    savefile = savedir + fname + '_v.pkl'
    print ('Saving velocity data at ' + savefile)
    f = open(savefile, 'wb')
    pickle.dump({'vr': vr_cart, 'vt': vt_cart, 'vp': vp_cart,\
            'vr_prime': vr_prime_cart, 'vt_prime': vt_prime_cart,\
            'vp_prime': vp_prime_cart}, f, protocol=4)
    f.close()

    # Compute the Cartesian velocities
    print ("Computing Cartesian interpolated velocities, vx, vy, vz")
    vx = -(vr_cart*np.sin(theta_u) + vt_cart*np.cos(theta_u))*np.cos(phi_u) +\
            vp_cart*np.sin(phi_u)
    vy = -(vr_cart*np.sin(theta_u) + vt_cart*np.cos(theta_u))*np.sin(phi_u) -\
            vp_cart*np.cos(phi_u)
    vz = vr_cart*np.cos(theta_u)  - vt_cart*np.sin(theta_u)

    vx_prime = -(vr_prime_cart*np.sin(theta_u) +\
            vt_prime_cart*np.cos(theta_u))*np.cos(phi_u) +\
            vp_prime_cart*np.sin(phi_u)
    vy_prime = -(vr_prime_cart*np.sin(theta_u) +\
            vt_prime_cart*np.cos(theta_u))*np.sin(phi_u) -\
            vp_prime_cart*np.cos(phi_u)
    vz_prime = vr_prime_cart*np.cos(theta_u)  - vt_prime_cart*np.sin(theta_u)

    # Write the Cartesian velocity data to dictionary and save it
    savefile = savedir + fname + '_vcart.pkl'
    print ('Saving velocity data at ' + savefile)
    f = open(savefile, 'wb')
    pickle.dump({'vx': vx, 'vy': vy, 'vz': vz,\
        'vx_prime': vx_prime, 'vy_prime': vy_prime, 'vz_prime': vz_prime},\
        f, protocol=4)
    f.close()

mag = get_parameter(dirname, 'magnetism')
if mag:
    if interpmag:
        br_3d = Spherical_3D(fname + '_0001', radatadir)
        bt_3d = Spherical_3D(fname + '_0002', radatadir)
        bp_3d = Spherical_3D(fname + '_0003', radatadir)
        br = br_3d.vals
        bt = bt_3d.vals
        bp = bp_3d.vals

        # br
        print ("Interpolating br")
        start = time.time()
        data = np.copy(br[:, ::-1, ::-1])
        data = np.append(data[:,:,:],np.expand_dims(data[0,:,:], axis=0),\
                axis=0)
        F_interp = RegularGridInterpolator(original_coords_flat, data,\
                bounds_error=False, fill_value=0.0)
        br_cart = F_interp(points).reshape((n_x, n_y, n_z))
        end = time.time()
        print ("Interpolating for br took %.1f sec" %(end - start))

        # bt
        print ("Interpolating bt")
        start = time.time()
        data = np.copy(bt[:, ::-1, ::-1])
        data = np.append(data[:,:,:],np.expand_dims(data[0,:,:], axis=0),\
                axis=0)
        F_interp = RegularGridInterpolator(original_coords_flat, data,\
                bounds_error=False, fill_value=0.0)
        bt_cart = F_interp(points).reshape((n_x, n_y, n_z))
        end = time.time()
        print ("Interpolating for bt took %.1f sec" %(end - start))

        # vp
        print ("Interpolating bp")
        start = time.time()
        data = np.copy(bp[:, ::-1, ::-1])
        data = np.append(data[:,:,:],np.expand_dims(data[0,:,:], axis=0),\
                axis=0)
        F_interp = RegularGridInterpolator(original_coords_flat, data,\
                bounds_error=False, fill_value=0.0)
        bp_cart = F_interp(points).reshape((n_x, n_y, n_z))
        end = time.time()
        print ("Interpolating for bp took %.1f sec" %(end - start))

        # br_prime
        print ("Interpolating br_prime")
        start = time.time()
        data = np.copy(br[:, ::-1, ::-1])
        data = data - np.mean(data, axis=0).reshape((1, nt, nr))
        data = np.append(data[:,:,:],np.expand_dims(data[0,:,:], axis=0),\
                axis=0)
        F_interp = RegularGridInterpolator(original_coords_flat, data,\
                bounds_error=False, fill_value=0.0)
        br_prime_cart = F_interp(points).reshape((n_x, n_y, n_z))
        end = time.time()
        print ("Interpolating for br_prime took %.1f sec" %(end - start))

        # bt_prime
        print ("Interpolating vt_prime")
        start = time.time()
        data = np.copy(bt[:, ::-1, ::-1])
        data = data - np.mean(data, axis=0).reshape((1, nt, nr))
        data = np.append(data[:,:,:],np.expand_dims(data[0,:,:], axis=0),\
                axis=0)
        F_interp = RegularGridInterpolator(original_coords_flat, data,\
                bounds_error=False, fill_value=0.0)
        bt_prime_cart = F_interp(points).reshape((n_x, n_y, n_z))
        end = time.time()
        print ("Interpolating for bt_prime took %.1f sec" %(end - start))

        # bp_prime
        print ("Interpolating bp_prime")
        start = time.time()
        data = np.copy(bp[:, ::-1, ::-1])
        data = data - np.mean(data, axis=0).reshape((1, nt, nr))
        data = np.append(data[:,:,:],np.expand_dims(data[0,:,:], axis=0),\
                axis=0)
        F_interp = RegularGridInterpolator(original_coords_flat, data,\
                bounds_error=False, fill_value=0.0)
        bp_prime_cart = F_interp(points).reshape((n_x, n_y, n_z))
        end = time.time()
        print ("Interpolating for bp_prime took %.1f sec" %(end - start))

        # Write the magnetic fields to dictionary and save it
        savefile = savedir + fname + '_B.pkl'
        print ('Saving magnetic field data at ' + savefile)
        f = open(savefile, 'wb')
        pickle.dump({'br': br_cart, 'bt': bt_cart, 'bp': bp_cart,\
                'br_prime': br_prime_cart, 'bt_prime': bt_prime_cart,\
                'bp_prime': bp_prime_cart}, f, protocol=4)
        f.close()

        # Compute the Cartesian velocities
        print ("Computing Cartesian interpolated magnetic fields, bx, by, bz")
        bx = -(br_cart*np.sin(theta_u) +\
                bt_cart*np.cos(theta_u))*np.cos(phi_u) + bp_cart*np.sin(phi_u)
        by = -(br_cart*np.sin(theta_u) +\
                bt_cart*np.cos(theta_u))*np.sin(phi_u) - bp_cart*np.cos(phi_u)
        bz = br_cart*np.cos(theta_u)  - bt_cart*np.sin(theta_u)

        bx_prime = -(br_prime_cart*np.sin(theta_u) +\
                bt_prime_cart*np.cos(theta_u))*np.cos(phi_u) +\
                bp_prime_cart*np.sin(phi_u)
        by_prime = -(br_prime_cart*np.sin(theta_u) +\
                bt_prime_cart*np.cos(theta_u))*np.sin(phi_u) -\
                bp_prime_cart*np.cos(phi_u)
        bz_prime = br_prime_cart*np.cos(theta_u) - bt_prime_cart*np.sin(theta_u)

        # Write the Cartesian velocity data to dictionary and save it
        savefile = savedir + fname + '_Bcart.pkl'
        print ('Saving Cartesian B field data at ' + savefile)
        f = open(savefile, 'wb')
        pickle.dump({'bx': bx, 'by': by, 'bz': bz,\
            'bx_prime': bx_prime, 'by_prime': by_prime, 'bz_prime': bz_prime},\
            f, protocol=4)
        f.close()

end_overall = time.time()
t_overall = end_overall - start_overall
time_mins = int(np.floor(t_overall/60.))
time_sec = int(t_overall - time_mins*60.)
print ("Interpolation took total %i min %i sec" %(time_mins, time_sec))
