"""
Reads in 3-D ball data and interpolates to a uniform cartesian grid for 3-D vis.

Usage:
    uniform_cartesian.py [options]

Options:
    --filename=<filename>   Name of file to read in [default: vol_00001.h5]
    --n=<n>                 n^3 resolution of uniform cartesian cube [default: 256]
    --svars=<svars>         Comma separated list of scalar variable names [default: ]
    --cartvvars=<cartvvars> Comma separated list of cartesian vector variable prefixes [default: ]
    --sphvvars=<sphvvars>   Comma separated list of spherical vector variable prefixes [default: ]
"""
from docopt import docopt
args = docopt(__doc__)

filename = args['--filename']
n = int(float(args['--n']))
svars = args['--svars']
sphvvars = args['--sphvvars']
cartvvars = args['--cartvvars']

import logging

import numpy as np
from scipy.spatial import Delaunay
from scipy.interpolate import LinearNDInterpolator
import h5py
import time
import netCDF4 as n4

logger = logging.getLogger(__name__)

#fields = args['--fields'].split(',')

f = h5py.File('{:s}'.format(filename), 'r')
r = f['r']/np.max(f['r'])
aspect = np.min(f['r'])/np.max(f['r'])
theta = f['theta']
phi = f['phi']

# [:,::-1,:] is to get theta into strictly increasing order
r = np.expand_dims(r, axis=0)[:,::-1,::-1]
theta = np.expand_dims(theta, axis=0)[:,::-1,::-1]
phi = np.expand_dims(phi, axis=1)[:,::-1,::-1]
# wrap phi points to ensure good interpolation properties in periodic direction
phi = np.append(phi[:,:,:],np.expand_dims(phi[0,:,:]+2*np.pi, axis=0), axis=0)

print("opening {}".format(filename))
print(" r:{}".format(r.shape))
print("th:{}".format(theta.shape))
print("ph:{}".format(phi.shape))

x = r*np.sin(theta)*np.cos(phi)
y = r*np.sin(theta)*np.sin(phi)
z = r*np.cos(theta)*np.ones_like(phi)

points = np.array([x.flatten(), y.flatten(), z.flatten()]).T
print("point cloud: {}".format(points.shape))

n_x = n_y = n_z = n
x_u = np.linspace(-1, 1, n_x)[:,None,None]
y_u = np.linspace(-1, 1, n_y)[None,:,None]
z_u = np.linspace(-1, 1, n_z)[None,None,:]

zero = np.zeros((n_x,n_y,n_z))
r_u= zero + (x_u**2 + y_u**2 + z_u**2)**0.5
theta_u = zero + np.arccos(z_u/r_u)
phi_u = zero + np.arctan2(y_u,x_u) + np.pi
#phi_u = zero + np.arctan2(y_u,x_u)
## Get phi value in range (0, 2*pi)
#phi_u[np.where(phi_u < 0.)] += (2*np.pi)

from scipy.interpolate import RegularGridInterpolator
original_coords_flat = (phi.flatten(), theta.flatten(), r.flatten())
new_coords_flat = (phi_u.flatten(), theta_u.flatten(), r_u.flatten())
points = np.array(list(zip(*new_coords_flat)))



dataset = n4.Dataset('vapor_field_data.nc', 'w', format='NETCDF4')
xset = dataset.createDimension('x', n)
yset = dataset.createDimension('y', n)
zset = dataset.createDimension('z', n)
tset = dataset.createDimension('t', None)
xs = dataset.createVariable('X', np.float64, ('z','y','x'))
ys = dataset.createVariable('Y', np.float64, ('z','y','x'))
zs = dataset.createVariable('Z', np.float64, ('z','y','x'))
#xs = dataset.createVariable('x', np.float64, ('x',))
#ys = dataset.createVariable('y', np.float64, ('y',))
#zs = dataset.createVariable('z', np.float64, ('z',))
ts = dataset.createVariable('t', np.float64, ('t',))
xs[:] = (x_u + zero).T
ys[:] = (y_u + zero).T
zs[:] = (z_u + zero).T
#xs[:] = x_u
#ys[:] = y_u
#zs[:] = z_u
ts[:] = 0
var_string = "X:Y:Z:"

data_dict = {}

if not (svars=='' or svars==None):
    for v in svars.split(','):
        print('Reading scalar variable {:s}'.format(v,))
        if v[-5:]=='_naxi':
            data = np.array(f[v[:-5]])[:,::-1,::-1]
            data = data - np.mean(data,axis=0,keepdims=True)
        else: data = np.array(f[v])[:,::-1,::-1]
        data = np.append(data[:,:,:],np.expand_dims(data[0,:,:], axis=0), axis=0)
        data_dict.update({v : data})
if not (sphvvars=='' or sphvvars==None):
    for v in sphvvars.split(','):
        print('Reading spherical vector variable {:s}'.format(v,))
        if v[-5:]=='_naxi':
            data = np.array(f[v[:-5]+'r'])[:,::-1,::-1]
            data = data - np.mean(data,axis=0,keepdims=True)
        else: data = np.array(f[v+'r'])[:,::-1,::-1]
        data_r = np.append(data[:,:,:],np.expand_dims(data[0,:,:], axis=0), axis=0)
        if v[-5:]=='_naxi':
            data = np.array(f[v[:-5]+'th'])[:,::-1,::-1]
            data = data - np.mean(data,axis=0,keepdims=True)
        else: data = np.array(f[v+'th'])[:,::-1,::-1]
        data_t = np.append(data[:,:,:],np.expand_dims(data[0,:,:], axis=0), axis=0)
        if v[-5:]=='_naxi':
            data = np.array(f[v[:-5]+'ph'])[:,::-1,::-1]
            data = data - np.mean(data,axis=0,keepdims=True)
        else: data = np.array(f[v+'ph'])[:,::-1,::-1]
        data_p = np.append(data[:,:,:],np.expand_dims(data[0,:,:], axis=0), axis=0)
        data_mag = np.sqrt(data_r**2 + data_t**2 + data_p**2)
        data_dict.update({v+'r':data_r,v+'th':data_t,v+'ph':data_p,v+'mag':data_mag})
if not (cartvvars=='' or cartvvars==None):
    for v in cartvvars.split(','):
        print('Reading cartesian vector variable {:s}'.format(v,))
        if v[-5:]=='_naxi':
            data = np.array(f[v[:-5]+'r'])[:,::-1,::-1]
            data = data - np.mean(data,axis=0,keepdims=True)
        else: data = np.array(f[v+'r'])[:,::-1,::-1]
        data_r = np.append(data[:,:,:],np.expand_dims(data[0,:,:], axis=0), axis=0)
        if v[-5:]=='_naxi':
            data = np.array(f[v[:-5]+'th'])[:,::-1,::-1]
            data = data - np.mean(data,axis=0,keepdims=True)
        else: data = np.array(f[v+'th'])[:,::-1,::-1]
        data_t = np.append(data[:,:,:],np.expand_dims(data[0,:,:], axis=0), axis=0)
        if v[-5:]=='_naxi':
            data = np.array(f[v[:-5]+'ph'])[:,::-1,::-1]
            data = data - np.mean(data,axis=0,keepdims=True)
        else: data = np.array(f[v+'ph'])[:,::-1,::-1]
        data_p = np.append(data[:,:,:],np.expand_dims(data[0,:,:], axis=0), axis=0)
        data_x = data_r*np.sin(theta)*np.cos(phi) + data_t*np.cos(theta)*np.cos(phi) - data_p*np.sin(phi)
        data_y = data_r*np.sin(theta)*np.sin(phi) + data_t*np.cos(theta)*np.sin(phi) + data_p*np.cos(phi)
        data_z = data_r*np.cos(theta) - data_t*np.sin(theta)
        data_mag = np.sqrt(data_r**2 + data_t**2 + data_p**2)
        data_dict.update({v+'x':data_x,v+'y':data_y,v+'z':data_z,v+'mag':data_mag})
    
for field in data_dict.keys():
    data = data_dict[field]
    print('Min '+field+' = ' +str(np.min(data)))
    print('Max '+field+' = ' +str(np.max(data)))
    print('StD '+field+' = ' +str(np.std(data)))
    F_interp = RegularGridInterpolator(original_coords_flat, data,
                                       bounds_error=False, fill_value=0)

    start_int = time.time()
    data_u = F_interp(points).reshape((n_x, n_y, n_z))
    data_u[np.where(r_u < aspect)] = 0.0 # Make sure data is zero outside shell boundaries
    data_u[np.where(r_u > 1.0)] = 0.0 # Make sure data is zero outside shell boundaries
    end_int = time.time()
    print("Interpolation took {:g} seconds for {}".format(end_int-start_int, field))
    datas = dataset.createVariable(field, np.float64, ('z','y','x'))
    datas[:] = data_u.T
    var_string += field + ':'

var_string = var_string[:-1]

dataset.close()

import subprocess
subprocess.call(["ncdfvdfcreate","-dims", "{:d}x{:d}x{:d}".format(n_x,n_y,n_z), "-vars", var_string, "-vars3d", var_string, "-extents", "-1:-1:-1:1:1:1", "vapor_field_data.nc", "vapor_field_data.vdf"])
subprocess.call(["ncdf2vdf", "-vars", var_string, "vapor_field_data.nc", "vapor_field_data.vdf"])
