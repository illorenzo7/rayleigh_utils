import numpy as np
import h5py
import sys, os
sys.path.append(os.environ['rapp'])
from rayleigh_diagnostics import Spherical_3D

snapshots = {'13200001': 3253.259, '13300001': 3270.732,\
        '13350001': 3382.610, '13450001': 3302.782,\
        '13550001': 3324.979}

snapshot_start = snapshots['13200001']

#fname = input('Filenumber: ')
dirname = sys.argv[1]
fname = os.listdir(dirname + '/Spherical_3D')[0][:-5]
print("Reading data for Spherical_3D/" + fname)
args = sys.argv[2:]
nargs = len(args)

justmag = False
Om_subtract_nHz = 0.0
for i in range(nargs):
    arg = args[i]
    if arg == '-mag':
        justmag = True
    elif arg == '-sub':
        Om_subtract_nHz = float(args[i+1])

try:
    snapshot_now = snapshots[fname]
    delta_Prot = snapshot_now - snapshot_start
except:
    delta_Prot = 0.0

br_3d = Spherical_3D(fname + '_0801')
bt_3d = Spherical_3D(fname + '_0802')
bp_3d = Spherical_3D(fname + '_0803')

br = br_3d.vals
bt = bt_3d.vals
bp = bp_3d.vals

nr = br_3d.nr
nt = br_3d.ntheta
nphi = br_3d.nphi

Om_subtract_Hz = Om_subtract_nHz*1e-9
Prot = 2*np.pi/8.61e-6 # rotation period in sec.
phi_deflection = ((delta_Prot*Prot*Om_subtract_Hz) % 1.0)*(2*np.pi)
nroll = int(phi_deflection/(2*np.pi)*nphi)
print ("Om_subtract (nHz): ", Om_subtract_nHz)
print ("Starting Prot: ", snapshot_start, " Current Prot: ",\
        snapshot_now)
print ("Delta Prot: ", delta_Prot)
print ("Longitude deflection: ", phi_deflection*180.0/np.pi)
print ("nroll: ", nroll)

if not justmag:
    vr_3d = Spherical_3D(fname + '_0001')
    vt_3d = Spherical_3D(fname + '_0002')
    vp_3d = Spherical_3D(fname + '_0003')

    s_3d = Spherical_3D(fname + '_0501')
    p_3d = Spherical_3D(fname + '_0502')

    vr = vr_3d.vals
    vt = vt_3d.vals
    vp = vp_3d.vals

    s = s_3d.vals
    p = p_3d.vals

# Deflect the fields through phi if want to get into particular
# rotating frame
br = np.roll(br, -nroll, axis=0)
bt = np.roll(bt, -nroll, axis=0)
bp = np.roll(bp, -nroll, axis=0)

if not justmag:
    vr = np.roll(vr, -nroll, axis=0)
    vt = np.roll(vt, -nroll, axis=0)
    vp = np.roll(vp, -nroll, axis=0)

    s = np.roll(p, -nroll, axis=0)
    p = np.roll(p, -nroll, axis=0)

print('Read dimensions (nphi,nth,nr)=({:d},{:d},{:d})'.format(nphi,nt,nr,))

r = br_3d.r
theta = br_3d.theta
phi = np.linspace(0,2*np.pi,nphi,endpoint=False)
phi = np.expand_dims(phi,axis=1)

f = h5py.File('vol_00001.h5','w')

hr = f.create_dataset('r',(1,nr,),dtype='f',data=r)
ht = f.create_dataset('theta',(nt,1,),dtype='f',data=theta)
hp = f.create_dataset('phi',(nphi,1,),dtype='f',data=phi)

if not justmag:
    hvr = f.create_dataset('vr',(nphi,nt,nr,),dtype='f',data=vr)
    hvth = f.create_dataset('vth',(nphi,nt,nr,),dtype='f',data=vt)
    hvph = f.create_dataset('vph',(nphi,nt,nr,),dtype='f',data=vp)
    print('Wrote velocity fields')

hBr = f.create_dataset('Br',(nphi,nt,nr,),dtype='f',data=br)
hBth = f.create_dataset('Bth',(nphi,nt,nr,),dtype='f',data=bt)
hBph = f.create_dataset('Bph',(nphi,nt,nr,),dtype='f',data=bp)
print('Wrote magnetic fields')

f.close()
