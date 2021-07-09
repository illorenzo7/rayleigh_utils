import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pylab as p 
import numpy as np
from cmd_util import *
import sys
from scipy.interpolate import RegularGridInterpolator as rgi
from numpy.random import rand
from diagnostic_reading import AzAverage, build_file_list, TimeAvg_AZAverages
from azavg_util import streamfunction

def help():
    print('plot_field_lines.py can (and should) be run with a number of options \n')
    print('--files=   A series of comma and/or colon separated integers which correspond to the desired iterations.\n  eg 100000,20000:10000:250000 \n  Default: Existing azavg file\n')
    print('--fname=   The name of the file to be written.\n  Default: mc_times.png\n')
    print('--rstar=   The radius of the star youre trying to model in cm.\n  Default: 2.588e10 \n')
    print('--rbcz=    The fractional radius corresponding to the bottom of your stars CZ, if relevant.\n  Default: 0 \n')
    print('--nlines=  The number of streamlines to calculate.\n  Default: 1000 \n')
    print('--rbnds=   A comma separated pair of fractional radii bounding the area to seed field lines in.\n  Default: full domain \n')
    print('--tbnds=   A comma separated pair of colatitudes in degrees bounding the area to seed field lines in.\n  Default: full domain \n')
    print('--smax=    The maximum streamline length permitted as a multiple of rstar.\n  Default: 2 pi\n')
    print('--sclose=  The distance required to consider a streamline closed as a multiple of rstar.\n Default: 1e-2\n')
    print('--ds=      The integration step length as a multiple of rstar.\n  Default: 1e-3\n')
    print('--nbins=   The number of bins to use in the histogram of circulation times.\n  Default: 50\n')
    print('--xlims=   A comma separated pair of time values in days to bound the xaxis with.\n  Default: Data min,max')
    print('--weightf= The weighting function to use for the distribution of streamline seeds.\n  Options: uniform, massflux\n  Default: uniform\n')
    print('--weightp= The power given to the weighting function, if non-uniform.\n  Default: 1\n')
    print('--snorm    A flag to weight the histogram by the inverse of the circulation length\n  Default: False')
    print('--help     Who fuckin knows when a code is this spaghetti?\n')
    sys.exit(0)

#Calculates a streamline of length s initiated at X0 for the vector field represented by the interpolating functions
def calcStreamline(X0,fnr,fnt,ds,smax,sclose):
    coords = np.zeros((2,1))
    coords[:,0]=X0
    time = 0
    started = False
    for k in range(int(np.ceil(smax/ds))):
        try:
            vr = fnr(coords[:,k])
            vt = fnt(coords[:,k])
        except ValueError:
            badc = coords[:,k]
            if badc[0]>np.max(theta) or badc[0]<np.min(theta): print('    Streamline went out of bounds (t={:.2f}) after {:d} iters'.format(badc[0],k))
            elif badc[1]>np.max(radius) or badc[1]<np.min(radius): print('    Streamline went out of bounds (r={:.2e}) after {:d} iters'.format(badc[1],k))
            else: print('    Streamline went out of bounds (r={:.2e},t={:.2f}) after {:d} iters'.format(badc[1],badc[0],k))
            return coords,0,0
        v = np.sqrt(vr**2+vt**2)
        time = time + ds/v
        dcoords = ds/v*np.array([vt/coords[1,-1],vr]).T
        new_coords = coords[:,-1]+dcoords
        coords = np.append(coords,new_coords.T,axis=1)
        dispx = coords[1,0]*np.cos(coords[0,0])-coords[1,-1]*np.cos(coords[0,-1])
        dispy = coords[1,0]*np.sin(coords[0,0])-coords[1,-1]*np.sin(coords[0,-1])
        dist = np.sqrt(dispx**2+dispy**2)
        if not started and dist>1.5*sclose: 
            started = True
        if started and dist<sclose:
            time = time + dist/v
            coords = np.append(coords,X0.reshape(2,1),axis=1)
            return coords,time,ds*(k+1)+dist
    print('    Streamline did not close and will be discarded')
    return coords,0,0


#Read and interpret all the arguments
args = sys.argv
opts = getOpt(args[1:],['fname=','files=','rstar=','rbcz=','nlines=','rbnds=','tbnds=','smax=','sclose=','ds=','nbins=','xlims=','weightf=','weightp=','snorm','help'])

if 'help' in opts: help()
if 'fname' in opts: fname = opts['fname']
else: fname = 'mc_times'
if fname[-4:] != '.png': fname = fname+'.png'
if 'files' in opts: 
    file_list = ['AZ_Avgs/'+convertNumber(int(x)) for x in parseList(opts['files'])]
    TimeAvg_AZAverages(file_list,'time_averaged_azavg')
else: print('No files chosen, using existing time averaged file')
if 'rstar' in opts: rstar = float(opts['rstar'])
else: rstar = 2.588e10
if 'rbcz' in opts: rbcz = float(opts['rbcz'])
else: rbcz = 0
if 'nlines' in opts: nlines = int(opts['nlines'])
else: nlines = 1000
if 'rbnds' in opts: rbnds = np.array([float(x) for x in opts['rbnds'].split(',')])*rstar
else: rbnds = [0,rstar]
if 'tbnds' in opts: tbnds = np.array([float(x) for x in opts['tbnds'].split(',')])*np.pi/180.
else: tbnds = [0,np.pi] 
if 'smax' in opts: smax = float(opts['smax'])*rstar
else: smax = 2*np.pi*rstar
if 'sclose' in opts: sclose = float(opts['sclose'])*rstar
else: sclose = 1e-2*rstar
if 'ds' in opts: ds = float(opts['ds'])*rstar
else: ds = 1e-3*rstar
if 'nbins' in opts: nbins = int(opts['nbins'])
else: nbins = 50
if 'xlims' in opts: xlims = [float(x) for x in opts['xlims'].split(',')]
else: xlims = None
if 'weightf' in opts: 
    weightf = opts['weightf']
    if not weightf in ['uniform','massflux']:
        print('Weighting function "{:s}" not recognized, using uniform'.format(weightf))
        weightf = 'uniform'
else: weightf = 'uniform'
if 'weightp' in opts: weightp = float(opts['weightp'])
else: weightp = 1.
if 'snorm' in opts: snorm = True
else: snorm = False

az = AzAverage(filename='time_averaged_azavg',path='')

# time index to grab (this file only has one time since it is an average)
ind = 0 

#Find the indices associated with the quantities we want to plot
vr_index   = az.lut[1]
vt_index   = az.lut[2]
rhovr_index   = az.lut[201]
rhovt_index   = az.lut[202]

#Grab quantities of interest from az.vals
n_r = az.nr
n_t = az.ntheta

vr = az.vals[:,:,vr_index,ind].reshape(n_t,n_r)[::-1,::-1]
vt = az.vals[:,:,vt_index,ind].reshape(n_t,n_r)[::-1,::-1]
rhovr = az.vals[:,:,rhovr_index,ind].reshape(n_t,n_r)[::-1,::-1]
rhovt = az.vals[:,:,rhovt_index,ind].reshape(n_t,n_r)[::-1,::-1]

radius = az.radius[::-1]
costheta = az.costheta[::-1]
theta = np.arccos(costheta)
try: 
    overlap_ind = np.where(radius[1:]==radius[:-1])[0][0]
    radius = np.append(radius[:overlap_ind],radius[overlap_ind+1:])
    vr = np.append(vr[:,:overlap_ind],vr[:,overlap_ind+1:],axis=1)
    vt = np.append(vt[:,:overlap_ind],vt[:,overlap_ind+1:],axis=1)
    rhovr = np.append(rhovr[:,:overlap_ind],rhovr[:,overlap_ind+1:],axis=1)
    rhovt = np.append(rhovt[:,:overlap_ind],rhovt[:,overlap_ind+1:],axis=1)
except IndexError: overlap_ind = None
psi = streamfunction(rhovr,rhovt,radius,costheta,order=0)

fnr = rgi((theta,radius),vr)
fnt = rgi((theta,radius),vt)
fnpsi = rgi((theta,radius),psi)

if rbnds[0]<np.min(radius): rbnds[0] = np.min(radius)
if rbnds[1]>np.max(radius): rbnds[1] = np.max(radius)
if tbnds[0]<np.min(theta): tbnds[0] = np.min(theta)
if tbnds[1]>np.max(theta): tbnds[1] = np.max(theta)

times = np.zeros(nlines)
lengths = np.zeros(nlines)

rout = np.max(radius)/rstar
rin = np.min(radius)/rstar
merid = np.linspace(-np.pi/2,np.pi/2,100)

f1 = p.figure(figsize=(18.,9.), dpi=300)
plt.rcParams.update({'font.size': 12})
if snorm: grid = plt.GridSpec(2,2)
else: grid = plt.GridSpec(1,3)
plt.subplot(grid[:,0])
plt.plot([0,0],[rin,rout],'k')
plt.plot([0,0],[-rin,-rout],'k')
plt.plot(rin*np.cos(merid),rin*np.sin(merid),'k')
plt.plot(rout*np.cos(merid),rout*np.sin(merid),'k')
plt.plot(rbcz*np.cos(merid),rbcz*np.sin(merid),'k--')

alpha = 100./nlines
if alpha>1: alpha = 1.

if weightf != 'uniform':
    if weightf == 'massflux': funcbase = np.sqrt(rhovr**2+rhovt**2)
    else: funcbase = np.ones_like(rhovr)
    func = rgi((theta,radius),funcbase**weightp)
    rs = np.linspace(rbnds[0],rbnds[1],100)
    ts = np.linspace(tbnds[0],tbnds[1],300)
    rrs,tts = np.meshgrid(0.5*(rs[1:]+rs[:-1]),0.5*(ts[1:]+ts[:-1]))
    drrs,dtts = np.meshgrid(rs[1:]-rs[:-1],ts[1:]-ts[:-1])
    dFs = func((tts,rrs))*rrs*drrs*dtts
    Frt = np.sum(np.sum(dFs))
    Ftn = np.cumsum(np.sum(dFs,axis=0))/Frt

for k in range(nlines):
    print('  Tracing line {:d}/{:d}'.format(k+1,nlines))
    if weightf == 'uniform': x = rand(2)*np.array([tbnds[1]-tbnds[0],rbnds[1]-rbnds[0]])+np.array([tbnds[0],rbnds[0]])
    else: 
        rind = np.argmin(np.abs(Ftn-rand(1)))
        Frn = np.cumsum(dFs[:,rind])/np.sum(dFs[:,rind])
        tind = np.argmin(np.abs(Frn-rand(1)))
        x = np.array([tts[tind,rind],rrs[tind,rind]])
    if fnpsi(x)>0: col = [1,0,0,alpha]
    else: col = [0,0,1,alpha]
    rs,times[k],lengths[k] = calcStreamline(x,fnr,fnt,ds,smax,sclose)
    if times[k]>0: plt.plot(rs[1,:]*np.cos(np.pi/2-rs[0,:])/rstar,rs[1,:]*np.sin(np.pi/2-rs[0,:])/rstar,color=col)

plt.xlim(0,1)
plt.ylim(-1,1)
plt.axis('equal')
plt.axis('off')
plt.title('Circulation Streamlines')

plt.subplot(grid[0,1:])
nztimes = times[np.where(times>0)[0]]/(24.*60.*60.)
timebins = np.logspace(np.log10(np.min(nztimes)),np.log10(np.max(nztimes)),nbins)
plt.hist(nztimes,bins=timebins,normed=True,weights=None)
if xlims != None: plt.xlim(xlims[0],xlims[1])
plt.gca().set_xscale('log')
plt.xlabel('Circulation time (days)')
plt.title('Unweighted Circulation Time Distribution')

if snorm:
    plt.subplot(grid[1,1:])
    histweights = 1./lengths[np.where(lengths>0)[0]]
    plt.hist(nztimes,bins=timebins,normed=True,weights=histweights)
    if xlims != None: plt.xlim(xlims[0],xlims[1])
    plt.gca().set_xscale('log')
    plt.xlabel('Circulation time (days)')
    plt.title('Weighted Circulation Time Distribution')

plt.tight_layout()
p.savefig(fname)

