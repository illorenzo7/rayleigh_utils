# Projection transforms, stolen from Connor Bice (09/04/20)
# Stolen on: 05/30/21
# Stolen by: Loren Matilsky

import numpy as np

def mollweide_transform(costheta, clon=0., shrinkage=1., precision=1.e-3): 
    # compute the spherical coordinates
    ntheta = len(costheta)
    nphi = 2*ntheta
    tt = np.arccos(costheta)
    lat = np.pi/2. - tt # these "latitudes" are in radians...
    lon = np.linspace(0., 2.*np.pi, nphi, endpoint=False)

    # shift the lons 
    lon = np.mod(lon-lon0,2*np.pi)-np.pi

    # make each array 2D

    #i = np.where(np.abs(lat)<np.pi/2)
    #print (
    #new_ts = lat[i]
    #old_ts = np.zeros_like(lat[i])
    #while(np.max(np.abs(new_ts-old_ts))>precision):
    #    old_ts = new_ts
    #    new_ts = old_ts - (2*old_ts + np.sin(2*old_ts) - np.pi*np.sin(lat[i]))/(2+2*np.cos(2*old_ts))
    #lat[i] = new_ts
    # compute the beta-angle for the projection 
    # (related to latitude by transcendental equation, which we need
    # to solve iteratively
    beta = np.zeros_like(lat)
    new_beta = np.copy(lat)
    while(np.max(np.abs(new_beta - beta)) > precision):
        beta = new_beta
        new_beta = beta - (2*beta + np.sin(2.*beta) -\
                np.pi*np.sin(lat))/(2.+2.*np.cos(2.*beta))
    # get a "meshgrid" from 1D arrays
    lon, beta = np.meshgrid(lon, beta, indexing='ij')
    xs = 2./np.pi*shrinkage*lon*np.cos(beta)
    ys = shrinkage*np.sin(beta)
    return xs, ys

def ortho_transform(r,lat,lon,lat0=0,lon0=0):
    xs = r*np.cos(lat)*np.sin(lon-lon0)
    ys = r*(np.cos(lat0)*np.sin(lat)-np.sin(lat0)*np.cos(lat)*np.cos(lon-lon0))
    cosc = np.sin(lat0)*np.sin(lat)+np.cos(lat0)*np.cos(lat)*np.cos(lon-lon0) #cosine of angular distance from center of view
    idx = np.where(cosc>=0) #these indices are on front of the globe, the rest should be clipped.
    return xs,ys,idx
