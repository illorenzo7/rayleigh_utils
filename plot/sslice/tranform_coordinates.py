# Projection transforms, stolen from Connor Bice (09/04/20)
# Stolen on: 05/30/21
# Stolen by: Loren Matilsky

import numpy as np

def mollweide_transform(r,lat,lon,lon0=0,precision=1e-3): 
    lon = np.mod(lon-lon0,2*np.pi)-np.pi
    i = np.where(np.abs(lat)<np.pi/2)
    new_ts = lat[i]
    old_ts = np.zeros_like(lat[i])
    while(np.max(np.abs(new_ts-old_ts))>precision):
        old_ts = new_ts
        new_ts = old_ts - (2*old_ts + np.sin(2*old_ts) - np.pi*np.sin(lat[i]))/(2+2*np.cos(2*old_ts))
    lat[i] = new_ts
    xs = 2/np.pi * r * lon * np.cos(lat)
    ys = r * np.sin(lat)
    return xs,ys

def ortho_transform(r,lat,lon,lat0=0,lon0=0):
    xs = r*np.cos(lat)*np.sin(lon-lon0)
    ys = r*(np.cos(lat0)*np.sin(lat)-np.sin(lat0)*np.cos(lat)*np.cos(lon-lon0))
    cosc = np.sin(lat0)*np.sin(lat)+np.cos(lat0)*np.cos(lat)*np.cos(lon-lon0) #cosine of angular distance from center of view
    idx = np.where(cosc>=0) #these indices are on front of the globe, the rest should be clipped.
    return xs,ys,idx
