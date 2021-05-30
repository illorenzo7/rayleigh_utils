# FULL ROUTINES BY CONNOR BICE (SENT TO LOREN + BRAD 09/04/2020)

import numpy as np
import matplotlib
from matplotlib import ticker

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

def mollweide(fig,ax,data,r_scale=1,cmap='RdYlBu_r',bounds=None,boundstype='std',boundsfactor=2,cbar=True,extcbar=False):
    #fig and ax should be the matplotlib figure and axes objects on which to plot
    #data is what is plotted, should have dimensions longitude x latitude
    #r_scale allows the option to scale the figure inward from an outer radial boundary
    #cmap is the desired colormap, should be a string rather than an actual colormap object
    #bounds allows for direct control of the endpoints of the colorbar. If not specified, ends
    # will be determined by boundstype, either std (boundsfactor * standard deviation) or minmax
    #cbar determines whether or not a colorbar will be plotted
    #extcbar is a flag to return the img object, allowing the user to control the colorbar more directly

    pix = fig.get_figwidth()*fig.get_dpi()
    if not bounds is None:
        minabs=bounds[0]
        maxabs=bounds[1]
    elif boundstype == 'std':
        maxabs = np.std(data)*boundsfactor
        minabs = -maxabs
    elif boundstype == 'minmax':
        maxabs = np.max(data)
        minabs = np.min(data)
    else:
        maxabs = 1
        minabs = -1
    levels = np.linspace(minabs,maxabs,101)
    norm = matplotlib.colors.Normalize(vmin=minabs,vmax=maxabs)
    cmap = matplotlib.cm.get_cmap(cmap,256)
    
    data = np.append(data,data[[0],:],axis=0)
    ntheta = len(data[0,:])
    nphi = len(data[:,0])

    for i in range(ntheta):
        for j in range(nphi):
            if (data[j,i] > maxabs): data[j,i] = maxabs
            elif (data[j,i] < minabs): data[j,i] = minabs

    lats = np.linspace(-np.pi/2,np.pi/2,ntheta)
    lons = np.linspace(0,2*np.pi-1e-6,nphi)
    LATS,LONS = np.meshgrid(lats,lons)
    xx,yy = mollweide_transform(1,LATS,LONS)
    img = ax.contourf(xx*r_scale,yy*r_scale,data,cmap=cmap,levels=levels,norm=norm)

    for lat in np.pi/6*np.arange(-2,3,1): #Plot some latitudes
        x,y = mollweide_transform(r_scale,np.zeros(101)+lat,np.linspace(0,2*np.pi-1e-6,101))
        ax.plot(x,y,'k:',linewidth=pix/1024)
    for lon in [np.pi/2,np.pi,3*np.pi/2]: #Plot some meridians
        x,y = mollweide_transform(r_scale,np.linspace(-np.pi/2,np.pi/2,101),np.zeros(101)+lon)
        ax.plot(x,y,'k:',linewidth=pix/1024)
    for rad in np.unique([r_scale,1]): #Plot the outline
        x,y = mollweide_transform(rad,np.pi/2*np.append(np.append(np.linspace(-1,1,101),np.linspace(1,-1,101)),-1),np.append(np.append(np.zeros(101),np.zeros(101)+2*np.pi-1e-6),0))
        ax.plot(x,y,'k',linewidth=pix/1024)
                    
    ax.axis('equal')
    ax.axis('off')
    if cbar: 
        cb = fig.colorbar(img, orientation='horizontal', shrink=0.25, aspect = 20, pad = 0.05)
        cb.ax.tick_params(labelsize = 8*pix/1024, width=0.5*pix/1024, length=3*pix/1024)
        tick_locator = ticker.MaxNLocator(nbins=5)
        cb.locator = tick_locator
        cb.update_ticks()
    if extcbar: return img

def orthographic(fig,ax,data,r_scale=1,cmap='RdYlBu_r',bounds=None,boundstype='std',boundsfactor=2,az=0,el=0,cbar=True,extcbar=False):
    #fig and ax should be the matplotlib figure and axes objects on which to plot
    #data is what is plotted, should have dimensions longitude x latitude
    #r_scale allows the option to scale the figure inward from an outer radial boundary
    #cmap is the desired colormap, should be a string rather than an actual colormap object
    #bounds allows for direct control of the endpoints of the colorbar. If not specified, ends
    # will be determined by boundstype, either std (boundsfactor * standard deviation) or minmax
    #az and el together determine the central longitude and lattitude to view from, respectively
    #cbar determines whether or not a colorbar will be plotted
    #extcbar is a flag to return the img object, allowing the user to control the colorbar more directly

    pix = fig.get_figwidth()*fig.get_dpi()
    az = az * np.pi/180
    el = el * np.pi/180
    if not bounds is None:
        minabs=bounds[0]
        maxabs=bounds[1]
    elif boundstype == 'std':
        maxabs = np.std(data)*boundsfactor
        minabs = -maxabs
    elif boundstype == 'minmax':
        maxabs = np.max(data)
        minabs = np.min(data)
    else:
        maxabs = 1
        minabs = -1
    levels = np.linspace(minabs,maxabs,101)
    norm = matplotlib.colors.Normalize(vmin=minabs,vmax=maxabs)
    cmap = matplotlib.cm.get_cmap(cmap,256)

    data = np.append(data,data[[0],:],axis=0)
    ntheta = len(data[0,:])
    nphi = len(data[:,0])

    lats = np.linspace(-np.pi/2,np.pi/2,ntheta)
    lons = np.linspace(0,2*np.pi-1e-6,nphi)
    LATS,LONS = np.meshgrid(lats,lons)
    xx,yy,idx = ortho_transform(1,LATS,LONS,el,az)
    
    for i in range(ntheta):
        for j in range(nphi):
            if (data[j,i] > maxabs): data[j,i] = maxabs
            elif (data[j,i] < minabs): data[j,i] = minabs
            thisidx = np.where(idx[0]==j)
            if not i in idx[1][thisidx]: data[j,i] = np.nan 

    img = ax.contourf(xx*r_scale,yy*r_scale,data,cmap=cmap,levels=levels,norm=norm)

    for lat in np.pi/6*np.arange(-2,3,1): #Plot some latitudes
        x,y,i = ortho_transform(r_scale,np.zeros(101)+lat,np.linspace(az-np.pi,az+np.pi,101),el,az)
        ax.plot(x[i],y[i],'k:',linewidth=pix/1024)
    for lon in np.pi/6*np.arange(0,12,1): #Plot some meridians
        x,y,i = ortho_transform(r_scale,np.linspace(-np.pi/2,np.pi/2,101),np.zeros(101)+lon,el,az)
        ax.plot(x[i],y[i],'k:',linewidth=pix/1024)
    for rad in np.unique([r_scale,1]): #Plot the outline
        ax.plot(rad*np.cos(lons),rad*np.sin(lons),'k',linewidth=pix/1024)

    ax.axis('equal')
    ax.axis('off')
    if cbar: 
        cb = fig.colorbar(img, orientation='horizontal', shrink=0.25, aspect = 20, pad = 0.05)
        cb.ax.tick_params(labelsize = 8*pix/1024, width=0.5*pix/1024, length=3*pix/1024)
        tick_locator = ticker.MaxNLocator(nbins=5)
        cb.locator = tick_locator
        cb.update_ticks()
    if extcbar: return img

def perspective(fig,data,cmap='RdYlBu_r',bounds=None,boundstype='std',boundsfactor=2,az=0,el=0,bgcol='w',bumpmap=None,usamp=1):
    #fig should be the matplotlib figure on which to plot. Does not support passed Axes
    #data is what is plotted, should have dimensions longitude x latitude
    #cmap is the desired colormap, should be a string rather than an actual colormap object
    #bounds allows for direct control of the endpoints of the colorbar. If not specified, ends
    # will be determined by boundstype, either std (boundsfactor * standard deviation) or minmax
    #az and el together determine the central longitude and lattitude to view from, respectively
    #bgcol is a matplotlib color character for the background
    #bumpmap encodes small multiplicative shifts to the plotting radius, and should have the same dimensions as data
    # e.g. 1+0.03*vr/max(abs(vr))
    #usamp is an integer undersample rate for plotting speed in high-res slices. (N)^2 -> (N / usamp)^2

    from mpl_toolkits.mplot3d import Axes3D
    if not bounds is None:
        minabs=bounds[0]
        maxabs=bounds[1]
    elif boundstype == 'std':
        maxabs = np.std(data)*boundsfactor
        minabs = -maxabs
    elif boundstype == 'minmax':
        maxabs = np.max(data)
        minabs = np.min(data)
    else:
        maxabs = 1
        minabs = -1
    norm = matplotlib.colors.Normalize(vmin=minabs,vmax=maxabs)
    cmap = matplotlib.cm.get_cmap(cmap,256)

    ntheta = len(data[0,:])
    nphi = len(data[:,0]) 
    lats = np.linspace(-np.pi/2,np.pi/2,ntheta)[::usamp]
    lons = np.append(np.linspace(0,2*np.pi,nphi,endpoint=False)[::usamp],2*np.pi)
    LATS,LONS = np.meshgrid(lats,lons)
    X = np.cos(LONS)*np.cos(LATS)
    Y = np.sin(LONS)*np.cos(LATS)
    Z = np.sin(LATS)
    if bumpmap is None: R = np.ones_like(X)
    else: R = np.append(bumpmap[::usamp,::usamp],bumpmap[[0],:],axis=0)

    data = data[::usamp,::usamp] 
    data = np.append(data,data[[0],:],axis=0)
    for i in range(ntheta):
        for j in range(nphi):
            if (data[j,i] > maxabs): data[j,i] = maxabs
            elif (data[j,i] < minabs): data[j,i] = minabs

    ax=fig.gca(projection='3d')
    ax.set_facecolor(bgcol)
    img = ax.plot_surface(R*X,R*Y,R*Z,facecolors=cmap(norm(data)),rstride=1,cstride=1,linewidth=0,antialiased=False, shade=True)
    ax.axis('equal')
    ax.axis('off')
    ax.view_init(el,az)
   

