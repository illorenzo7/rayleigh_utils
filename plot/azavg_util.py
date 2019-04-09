###########################################################################
# This file contains utilities useful for plotting the data from azimuthal
# averages. (Files in the directory AZ_Avgs). 
# Written by Nick Featherstone
# First modified by Loren Matilsky, 03/09/2018

import numpy as np
import matplotlib as mpl
mpl.use('TkAgg')
import matplotlib.pyplot as plt
plt.rcParams['mathtext.fontset'] = 'dejavuserif'
csfont = {'fontname':'DejaVu Serif'}
plt.rcParams['contour.negative_linestyle'] = 'solid'

def fmt(x, pos):
    a, b = '{:.1e}'.format(x).split('e')
    b = int(b)
    return r'${} \times 10^{{{}}}$'.format(a, b)

def get_lims(arr, boundstype='minmax', caller_minmax=(-10.,10.)): 
                    # (Note the (-10, 10) are just placeholders so the
                    # calling function doesn't need to pass a required
                    # min/max tuple that is never used

    min_val, max_val = np.min(arr), np.max(arr)
    if (boundstype == 'manual'):
        min_val, max_val = caller_minmax
    return min_val, max_val

def plot_azav(fig, axis, field, radius, costheta, sintheta,
        mycmap=plt.cm.RdYlBu_r, units = r'$\rm{m}\ \rm{s}^{-1}$', 
        boundstype = 'minmax', nlevs = 10, caller_minmax=None, 
        plotcontours=True, plotfield=True,
        norm=None, levels=None):
    '''Takes a figure with a subplot (axis) of aspect ratio 1x2 and adds
    a plot in the meridional plane to the axis, with colorbar in the "cavity"
    of the meridional plane'''

    if (boundstype == 'minmax'):
        mini, maxi = get_lims(field)
    elif (boundstype == 'manual'):
        mini, maxi = caller_minmax

    # Get the exponent to use for scientific notation
    extent = np.max((np.abs(mini), np.abs(maxi)))
    exp = int(np.floor(np.log10(extent)))
    divisor = 10**exp
    
    # Normalize field by divisor
    field /= divisor
    mini /= divisor
    maxi /= divisor
   
    # Get the position of the axes on the figure
    pos = axis.get_position().get_points()
    axis_left, axis_bottom = pos[0]
    axis_right, axis_top = pos[1]
    axis_width = axis_right - axis_left
    axis_height = axis_top - axis_bottom
    axis_aspect = axis_height/axis_width
   
    # Set the colorbar axis to be in the "cavity" of the meridional plane
    # The colorbar height is set by making sure it "fits" in the cavity
    chi = np.min(radius)/np.max(radius)
    cavity_height = axis_height*chi
    cbaxis_center_x = axis_left + 0.3*axis_width
    cbaxis_center_y = axis_bottom + axis_height/2
    cbaxis_aspect = 10
    cbaxis_height = 0.5*cavity_height
    cbaxis_width = cbaxis_height/cbaxis_aspect / axis_aspect
    
    cbaxis_left = cbaxis_center_x - cbaxis_width/2
    cbaxis_bottom = cbaxis_center_y - cbaxis_height/2
    
    #Modified version of Antoine Strukarek's routine
    r = radius/np.max(radius)
    n_r = len(r)
    n_t = len(costheta)
    rtmp = r.reshape(1, n_r)
    cthtmp = costheta.reshape(n_t, 1)
    sthtmp = sintheta.reshape(n_t, 1)
    xr = np.dot(cthtmp, rtmp)
    yr = np.dot(sthtmp, rtmp)
    
    # Specify linewidths to be used in the meridional plane, one for the 
    # boundary (lw) and one for the contours (contour_lw)
    lw = 1
    contour_lw = .2
    
    if (plotfield):
        plt.sca(axis)
        plt.pcolormesh(yr,xr,field, vmin=mini, vmax=maxi,\
                cmap=mycmap, norm=norm)
     
        cbaxes = fig.add_axes([cbaxis_left, cbaxis_bottom,\
                       cbaxis_width, cbaxis_height])
        cbar = plt.colorbar(cax=cbaxes)

        fsize = 8 # fontsize for colorbar ticks and labels
        cbaxes.tick_params(labelsize=fsize)
#        cbar.set_label(units, rotation=270, labelpad=25, fontsize=18)
        cbar.ax.tick_params(labelsize=fsize)   #font size for the ticks
        ticks = np.array([mini, 0, maxi])
        ticklabels = []
        for i in range(len(ticks)):
            ticklabels.append(str(round(ticks[i],1)))
        ticks = np.array(ticks)
        cbar.set_ticks(ticks)
        cbar.set_ticklabels(ticklabels)

        # Put the units and exponent to left of colorbar
        cbar_label = (r'$10^{%i}\ $' %exp) + units
        fig.text(cbaxis_left - 0.3*cbaxis_width, cbaxis_center_y, cbar_label,\
            ha='right', va='center', rotation=90, fontsize=fsize)

    # Plot the boundary of the meridional plane
    plt.sca(axis)
    plt.plot(r[0]*sintheta,r[0]*costheta,'k', linewidth=lw)
    plt.plot(r[n_r-1]*sintheta,r[n_r-1]*costheta,'k', linewidth=lw)
    plt.plot([0,0], [r[n_r-1],r[0]], 'k', linewidth=lw)
    plt.plot([0,0], [-r[n_r-1],-r[0]], 'k', linewidth=lw)

    # Set axis ranges to be just outside the boundary lines
    lilbit = 0.01
    axis.set_xlim((-lilbit, 1 + lilbit))
    axis.set_ylim((-1 - lilbit, 1 + lilbit))
    axis.axis('off') 

    if (plotcontours):
        if levels is None:
            levs=mini+np.linspace(1,nlevs,nlevs)/float(nlevs)*(maxi-mini)
        else: # the caller specified specific contour levels to plot!
            levs=np.array(levels)
        plt.contour(yr,xr,field,colors='k',levels=levs, linewidths=contour_lw)

def streamfunction(vr,vt,r,cost,order=0):
    """------------------------------------------------------------
    This routine takes as input a divergenceless axisymmetric 
    vector field in spherical coordinates and computes from 
    it a streamfunction (a.k.a. a flux flunction).  The grid
    is decribed by r and costheta and can be non-uniform.
   ------------------------------------------------------------
    INPUTS:
   
    Vr, Vtheta = the 2-d vector velocity (or magnetic) field.
                 Dimensions are (N_Theta,N_R)
    r,cost     = the radius and cos(colatitude) of the grid.
                 r is assumed to vary from rmax to rmin and 
                 costheta from  1 to -1 (i.e. 90 degrees
                 to -90 degrees in latitude).
                 Dimensions are r(N_R), costheta(N_Theta)
    order      = If greater than zero, integration begins at the
                 outer shell and the north pole and proceeds
                 inward and southward.  If less than zero,
                 integration begins at the inner shell and 
                 south pole and proceeds upward and northward.
                 If equal to zero, both are done and an average
                 is taken.
   ------------------------------------------------------------
    OUTPUTS:
   
    psi = the streamfunction
   ------------------------------------------------------------
    """

    n_t,n_r=np.shape(vr)
    dtheta = np.zeros(n_t)
    dr     = np.zeros(n_r)

    psi = np.zeros((n_t,n_r))

    dpsi_dr = np.zeros((n_t,n_r))
    dpsi_dt = np.zeros((n_t,n_r))

    theta = np.arccos(cost)
    sint  = np.sqrt(1.0-cost**2)

    for i in range(n_t):
        dpsi_dr[i,:] = -r*sint[i]*vt[i,:]
        dpsi_dt[i,:] = r*r*sint[i]*vr[i,:]

    if (order >= 0):
        # double precision accumulation
        dtheta[1:n_t] = theta[1:n_t]-theta[0:n_t-1]
        dr[1:n_r] = r[1:n_r]-r[0:n_r-1]

        dtheta[0]=0 
        dr[0]=0

        for i in range(1,n_r):
            psi[1:n_t,i] = psi[1:n_t,i-1] + dpsi_dr[1:n_t,i]*dr[i]
        for i in range(1,n_t):
            psi[i,1:n_r] = psi[i-1,1:n_r] + dpsi_dt[i,1:n_r]*dtheta[i]

    if (order <= 0):
        psi2=np.zeros((n_t,n_r))
        
        dtheta[0:n_t-1] = theta[0:n_t-1]-theta[1:n_t]
        dr[0:n_r-1] = r[0:n_r-1]-r[1:n_r]
        
        dtheta[n_t-1]=0 
        dr[n_r-1]=0
        
        for i in range (n_r-2, -1, -1): 
            psi[0:n_t-1,i] = psi[0:n_t-1,i+1] + dpsi_dr[0:n_t-1,i]*dr[i]
        for i in range(n_t-2, -1, -1):
            psi[i,0:n_r-1] = psi[i+1,0:n_r-1] + dpsi_dt[i,0:n_r-1]*dtheta[i]
        
        if (order < 0):
            return psi2
        else:
            psi=0.5*(psi+psi2)
            
    return psi
