# This file contains utilities useful for plotting the data from azimuthal
# averages,
# files in the directory AZ_Avgs.
# Written by Nick Featherstone
# First modified by Loren Matilsky, 03/09/2018

import numpy as np
import matplotlib as mpl
from matplotlib import ticker
mpl.use('TkAgg')
import matplotlib.pyplot as plt
from matplotlib import colors
plt.rcParams['mathtext.fontset'] = 'dejavuserif'
csfont = {'fontname':'DejaVu Serif'}
plt.rcParams['contour.negative_linestyle'] = 'solid'
import sys, os
sys.path.append(os.environ['rapl'])
sys.path.append(os.environ['raco'])
from common import *
from plotcommon import *

def plot_azav(field, rr, cost, rbcz=None, minmaxrz=None, linthreshrz=None, linscalerz=None, cmaprz=None, rvals=np.array([]), rstyles=[], rcolors=[], rlw=1.0, plotlatlines=True, latvals=np.array([]), latstyles=[], latcolors=[], latlw=1.0, **kwargs):
    # **kwargs corresponds to my_contourf

    # make copy of field
    field = np.copy(field)
    field_full = np.copy(field)

    # grid info
    nt, nr = len(cost), len(rr)
    zeros = np.zeros((nt, nr))
    rr_2d = rr.reshape((1, nr)) + zeros
    rmax = np.max(rr_2d)
    cost_2d = cost.reshape((nt, 1)) + zeros
    sint_2d = np.sqrt(1.0 - cost_2d**2.0)
    tt_lat = 180.0/np.pi*(np.pi/2.0 - np.arccos(cost_2d))
    # tt_lat is 2D here!

    # use these to plot lines/boundaries
    xx_full = rr_2d*sint_2d/rmax
    yy_full = rr_2d*cost_2d/rmax
    rr_full = rr_2d/rsun
    tt_lat_full = np.copy(tt_lat)

    if rbcz is None: # just plotting 1 domain
        xx = xx_full
        yy = yy_full
        field = field_full
        rr_cz = rr_full
    else: # plotting two domains
        irbcz = np.argmin(np.abs(rr/rsun - rbcz))
        fieldrz = field[:, irbcz+1:]
        field = field[:, :irbcz+1]
        xx = (rr_2d*sint_2d)[:, :irbcz+1]/ro
        yy = (rr_2d*cost_2d)[:, :irbcz+1]/ro
        xxrz = (rr_2d*sint_2d)[:, irbcz+1:]/ro
        yyrz = (rr_2d*cost_2d)[:, irbcz+1:]/ro
        rr_cz = rr_full[:, :irbcz+1]
        rr_rz = rr_full[:, irbcz+1:]
        tt_lat_rz = tt_lat[:, irbcz+1:]
        tt_lat = tt_lat[:, :irbcz+1]

    if not rbcz is None: # plot the RZ field first (before "showing"
        # the plot, potentially)
        # will need to change some kwargs:
        kwargsrz = dict(kwargs)
        kwargsrz['showplot'] = False
        kwargsrz['minmax'] = minmaxrz
        kwargsrz['linthresh'] = linthreshrz
        kwargsrz['linscale'] = linscalerz
        if cmaprz is None:
            if 'logscale' in kwargs.keys() or 'posdef' in kwargs.keys():
                if kwargs['logscale']:
                    cmaprz = 'Blues'
                elif kwargs['posdef']:
                    cmaprz = 'cividis'
                else: # this is the default, and symlog
                    cmaprz = 'PuOr_r'    
            else:
                cmaprz = 'PuOr_r'    
        kwargsrz['cmap'] = cmaprz
        out = my_contourf(xxrz, yyrz, fieldrz, func1=rr_rz, func2=tt_lat_rz, cbar_no=2, **kwargsrz)

    # regardless, plot the CZ field
    if not rbcz is None and (not 'fig' in kwargs.keys() or\
            not 'ax' in kwargs.keys()): # in this case, need to use the
        # already-generated fig, ax from the RZ plot
        kwargs['fig'], kwargs['ax'] = out
    fig, ax = my_contourf(xx, yy, field, func1=rr_cz, func2=tt_lat, **kwargs)

    # potentially plot coordinate lines
    if plotlatlines:
        if latvals == []:
            latvals = np.arange(-60.0, 89.0, 30.0)
    else:
        latvals = []

    my_contourf(xx_full, yy_full, field_full, fig=fig, ax=ax, plotfield=False, func1=rr_full, vals1=rvals, linestyles1=rstyles, colors1=rcolors, lw1=rlw, func2=tt_lat_full, vals2=latvals, linestyles2=latstyles, colors2=latcolors, lw2=latlw)

def plot_azav_half(field, rr, cost, sym='even', fig=None, ax=None,\
        cmap='RdYlBu_r', units='', minmax=None, posdef=False, logscale=False,\
        symlog=False, linthresh=None, linscale=None, plotcontours=True,\
        plotfield=True, nlevs=10, levels=None, plotlatlines=False, rvals=[],\
        norm=None, fontsize=10, showplot=False, plot_cbar=True):
	
    '''Takes a figure with a subplot (axis) of aspect ratio 1x1 (or
    generates default axes if they are not provided) and adds
    a plot of the upper meridional plane to the axis, averaging the full
    plane with either even or odd symmetry
    '''

    # Get default bounds if not specified
    # Do this before "folding" the field
    if minmax is None:
        field_cut = trim_field(np.copy(field), rr, cost)
        minmax = get_satvals(field_cut, posdef=posdef, logscale=logscale,\
                symlog=symlog)
        
    # Grid info -- make "folded" grid
    nr = len(rr)
    nt = len(cost)
    ri, ro = np.min(rr), np.max(rr)
    
    # "Fold" the field in half
    # Please don't try this with odd N_theta!
    it_half = int(nt/2)
    
    # "Half grid" quantities
    cost = cost[it_half:]
    lats = (np.pi/2. - np.arccos(cost))*180./np.pi
    sint = np.sqrt(1. - cost**2.)
    
    # Calculate the grid on which to plot
    rr_2d = rr.reshape((1, nr))
    cost_2d = cost.reshape((int(nt/2), 1))
    sint_2d = sint.reshape((int(nt/2), 1))
    xx = rr_2d*sint_2d/ro
    zz = rr_2d*cost_2d/ro    
    # Field is ordered from theta=pi (south pole) to theta=0 (north pole)
    # Average the Northern and Southern hemispheres together (must flip 
    # the Southern hemisphere with respect to latitude, then add or subtract it
    if sym=='even':
        field = 0.5*(field[it_half:, :] +\
                np.flip(field[:it_half, :], axis=0))
    elif sym=='odd':
        field = 0.5*(field[it_half:, :] -\
                np.flip(field[:it_half, :], axis=0))

    # Get the exponent to use for scientific notation
    # and normalize the field by 10**exp
    if not logscale:
        maxabs = max(np.abs(minmax[0]), np.abs(minmax[1]))
        exp = get_exp(maxabs)
        divisor = 10.**exp
        
        # Normalize field by divisor
        field /= divisor
        minmax = minmax[0]/divisor, minmax[1]/divisor
    
    # Create a default set of figure axes if they weren't already
    # specified by user
    if fig is None or ax is None:
        fig, ax = default_axes_1by1()
        showplot = True # probably in this case the user just
        # ran plot_azav from the command line wanting to view
        # view the plot
    
    # Specify linewidths to be used in the meridional plane, one for the 
    # boundary (lw) and one for the contours (contour_lw)
    lw = 1
    contour_lw = .2
    
    if (plotfield):
        plt.sca(ax)
#        levs = np.linspace(minmax[0], minmax[1], 100)
#        im = ax.contourf(xx, zz, field, cmap='RdYlBu_r',\
#        levels=levs)
#        plt.sca(ax)
        if logscale:
            plt.pcolormesh(xx, zz, field, cmap='Greys',\
                    norm=colors.LogNorm(vmin=minmax[0], vmax=minmax[1]))
        else:
            if posdef:
                cmap = 'plasma'
            else:
                cmap = 'RdYlBu_r'
            plt.pcolormesh(xx, zz, field, vmin=minmax[0], vmax=minmax[1],\
                    cmap=cmap)
        
        # Get the position of the axes on the figure
        ax_left, ax_right, ax_bottom, ax_top = axis_range(ax)
        ax_width = ax_right - ax_left
        ax_height = ax_top - ax_bottom
        ax_aspect = ax_height/ax_width
       
        # Set the colorbar axis to be in the "cavity" of the meridional plane
        # The colorbar height is set by making sure it "fits" in the cavity
        chi = ri/ro # aspect ratio of the shell
        cavity_height = ax_height*chi
        cbax_aspect = 10
        cbax_height = 0.7*cavity_height
        cbax_width = cbax_height/cbax_aspect / ax_aspect
        
        cbax_left = ax_left + 0.1*ax_width
        cbax_bottom = ax_bottom + 0.1*ax_height
        
        cbaxes = fig.add_axes([cbax_left, cbax_bottom,\
                       cbax_width, cbax_height])
        cbar = plt.colorbar(cax=cbaxes)

        cbaxes.tick_params(labelsize=fontsize)
        cbar.ax.tick_params(labelsize=fontsize)   #font size for the ticks

        if not logscale:
            if posdef:
                mid = (minmax[0] + minmax[1])/2.
                ticks = np.array([minmax[0], mid, minmax[1]])
            else:
                ticks = np.array([minmax[0], 0., minmax[1]])
            ticklabels = []
            for i in range(len(ticks)):
                ticklabels.append(str(round(ticks[i], 1)))
            ticks = np.array(ticks)
            cbar.set_ticks(ticks)
            cbar.set_ticklabels(ticklabels)
            cbar_label = (r'$\times10^{%i}\ $' %exp) + units
        else:
            cbar_label = units
        # Put the units and exponent to left of colorbar
        fig.text(cbax_left - 0.3*cbax_width, cbax_bottom + cbax_height/2.,\
                cbar_label, ha='right', va='center', rotation=90,\
                fontsize=fontsize)

    # Plot the boundary of the meridional plane
    plt.sca(ax)
    # outer boundary
    plt.plot(rr[0]/ro*sint, rr[0]/ro*cost, 'k', linewidth=lw)
    # inner boundary
    plt.plot(rr[-1]/ro*sint, rr[-1]/ro*cost, 'k', linewidth=lw) 
    # polar "edge"
    plt.plot([0,0], [rr[-1]/ro, rr[0]/ro], 'k', linewidth=lw) 
    # equatorial "edge"
    plt.plot([rr[-1]/ro, rr[0]/ro], [0,0], 'k', linewidth=lw) 

    # Plot latitude lines, if desired
    if plotlatlines: 
        lats = (np.pi/2. - np.arccos(cost))*180./np.pi
        lats_to_plot = np.arange(0., 90., 15.) 
        for lat in lats_to_plot:
            it = np.argmin(np.abs(lats - lat))
            xx_in, zz_in = rr[-1]/ro*sint[it], rr[-1]/ro*cost[it]
            xx_out, zz_out = rr[0]/ro*sint[it], rr[0]/ro*cost[it]
            plt.sca(ax)
            plt.plot([xx_in, xx_out], [zz_in, zz_out], 'k',\
                    linewidth=contour_lw)

    # Set axis ranges to be just outside the boundary lines
    lilbit = 0.01
    ax.set_xlim((-lilbit, 1 + lilbit))
    ax.set_ylim((-lilbit, 1 + lilbit))
    ax.axis('off') 

    # Plot contours in the meridional plane, if desired
    if plotcontours:
        # Determine the contour levels
        if levels is None:
            if logscale:
                levels = np.logspace(minexp, maxexp, nlevs)
            else:
                levels = np.linspace(minmax[0], minmax[1], nlevs)
        else: # the caller specified specific contour levels to plot!
            levels = np.array(levels)
        # Determine how to color the contours
        if logscale:
            contour_color = 'r'
        elif posdef:
            contour_color = 'w'
        else:
            contour_color = 'k'
        # plot the contours
        plt.sca(ax)
        plt.contour(xx, zz, field, colors=contour_color, levels=levels,\
                linewidths=contour_lw)
    if showplot:
        plt.show()

def plot_quiver(vr, vt, rr, cost, fig=None, ax=None, minmax=None,\
        plotlatlines=False, rvals=[], fontsize=10,\
        showplot=False, scale=None, plot_poles=False, nsample_t=20,\
        nsample_r=10, scale_by_mag=True, plotboundary=True):

    ''' Takes (or creates) set of axes with physical aspect ratio 1x2
    and adds a vector (quiver) plot of [vx, vy] in the meridional plane to 
    the axes.'''

    # First things first, make sure Python does not modify any of the 
    # arrays it was passed (shouldn't fucking have to do this)
    vmag = np.sqrt(np.std(vr**2 + vt**2))
    vr = np.copy(vr)/vmag
    vt = np.copy(vt)/vmag
    rr = np.copy(rr)
    cost = np.copy(cost)
    sint = np.sqrt(1. - cost**2.)

    # Derivative grid info
    ri, ro = np.min(rr), np.max(rr)
    nr = len(rr)
    nt = len(cost)

    # Create a default set of figure axes if they weren't already
    # specified by user
    if fig is None or ax is None:
        fig, ax = default_axes_1by2()
        showplot = True # probably in this case the user just
        # ran plot_azav from the command line wanting to view
        # view the plot
   
    # Get the position of the axes on the figure
    pos = ax.get_position().get_points()
    ax_left, ax_bottom = pos[0]
    ax_right, ax_top = pos[1]
    
    # Calculate the grid on which to plot
    rr2d = rr.reshape((1, nr))
    tt = np.arccos(cost)
    cost2d = cost.reshape((nt, 1))
    sint2d = sint.reshape((nt, 1))
    xx = rr2d*sint2d/ro
    zz = rr2d*cost2d/ro

    # Compute the vector field in cylindrical coordinates
    vx = sint2d*vr + cost2d*vt
    vz = cost2d*vr - sint2d*vt

    if not plot_poles:
        it_75_south = np.argmin(np.abs(tt - 11*np.pi/12))
        it_75_north = np.argmin(np.abs(tt - np.pi/12))
        xx = xx[it_75_south:it_75_north]
        zz = zz[it_75_south:it_75_north]
        vx = vx[it_75_south:it_75_north]
        vz = vz[it_75_south:it_75_north]
    if not scale_by_mag:
        v_mag = np.sqrt(vx**2 + vz**2)
        vx /= (20*v_mag)
        vz /= (20*v_mag)

    # Specify linewidths to be used in the meridional plane to plot the 
    # boundaries
    lw = 1.
    contour_lw = 0.2
   
    # Make quiver plot
    plt.sca(ax)
    # Scale down the arrays
    nskip_t, nskip_r = nt//nsample_t, nr//nsample_r
    xx = xx[::nskip_t, ::nskip_r]
    zz = zz[::nskip_t, ::nskip_r]
    vx = vx[::nskip_t, ::nskip_r]
    vz = vz[::nskip_t, ::nskip_r]

    plt.quiver(xx, zz, vx, vz, scale=scale)
   
    # Plot the boundary of the meridional plane
    if plotboundary:
        plt.plot(rr[0]/ro*sint, rr[0]/ro*cost, 'k', linewidth=lw)
        plt.plot(rr[-1]/ro*sint, rr[-1]/ro*cost, 'k', linewidth=lw)
        plt.plot([0.,0.], [rr[-1]/ro, rr[0]/ro], 'k', linewidth=lw)
        plt.plot([0.,0.], [-rr[-1]/ro, -rr[0]/ro], 'k', linewidth=lw)

    # Plot latitude lines, if desired
    if plotlatlines: 
        lats = (np.pi/2. - np.arccos(cost))*180./np.pi
        lats_to_plot = np.arange(-75., 90., 15.) 
        for lat in lats_to_plot:
            it = np.argmin(np.abs(lats - lat))
            xx_in, zz_in = rr[-1]/ro*sint[it], rr[-1]/ro*cost[it]
            xx_out, zz_out = rr[0]/ro*sint[it], rr[0]/ro*cost[it]
            plt.sca(ax)
            plt.plot([xx_in, xx_out], [zz_in, zz_out], 'k',\
                    linewidth=contour_lw)

    # Plot various radii, if desired
    # rvals must be given normalized to outer rr
    for rval in rvals: 
        plt.sca(ax)
        rval /= ro
        plt.plot(rval*sint, rval*cost, 'k--', linewidth=.7)

    # Set ax ranges to be just outside the boundary lines
    lilbit = 0.01
    ax.set_xlim((-lilbit, 1 + lilbit))
    ax.set_ylim((-1 - lilbit, 1 + lilbit))
    ax.axis('off') 

    if showplot:
        plt.show()

def streamfunction(vr,vt,r,cost,order=0):
    """------------------------------------------------------------
    This routine takes as input a divergenceless axymmetric 
    vector field in spherical coordinates and computes from 
    it a streamfunction (a.k.a. a flux flunction).  The grid
    is decribed by r and cost and can be non-uniform.
   ------------------------------------------------------------
    INPUTS:
   
    Vr, Vtheta = the 2-d vector velocity (or magnetic) field.
                 Dimensions are (N_Theta,N_R)
    r,cost     = the rr and cos(colatitude) of the grid.
                 r is assumed to vary from rmax to rmin and 
                 cost from  1 to -1 (i.e. 90 degrees
                 to -90 degrees in latitude).
                 Dimensions are r(N_R), cost(N_Theta)
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

def plot_azav_grid(terms, rr, cost, maintitle=None, titles=None, fig_width_inches=None, sub_width_inches=None, ncol=6, cmap='RdYlBu_r', units='', minmax=None, posdef=False, logscale=False, symlog=False, plotcontours=True, plotfield=True, nlevs=10, levels=None, plotlatlines=False, rvals=[], fontsize=default_labelsize, cbar_aspect=1.0/20.0, showplot=False, plot_cbar=True, lw_scaling=1.0, cbar_scaling=1.0, linthresh=None, linscale=None, plotboundary=True, rbcz=None, minmaxrz=None, linthreshrz=None, linscalerz=None, nosci=False, cbar_prec=1, latav=False, tw=None, totsig=None):

    # may need to make units an array
    nplots = len(terms)
    if isinstance(units, str): # units is just one label for everybody
        units = np.array([units]*nplots)

    # possibly sum some terms, based on totsig
    if not totsig is None:
        if totsig == 'sumrow':
            ncol_loc = ncol
            totsig = np.ones(ncol)
        else:
            ncol_loc = len(totsig)
        nrow_loc = len(terms)//ncol_loc

        iterm = 0
        titles = titles.tolist()
        units = units.tolist()
        for irow in range(nrow_loc):
            tot_term = np.zeros_like(terms[0])
            for icol in range(ncol_loc):
                tot_term += terms[iterm]*totsig[icol]
                iterm += 1

            # insert the tot_term at the correct place
            terms.insert(iterm, tot_term)
            titles.insert(iterm, 'tot')
            units.insert(iterm, units[iterm - 2])
            iterm += 1
            nplots += 1
        # need to made ncol 1 bigger to include the tot term
        ncol += 1

    # set up figure + axes
    sub_margin_bottom_inches = 0.75*(2.0 - (rbcz is None)) 
    if fig_width_inches is None and sub_width_inches is None:
        sub_width_inches = 2.0
    sub_aspect = 2.0

    # Generate the figure of the correct dimensions
    fig, axs, fpar = make_figure(nplots, sub_width_inches, sub_aspect, ncol=ncol, sub_margin_bottom_inches=sub_margin_bottom_inches)

    # possibly latitudinal average figure as well
    if latav:
        if rbcz is None:
            sub_margin_right_inches = None
        else:
            sub_margin_right_inches = default_margin_xlabel
        av_fig, av_axs, av_fpar = make_figure(nplots, ncol=ncol, sub_margin_top_inches=default_margin_label, sub_margin_left_inches=default_margin_ylabel, sub_margin_right_inches=sub_margin_right_inches, sub_margin_bottom_inches=default_margin_xlabel)
        xlabel = r'$r/R_\odot$'
        nt = len(cost)
        if tw is None: # just average everything unweighted
            tw = 1.0/nt + np.zeros(nt)
        tw_2d = tw.reshape((nt, 1))

    # plot all the terms
    for iplot in range(nplots):
        icol = iplot%fpar['ncol']
        irow = iplot//fpar['ncol']
        ax = axs[irow, icol]
        plot_azav(terms[iplot], rr, cost, fig=fig, ax=ax, cmap=cmap, units=units[iplot], minmax=minmax, posdef=posdef, logscale=logscale, symlog=symlog, plotcontours=plotcontours, plotfield=plotfield, nlevs=nlevs, levels=levels, plotlatlines=plotlatlines, rvals=rvals, fontsize=fontsize, cbar_aspect=cbar_aspect, showplot=showplot, plot_cbar=plot_cbar, lw_scaling=lw_scaling, cbar_scaling=cbar_scaling, linthresh=linthresh, linscale=linscale, plotboundary=plotboundary, rbcz=rbcz, minmaxrz=minmaxrz, linthreshrz=linthreshrz, linscalerz=linscalerz, nosci=nosci, cbar_prec=cbar_prec)

        if not titles is None:
            title_loc = '(' + letters[iplot] + ') ' + titles[iplot]
        else:
            title_loc = '(' + letters[iplot] + ')'
        ax.set_title(title_loc, loc='left', va='bottom', **csfont, fontsize=fontsize)

        # possibly plot the lat. average
        if latav:
            av_term = np.sum(terms[iplot]*tw_2d, axis=0)
            av_ax = av_axs[irow, icol]
            lineplot(rr/rsun, av_term, ax=av_ax, xlabel=xlabel, ylabel=units[iplot], title=titles[iplot], xcut=rbcz, xvals=rvals, yminmax=minmax, yminmax2=minmaxrz)

    # Put the main title in upper left
    fig.text(fpar['margin_left'] + fpar['sub_margin_left'], 1.0 - fpar['margin_top'], maintitle, ha='left', va='bottom', fontsize=fontsize, **csfont)

    if latav:
        av_fig.text(av_fpar['margin_left'] + av_fpar['sub_margin_left'], 1.0 - av_fpar['margin_top'], maintitle + ' (lat. avg.)', ha='left', va='bottom', fontsize=fontsize, **csfont)

    if latav:
        return fig, av_fig
    else:
        return fig
