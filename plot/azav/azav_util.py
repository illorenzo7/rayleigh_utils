# This file contains utilities useful for plotting the data from azimuthal
# averages,
# files in the directory AZ_Avgs.
# Written by Nick Featherstone
# First modified by Loren Matilsky, 03/09/2018

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
plt.rcParams['contour.negative_linestyle'] = 'solid'
import sys, os
sys.path.append(os.environ['rapl'])
sys.path.append(os.environ['raco'])
from common import *
from plotcommon import *

# default azav fig dimensions
azav_fig_dimensions = dict({'sub_width_inches': 2, 'sub_aspect': 2, 'sub_margin_left_inches': default_margin, 'sub_margin_top_inches': 1/4, 'sub_margin_bottom_inches': 1/2, 'margin_top_inches': 1})

# plot_azav needs my_contourf args, then some
plot_azav_kwargs_default = dict({'rbcz': None, 'minmaxrz': None, 'cmaprz': None, 'rvals': np.array([]), 'plotlatlines': True, 'latvals': np.arange(-60., 90., 30.), 'plotboundary': True,\
        'linestyles1': np.array(['-']), 'linewidths1': np.array([default_lw]), 'linecolors1': np.array(['k']),\
       'linestyles2': np.array(['-']), 'linewidths2': np.array([default_lw]), 'linecolors2': np.array(['k'])})

# add in my_contourf stuff
my_contourf_kwargs_default['cbar_aspect'] = 1/8
plot_azav_kwargs_default.update(my_contourf_kwargs_default)

def plot_azav(field, rr, cost, fig, ax,  **kwargs):
    find_bad_keys(plot_azav_kwargs_default, kwargs, 'plot_azav')
    kw = update_dict(plot_azav_kwargs_default, kwargs)
    kw_my_contourf = update_dict(my_contourf_kwargs_default, kwargs)

    # make copy of field
    field_full = np.copy(field)

    # grid info
    nt, nr = len(cost), len(rr)
    zeros = np.zeros((nt, nr))
    rr_2d = rr.reshape((1, nr)) + zeros
    rmax = np.max(rr_2d)
    cost_2d = cost.reshape((nt, 1)) + zeros
    sint_2d = np.sqrt(1.0 - cost_2d**2.0)
    tt_lat = 180.0/np.pi*(np.pi/2.0 - np.arccos(cost))

    # use these to plot lines/boundaries
    xx_full = rr_2d*sint_2d/rmax
    yy_full = rr_2d*cost_2d/rmax

    if kw.rbcz is None: # just plotting 1 domain
        xx = xx_full
        yy = yy_full
        field = field_full
    else: # plotting two domains
        irbcz = np.argmin(np.abs(rr/rsun - kw.rbcz))

        field = field[:, :irbcz+1]
        xx = (rr_2d*sint_2d)[:, :irbcz+1]/rmax
        yy = (rr_2d*cost_2d)[:, :irbcz+1]/rmax

        fieldrz = field_full[:, irbcz+1:]
        xxrz = (rr_2d*sint_2d)[:, irbcz+1:]/rmax
        yyrz = (rr_2d*cost_2d)[:, irbcz+1:]/rmax

    # plot the CZ field
    my_contourf(xx, yy, field, fig, ax, **kw_my_contourf)

    # possibly plot RZ field
    if not kw.rbcz is None: 
        # will need to change some contourf kwargs:
        kw_my_contourf.minmax = kw.minmaxrz
        kw_my_contourf.allticksoff = False # no need to turn off ticks twice
        if kw.cmaprz is None:
            if kw.posdef: 
                kw_my_contourf.cmap = 'cividis'
            else:
                kw_my_contourf.cmap = 'PuOr_r'    
        else:
            kw_my_contourf.cmap = kw.cmaprz
        kw_my_contourf.cbar_no = 2
        my_contourf(xxrz, yyrz, fieldrz, fig, ax, **kw_my_contourf)

    # potentially plot coordinate lines
    if not kw.plotlatlines:
        kw.latvals = np.array([])
    for ind in [1, 2]:
        if ind == 1:
            vals = make_array(kw.rvals, tolist=True)
            linecolors = kw.linecolors1
            linestyles = kw.linestyles1
            linewidths = kw.linewidths1
        if ind == 2:
            vals = make_array(kw.latvals, tolist=True)
            linecolors = kw.linecolors2
            linestyles = kw.linestyles2
            linewidths = kw.linewidths2

        # make lists from everything that needs to be
        linecolors = make_array(linecolors, tolist=True, length=len(vals))
        linestyles = make_array(linestyles, tolist=True, length=len(vals))
        linewidths = make_array(linewidths, tolist=True, length=len(vals))

        if kw.plotboundary:
            if ind == 1:
                to_add = [np.min(rr/rsun), np.max(rr/rsun)]
            if ind == 2:
                to_add = [np.min(tt_lat), np.max(tt_lat)]
            vals = [to_add[0]] + vals + [to_add[1]]
            linewidths = [default_lw] + linewidths + [default_lw] 
            linestyles = ['-'] + linestyles + ['-'] # make boundary lines solid
            linecolors = ['k'] + linecolors + ['k'] # make boundary lines black

        # make vals an array again
        for i in range(len(vals)):
            val = vals[i]
            if ind == 1:
                irval = np.argmin(np.abs(rr/rsun - val))
                xline, yline = xx_full[:, irval], yy_full[:, irval]
            if ind == 2:
                ilatval = np.argmin(np.abs(tt_lat - val))
                xline, yline = xx_full[ilatval, :], yy_full[ilatval, :]
            ax.plot(xline, yline, linewidth=linewidths[i], linestyle=linestyles[i], color=linecolors[i])

plot_azav_half_kwargs_default = dict(plot_azav_kwargs_default)
plot_azav_half_kwargs_default['sym'] = 'even'
plot_azav_half_kwargs_default['latvals'] = np.arange(0., 90., 30.)

def plot_azav_half(field, rr, cost, fig, ax,  **kwargs):
    '''Takes a figure with a subplot (axis) of aspect ratio 1x1 (or
    generates default axes if they are not provided) and adds
    a plot of the upper meridional plane to the axis, averaging the full
    plane with either even or odd symmetry
    '''
    find_bad_keys(plot_azav_half_kwargs_default, kwargs, 'plot_azav_half')
    kw = update_dict(plot_azav_half_kwargs_default, kwargs)
    kw_my_contourf = update_dict(my_contourf_kwargs_default, kwargs)

    # make copy of field
    field_full = np.copy(field)

    # may need to symmetrize array
    if not kw.sym is None:
        print('got here')

    # grid info
    nt, nr = len(cost), len(rr)
    zeros = np.zeros((nt, nr))
    rr_2d = rr.reshape((1, nr)) + zeros
    rmax = np.max(rr_2d)
    cost_2d = cost.reshape((nt, 1)) + zeros
    sint_2d = np.sqrt(1.0 - cost_2d**2.0)
    tt_lat = 180.0/np.pi*(np.pi/2.0 - np.arccos(cost))

    # use these to plot lines/boundaries
    xx_full = rr_2d*sint_2d/rmax
    yy_full = rr_2d*cost_2d/rmax

    if kw.rbcz is None: # just plotting 1 domain
        xx = xx_full
        yy = yy_full
        field = field_full
    else: # plotting two domains
        irbcz = np.argmin(np.abs(rr/rsun - kw.rbcz))

        field = field[:, :irbcz+1]
        xx = (rr_2d*sint_2d)[:, :irbcz+1]/rmax
        yy = (rr_2d*cost_2d)[:, :irbcz+1]/rmax

        fieldrz = field_full[:, irbcz+1:]
        xxrz = (rr_2d*sint_2d)[:, irbcz+1:]/rmax
        yyrz = (rr_2d*cost_2d)[:, irbcz+1:]/rmax

    # plot the CZ field
    my_contourf(xx, yy, field, fig, ax, **kw_my_contourf)

    # possibly plot RZ field
    if not kw.rbcz is None: 
        # will need to change some contourf kwargs:
        kw_my_contourf.minmax = kw.minmaxrz
        kw_my_contourf.allticksoff = False # no need to turn off ticks twice
        if kw.cmaprz is None:
            if kw.posdef: 
                kw_my_contourf.cmap = 'cividis'
            else:
                kw_my_contourf.cmap = 'PuOr_r'    
        else:
            kw_my_contourf.cmap = kw.cmaprz
        kw_my_contourf.cbar_no = 2
        my_contourf(xxrz, yyrz, fieldrz, fig, ax, **kw_my_contourf)

    # potentially plot coordinate lines
    if not kw.plotlatlines:
        kw.latvals = np.array([])
    for ind in [1, 2]:
        if ind == 1:
            vals = make_array(kw.rvals, tolist=True)
            linecolors = kw.linecolors1
            linestyles = kw.linestyles1
            linewidths = kw.linewidths1
        if ind == 2:
            vals = make_array(kw.latvals, tolist=True)
            linecolors = kw.linecolors2
            linestyles = kw.linestyles2
            linewidths = kw.linewidths2

        # make lists from everything that needs to be
        linecolors = make_array(linecolors, tolist=True, length=len(vals))
        linestyles = make_array(linestyles, tolist=True, length=len(vals))
        linewidths = make_array(linewidths, tolist=True, length=len(vals))

        if kw.plotboundary:
            if ind == 1:
                to_add = [np.min(rr/rsun), np.max(rr/rsun)]
            if ind == 2:
                to_add = [np.min(tt_lat), np.max(tt_lat)]
            vals = [to_add[0]] + vals + [to_add[1]]
            linewidths = [default_lw] + linewidths + [default_lw] 
            linestyles = ['-'] + linestyles + ['-'] # make boundary lines solid
            linecolors = ['k'] + linecolors + ['k'] # make boundary lines black

        # make vals an array again
        for i in range(len(vals)):
            val = vals[i]
            if ind == 1:
                irval = np.argmin(np.abs(rr/rsun - val))
                xline, yline = xx_full[:, irval], yy_full[:, irval]
            if ind == 2:
                ilatval = np.argmin(np.abs(tt_lat - val))
                xline, yline = xx_full[ilatval, :], yy_full[ilatval, :]
            ax.plot(xline, yline, linewidth=linewidths[i], linestyle=linestyles[i], color=linecolors[i])

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

plot_azav_grid_kwargs_default = dict({'maintitle': None, 'titles': None, 'shav': False, 'tw': None, 'totsig': None})
plot_azav_grid_kwargs_default.update(plot_azav_kwargs_default)
make_figure_kwargs_default.update(azav_fig_dimensions)
plot_azav_grid_kwargs_default.update(make_figure_kwargs_default)
# need a nother make_figure_kwargs for the shav plot (possibly)
make_figure_kwargs_default_shav = dict({**make_figure_kwargs_default})
make_figure_kwargs_default_shav.update(lineplot_fig_dimensions)
for key, val in make_figure_kwargs_default_shav.items():
    plot_azav_grid_kwargs_default[key + '_shav'] = val

def plot_azav_grid(terms, rr, cost, **kwargs):
    find_bad_keys(plot_azav_grid_kwargs_default, kwargs, 'plot_azav')

    kw = update_dict(plot_azav_grid_kwargs_default, kwargs)
    kw_plot_azav = update_dict(plot_azav_kwargs_default, kwargs)
    kw_make_figure = update_dict(make_figure_kwargs_default, kwargs)
    kw_make_figure_shav = dict({})
    for key in make_figure_kwargs_default_shav:
        kw_make_figure_shav[key] = kw[key + '_shav']
    kw_make_figure_shav = dotdict(kw_make_figure_shav)

    # possibly sum some terms, based on totsig
    nplots = len(terms)
    if not kw.totsig is None:
        if kw.ncol is None:
            kw.ncol = nplots
        nrow = len(terms)//kw.ncol

        # totsig works by summing over rows
        # some rows may get summed differently (but try to avoid this)
        if np.isscalar(kw.totsig):
            if kw.totsig == 'sumrow':
                ncol_loc = kw.ncol
                kw.totsig = np.ones(len(terms))
        if len(kw.totsig) < len(terms): # must repeat according to to rows
            # first make sure totsig is the length of a row
            if len(kw.totsig) == kw.ncol:
                kw.totsig = kw.totsig.tolist()
                kw.totsig *= nrow
                kw.totsig = np.array(kw.totsig)
            else:
                print ('ERROR: (len(totsig) = %i) < (len(terms) = %i)' %(len(kw.totsig), len(terms)))
                print ('but (len(totsig) = %i) != (ncol = %i)' %(len(kw.totsig), kw.ncol))
                print ('exiting')
                sys.exit()

        count = 0
        iterm = 0
        kw.titles = kw.titles.tolist()
        for irow in range(nrow):
            tot_term = np.zeros_like(terms[0])
            for icol in range(kw.ncol):
                tot_term += terms[iterm]*kw.totsig[count]
                iterm += 1
                count += 1

            # insert the tot_term at the correct place
            terms.insert(iterm, tot_term)
            kw.titles.insert(iterm, 'tot')
            iterm += 1
            nplots += 1

        # need to made ncol 1 bigger to include the tot term
        kw.ncol += 1

    # make plot
    kw_make_figure.nplots = nplots
    kw_make_figure.ncol = kw.ncol
    if not kw.rbcz is None:
        kw_make_figure.sub_margin_bottom_inches *= 2
    fig, axs, fpar = make_figure(**kw_make_figure)

    # possibly latitudinal average figure as well
    if kw.shav:
        kw_make_figure_shav.nplots = nplots
        kw_make_figure_shav.ncol = kw.ncol
        if not kw.rbcz is None:
            kw_make_figure_shav.sub_margin_bottom_inches *= 2

        if kw.rbcz is None:
            kw_make_figure_shav.sub_margin_right_inches = default_margin
        else:
            kw_make_figure_shav.sub_margin_right_inches = default_margin_ylabel

        av_fig, av_axs, av_fpar = make_figure(**kw_make_figure_shav)
        xlabel = r'$r/R_\odot$' 
        nt = len(cost)
        if kw.tw is None: # just average everything unweighted
            kw.tw = 1.0/nt + np.zeros(nt)
        tw_2d = kw.tw.reshape((nt, 1))

    # plot all the terms
    for iplot in range(nplots):
        icol = iplot%fpar['ncol']
        irow = iplot//fpar['ncol']
        ax = axs[irow, icol]
        plot_azav(terms[iplot], rr, cost, fig, ax, **kw_plot_azav)

        if not kw.titles is None:
            title_loc = '(' + letters[iplot] + ') ' + kw.titles[iplot]
        else:
            title_loc = '(' + letters[iplot] + ')'
        ax.set_title(title_loc, loc='left', va='bottom', fontsize=default_titlesize)

        # possibly plot the lat. average
        if kw.shav:
            av_term = np.sum(terms[iplot]*tw_2d, axis=0)
            av_ax = av_axs[irow, icol]
            lineplot(rr/rsun, [av_term], av_ax, xlabel=xlabel, title=title_loc, xcut=kw.rbcz, xvals=kw.rvals, minmax=kw.minmax, minmax2=kw.minmaxrz, plotleg=False)

    # Put the main title in upper left
    fig.text(fpar['margin_left'] + fpar['sub_margin_left'], 1.0 - fpar['margin_top'], kw.maintitle, ha='left', va='bottom', fontsize=default_titlesize)

    if kw.shav:
        av_fig.text(av_fpar['margin_left'] + av_fpar['sub_margin_left'], 1.0 - av_fpar['margin_top'], kw.maintitle + ' (Shell Averaged)', ha='left', va='bottom', fontsize=default_titlesize)

    if kw.shav:
        return fig, av_fig
    else:
        return fig
