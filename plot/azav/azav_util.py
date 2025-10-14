# Author: Loren Matilsky, adapted from Nick Featherstone
# Created: 12/19/2022
#
# Description: Module for routines useful for plotting AZ_Avgs data

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
plt.rcParams['contour.negative_linestyle'] = 'solid'
import sys, os
sys.path.append(os.environ['rapl'])
sys.path.append(os.environ['raco'])
from common import *
from plotcommon import *
from grid_util import compute_theta_grid

# default azav fig dimensions
azav_fig_dimensions = dict({'sub_aspect': 2, 'sub_width_inches': 2, 'sub_margin_left_inches': 3/4,  'sub_margin_right_inches': 1/2, 'sub_margin_top_inches': 1/4, 'sub_margin_bottom_inches': 1, 'margin_top_inches': 1})

# plot_azav needs my_contourf args, then some
kw_plot_azav_default = dict(
        {'rcut': None, 'order': 'down', 'minmax2': None, 'cmap2': None, 'rvals': None, 'plotlatlines': True, 'latvals': np.arange(-60., 90., 30.), 'plotboundary': True,
        'linestyles1': np.array(['--']), 'linewidth': default_lw, 'linecolors1': np.array(['k']),
       'linestyles2': np.array(['--']), 'linecolors2': np.array(['k']),
       'halfplane': False, 'sym': False, 'antisym': False, 'fontsize': default_labelsize, 'plotaxis': True,\
        'modrms': False
       })

# add in my_contourf stuff
kw_plot_azav_default.update(kw_my_contourf_default)

def plot_azav(field, rr, cost, fig, ax,  **kw_in):
    find_bad_keys(kw_plot_azav_default, kw_in, 'plot_azav')
    kw = update_dict(kw_plot_azav_default, kw_in)
    kw_my_contourf = update_dict(kw_my_contourf_default, kw_in)

    # make copies of field and costheta
    # (effectively stop Python from passing arrays by reference)
    field = np.copy(field)
    cost = np.copy(cost)

    # if "half plane" only plot northern hemisphere
    if kw.halfplane:
        iteq = len(cost)//2
        if kw.sym:
            field = 0.5*(field[iteq:, :] + field[:iteq, :][::-1, :])
        elif kw.antisym:
            field = 0.5*(field[iteq:, :] - field[:iteq, :][::-1, :])
        else:
            field = field[iteq:, :]

        cost = cost[iteq:]

    if kw.plotaxis: # make cbar offset bigger
        if kw_my_contourf.cbar_offset is None:
            kw_my_contourf.cbar_offset = 1/2

    # grid info
    nt, nr = len(cost), len(rr)
    zeros = np.zeros((nt, nr))
    rr_2d = rr.reshape((1, nr)) + zeros
    rmax = np.max(rr_2d)
    cost_2d = cost.reshape((nt, 1)) + zeros
    sint_2d = np.sqrt(1.0 - cost_2d**2.0)
    tt_lat = 180.0/np.pi*(np.pi/2.0 - np.arccos(cost))

    # possibly remove the rms field
    if kw.modrms:
        tt, tw = compute_theta_grid(nt)
        field_rms = np.sqrt(np.sum(field**2*tw.reshape((nt, 1)), axis=0))
        field /= field_rms.reshape((1, nr))

    # use these to plot lines/boundaries
    xx_full = rr_2d*sint_2d/rmax
    yy_full = rr_2d*cost_2d/rmax

    if kw.rcut is None: # just plotting 1 domain
        xx = xx_full
        yy = yy_full
    else: # plotting 2 domains
        ircut = np.argmin(np.abs(rr - kw.rcut))

        field_full = np.copy(field)

        if kw.order == 'down':
            field = field_full[:, :ircut+1]
            xx = (rr_2d*sint_2d)[:, :ircut+1]/rmax
            yy = (rr_2d*cost_2d)[:, :ircut+1]/rmax

            field2 = field_full[:, ircut+1:]
            xx2 = (rr_2d*sint_2d)[:, ircut+1:]/rmax
            yy2 = (rr_2d*cost_2d)[:, ircut+1:]/rmax

        else:
            field = field_full[:, ircut+1:]
            xx = (rr_2d*sint_2d)[:, ircut+1:]/rmax
            yy = (rr_2d*cost_2d)[:, ircut+1:]/rmax

            field2 = field_full[:, :ircut+1]
            xx2 = (rr_2d*sint_2d)[:, :ircut+1]/rmax
            yy2 = (rr_2d*cost_2d)[:, :ircut+1]/rmax


    # plot the first (upper) field
    kw_my_contourf.cbar_no = 1
    my_contourf(xx, yy, field, fig, ax, **kw_my_contourf)

    # possibly plot a second (lower) field
    if not kw.rcut is None: 
        # will need to change some contourf kw:
        kw_my_contourf.minmax = kw.minmax2
        kw_my_contourf.allticksoff = False # no need to turn off ticks twice
        if kw.cmap2 is None:
            if kw.posdef: 
                kw_my_contourf.cmap = 'cividis'
            else:
                kw_my_contourf.cmap = 'PuOr_r'    
        else:
            kw_my_contourf.cmap = kw.cmap2
        kw_my_contourf.cbar_no = 2
        my_contourf(xx2, yy2, field2, fig, ax, **kw_my_contourf)

    # potentially plot coordinate lines
    if not kw.plotlatlines:
        kw.latvals = np.array([])
    if kw.rvals is None:
        kw.rvals = np.array([])

    for ind in [1, 2]:
        if ind == 1:
            vals = make_array(kw.rvals, tolist=True)
            linecolors = kw.linecolors1
            linestyles = kw.linestyles1
        if ind == 2:
            if kw.halfplane:
                # exclude latvals close to boundary, which now includes equator
                tol = 5. # exclude latvals within 5 degrees
                latvals_new = []
                for latval in kw.latvals:
                    if latval > 5.:
                        latvals_new.append(latval)
                kw.latvals = np.array(latvals_new)

            vals = make_array(kw.latvals, tolist=True)
            linecolors = kw.linecolors2
            linestyles = kw.linestyles2

        # make lists from everything that needs to be
        linecolors = make_array(linecolors, tolist=True, length=len(vals))
        linestyles = make_array(linestyles, tolist=True, length=len(vals))

        if kw.plotboundary:
            if ind == 1:
                to_add = [np.min(rr), np.max(rr)]
            if ind == 2:
                to_add = [np.min(tt_lat), np.max(tt_lat)]
            vals = [to_add[0]] + vals + [to_add[1]]
            linestyles = ['-'] + linestyles + ['-'] # make boundary lines solid
            linecolors = ['k'] + linecolors + ['k'] # make boundary lines black

        # make vals an array again
        for i in range(len(vals)):
            val = vals[i]
            if ind == 1:
                irval = np.argmin(np.abs(rr - val))
                xline, yline = xx_full[:, irval], yy_full[:, irval]
            if ind == 2:
                ilatval = np.argmin(np.abs(tt_lat - val))
                xline, yline = xx_full[ilatval, :], yy_full[ilatval, :]
            ax.plot(xline, yline, linewidth=kw.linewidth, linestyle=linestyles[i], color=linecolors[i])
        
    if kw.plotaxis: # plot x and z coordinate axes
        ax.axis('on')
        plt.sca(ax)
        plt.minorticks_on()
        plt.tick_params(top=True, right=True, direction='in', which='both',\
                labelsize=kw.fontsize)
        plt.xlabel(r'$x/r_{\rm{out}}$', fontsize=kw.fontsize)
        plt.ylabel(r'$z/r_{\rm{out}}$', fontsize=kw.fontsize)
    else:
        ax.axis('off')

    # reset the axis limits (contourf has already done its thing;
    # but for this one we either want 1 x 1 (half-plane) or 1 x 2
    lilbit = 0.01
    ax.set_xlim(-lilbit, 1.0 + lilbit)
    if kw.halfplane:
        ax.set_ylim(-lilbit, 1.0 + lilbit)
    else:
        ax.set_ylim(-1.0 -lilbit, 1.0 + lilbit)

    if kw.modrms:
        rmin = np.min(rr)
        shell_depth = rmax - rmin
        height = (rr - rmin)/shell_depth
        ymin, ymax = ax.get_ylim()
        ymin -= 1.0
        ax.set_ylim(ymin, ymax)

        ax2 = ax.twinx()
        ax2.plot(height[1:-1], field_rms[1:-1])
        ax2.set_yscale('log')
        ymin = np.min(field_rms[1:-1])
        ymax = np.max(field_rms[1:-1])
        Dy = ymax/ymin
        ax2.set_ylim(ymin, ymax*Dy**2)
        ax.set_xlabel(r'$x/r_{\rm{out}}$' + ' (or height)', fontsize=kw.fontsize)

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
            psi2[0:n_t-1,i] = psi2[0:n_t-1,i+1] + dpsi_dr[0:n_t-1,i]*dr[i]
        for i in range(n_t-2, -1, -1):
            psi2[i,0:n_r-1] = psi2[i+1,0:n_r-1] + dpsi_dt[i,0:n_r-1]*dtheta[i]
        
        if (order < 0):
            return psi2
        else:
            psi=0.5*(psi+psi2)
            
    return psi

kw_plot_azav_grid_default = dict({'maintitle': None, 'titles': None, 'shav': False, 'sub': False, 'tw': None, 'totsig': None, 'minmaxs': None, 'iplots': None, 'groupname': None})
kw_plot_azav_grid_default.update(kw_plot_azav_default)
kw_make_figure_az = kw_make_figure_default.copy()
kw_make_figure_az.update(azav_fig_dimensions)
kw_plot_azav_grid_default.update(kw_make_figure_az)
# need another make_figure_kw for the shav plot (possibly)
kw_make_figure_shav = kw_make_figure_az.copy()
kw_make_figure_shav.update(lineplot_fig_dimensions)
for key, val in kw_make_figure_shav.items():
    kw_plot_azav_grid_default[key + '_shav'] = val

def plot_azav_grid(terms, rr, cost, **kw_in):
    find_bad_keys(kw_plot_azav_grid_default, kw_in, 'plot_azav')

    kw = update_dict(kw_plot_azav_grid_default, kw_in)
    kw_plot_azav = update_dict(kw_plot_azav_default, kw_in)
    kw_make_figure = update_dict(kw_make_figure_az, kw_in)
    kw_make_figure_shav = dict({})
    for key in kw_make_figure_shav:
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
        if not isinstance(kw.titles, list):
            kw.titles = kw.titles.tolist()
        for irow in range(nrow):
            tot_term = np.zeros_like(terms[0])
            for icol in range(kw.ncol):
                tot_term += terms[iterm]*kw.totsig[count]
                # if totsig != 0,
                # also update the term itself with the proper sign
                if kw.totsig[count] != 0:
                    terms[iterm] *= kw.totsig[count]

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
    if not kw.rcut is None:
        kw_make_figure.sub_margin_bottom_inches *= 2
    fig, axs, fpar = make_figure(**kw_make_figure)

    # may need these theta weights
    nt = len(cost)
    nr = len(rr)
    if kw.tw is None: # just average everything unweighted
        kw.tw = 1.0/nt + np.zeros(nt)
    tw_2d = kw.tw.reshape((nt, 1))

    # possibly latitudinal average figure as well
    if kw.shav:
        kw_make_figure_shav.nplots = nplots
        kw_make_figure_shav.ncol = kw.ncol
        if not kw.rcut is None:
            kw_make_figure_shav.sub_margin_bottom_inches *= 2

        if kw.rcut is None:
            kw_make_figure_shav.sub_margin_right_inches = default_margin
        else:
            kw_make_figure_shav.sub_margin_right_inches = default_margin_ylabel

        av_fig, av_axs, av_fpar = make_figure(**kw_make_figure_shav)
        xlabel = r'$r/R_\odot$' 

    # plot all the terms
    count = 0 # for minmaxs, if specified
    for iplot in range(nplots):
        icol = iplot%fpar['ncol']
        irow = iplot//fpar['ncol']
        ax = axs[irow, icol]

        # maybe need spherical avg
        if kw.shav or kw.sub:
            av_term = np.sum(terms[iplot]*tw_2d, axis=0)
            if kw.sub:
                terms[iplot] -= av_term.reshape((1, nr))

        # check if we need minmax from minmaxs
        if not kw.iplots is None:
            if iplot in make_array(kw.iplots): # these are the exceptions to default
                kw_plot_azav.minmax = kw.minmaxs[2*count:2*(count+1)]
                count += 1
            else: # this is default
                kw_plot_azav.minmax = kw.minmax

        # check if groupname is 'v':
        if kw.groupname == 'v' and kw.plotcontours: # still eliminate contours on v_r and v_theta
            if iplot in [0, 1]:
                kw_plot_azav.plotcontours = False
            else:
                kw_plot_azav.plotcontours = True

        plot_azav(terms[iplot], rr, cost, fig, ax, **kw_plot_azav)

        if not kw.titles is None:
            title_loc = '(' + letters[iplot] + ') ' + kw.titles[iplot]
        else:
            title_loc = '(' + letters[iplot] + ')'
        ax.set_title(title_loc, loc='left', va='bottom', fontsize=default_titlesize)

        # possibly plot the lat. average
        if kw.shav:
            av_ax = av_axs[irow, icol]
            lineplot(rr, [av_term], av_ax, xlabel=xlabel, title=title_loc, xcut=kw.rcut,  minmax=kw.minmax, minmax2=kw.minmax2, plotleg=False)

    # Put the main title in upper left
    # check if there is an rcut
    if not kw.rcut is None:
        kw.maintitle += ('\nrcut = %1.3e' %kw.rcut)
    fig.text(fpar['margin_left'] + fpar['sub_margin_left'], 1.0 - fpar['margin_top'], kw.maintitle, ha='left', va='bottom', fontsize=default_titlesize)

    if kw.shav:
        av_fig.text(av_fpar['margin_left'] + av_fpar['sub_margin_left'], 1.0 - av_fpar['margin_top'], kw.maintitle + ' (Shell Averaged)', ha='left', va='bottom', fontsize=default_titlesize)

    if kw.shav:
        return fig, av_fig
    else:
        return fig
