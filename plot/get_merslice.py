# Author: Loren Matilsky
# Created: 02/02/2020
#
# Extremely long and boring script to find fundamental fluid quantities
# or derivative quantities from a meridional slice. 
# Takes a meridional slice [mer] and
# varname in [vr, vt, vp, vl, vz, ...]
# returns the slice for the variable as an array of shape 
# (nphi, ntheta, nr) here nphi is the number of longitudes sampled by
# the meridional slice
import numpy as np
import sys, os
sys.path.append(os.environ['rapp'])
sys.path.append(os.environ['raco'])
from varprops import var_indices, var_indices_old
from common import *
from get_eq import get_eq
from rayleigh_diagnostics import GridInfo

def get_merslice(mer, varname, dirname=None, az=None, sh=None, old=False,\
        j=0):
    # Given a Meridional_Slices object (nphi, ntheta, nr, nq, niter), 
    # return the field (nphi, ntheta, nr) associated with [varname] 
    # take the jth time slice (default j=0).
    # For "primed" variables, the routine will use provided AZ_Avgs object 
    # [az], or else estimate the longitudinal average from [mer].
    # for "sphericallly primed" variables, the routine will use provided
    # Shell_Avgs object [sh], or else estimate spherical average from [mer].
    # If [az] or [sh] are not provided for _prime or _prime_sph vars,
    # will use integration weights (from GridInfo) to subtract off various 
    # averages 
    # reference-state parameters (get_eq) come from [dirname]
    # _prime refers to az-avg subtracted
    # _prime_sph refers to sph-avg subtracted

    # Get shape of grid
    nphi = mer.nphi
    nt = mer.ntheta
    nr = mer.nr

    # Get sine and cosine if needed
    if 'l' in varname or 'z' in varname:
        sint = (mer.sintheta).reshape((1, nt, 1))
        cost = (mer.costheta).reshape((1, nt, 1))        
   
    # get reference-state stuff if needed
    if not dirname is None:
        eq = get_eq(dirname)
        ref_rho = (eq.density).reshape((1, 1, nr))
        ref_T = (eq.temperature).reshape((1, 1, nr))
        ref_P = (eq.pressure).reshape((1, 1, nr))
  
    # Get raw values associated with jth time slice
    vals = mer.vals[:, :, :, :, j]

    # Get quantity indexes (new or old)---do I actually ever need the 
    # old ones again?
    qind_vr = var_indices['vr']
    qind_vt = var_indices['vt']
    qind_vp = var_indices['vp']
    
    qind_s = var_indices['s']
    qind_p = var_indices['p']
    
    qind_omr = var_indices['omr']
    qind_omt = var_indices['omt']
    qind_omp = var_indices['omp']
    
    qind_br = var_indices['br']
    qind_bt = var_indices['bt']
    qind_bp = var_indices['bp']
    
    if old:
        qind_vr = var_indices_old['vr']
        qind_vt = var_indices_old['vt']
        qind_vp = var_indices_old['vp']
        
        qind_omr = var_indices_old['omr']
        qind_omt = var_indices_old['omt']
        qind_omp = var_indices_old['omp']
        
        qind_s = var_indices_old['s']
        qind_p = var_indices_old['p']       
               
        qind_br = var_indices_old['br']
        qind_bt = var_indices_old['bt']
        qind_bp = var_indices_old['bp']
    
    # Spherical velocities
    if varname == 'vr':
        ind_vr = mer.lut[qind_vr]
        merslice = vals[:, :, :, ind_vr]/100. # measure velocities in m/s
    elif varname == 'vt':
        ind_vt = mer.lut[qind_vt]
        merslice = vals[:, :, :, ind_vt]/100.
    elif varname == 'vp':
        ind_vp = mer.lut[qind_vp]
        merslice = vals[:, :, :, ind_vp]/100.
    elif varname == 'vr_prime':
        ind_vr = mer.lut[qind_vr]
        vr_slice = vals[:, :, :, ind_vr]/100.
        if az is None:
            print("CAUTION: for %s, getting azav from merslice" %varname)
            vr_az = np.mean(vr_slice, axis=0)
        else:
            vr_az = az.vals[:, :, az.lut[qind_vr], j]/100.
        merslice = vr_slice - vr_az.reshape((1, nt, nr))
    elif varname == 'vt_prime':
        ind_vt = mer.lut[qind_vt]
        vt_slice = vals[:, :, :, ind_vt]/100.
        if az is None:
            print("CAUTION: for %s, getting azav from merslice" %varname)
            vt_az = np.mean(vt_slice, axis=0)
        else:
            vt_az = az.vals[:, :, az.lut[qind_vt], j]/100.
        merslice = vt_slice - vt_az
    elif varname == 'vp_prime':
        ind_vp = mer.lut[qind_vp]
        vp_slice = vals[:, :, :, ind_vp]/100.
        if az is None:
            print("CAUTION: for %s, getting azav from merslice" %varname)
            vp_az = np.mean(vp_slice, axis=0)
        else:
            vp_az = az.vals[:, :, az.lut[qind_vp], j]/100.
        merslice = vp_slice - vp_az
            
    # Cylindrical-coordinate velocities
    elif varname == 'vl':
        ind_vr, ind_vt = mer.lut[qind_vr], mer.lut[qind_vt]
        vr_slice = vals[:, :, :, ind_vr]/100.
        vt_slice = vals[:, :, :, ind_vt]/100.
        merslice = vr_slice*sint + vt_slice*cost       
    elif varname == 'vz':
        ind_vr, ind_vt = mer.lut[qind_vr], mer.lut[qind_vt]
        vr_slice = vals[:, :, :, ind_vr]/100.
        vt_slice = vals[:, :, :, ind_vt]/100.
        merslice = vr_slice*cost - vt_slice*sint   
    elif varname == 'vl_prime':
        ind_vr, ind_vt = mer.lut[qind_vr], mer.lut[qind_vt]        
        vr_slice = vals[:, :, :, ind_vr]/100.
        vt_slice = vals[:, :, :, ind_vt]/100.    
        vl_slice = vr_slice*sint + vt_slice*cost
        if az is None:
            print("CAUTION: for %s, getting azav from merslice" %varname)
            vl_az = np.mean(vl_slice, axis=0)
        else:
            vr_az = az.vals[:, :, az.lut[qind_vr], j]/100.
            vt_az = az.vals[:, :, az.lut[qind_vt], j]/100.
            vl_az = vr_az*sint + vt_az*cost
        merslice = vl_slice - vl_az
    elif varname == 'vz_prime':
        ind_vr, ind_vt = mer.lut[qind_vr], mer.lut[qind_vt]        
        vr_slice = vals[:, :, :, ind_vr]/100.
        vt_slice = vals[:, :, :, ind_vt]/100.    
        vz_slice = vr_slice*cost - vt_slice*sint  
        if az is None:
            print("CAUTION: for %s, getting azav from merslice" %varname)
            vz_az = np.mean(vz_slice, axis=0)
        else:
            vr_az = az.vals[:, :, az.lut[qind_vr], j]/100.
            vt_az = az.vals[:, :, az.lut[qind_vt], j]/100.
            vz_az = vr_az*sint - vt_az*cost
        merslice = vz_slice - vz_az
            
    # Spherical vorticities
    elif varname == 'omr':
        ind_omr = mer.lut[qind_omr]
        merslice = vals[:, :, :, ind_omr]
    elif varname == 'omt':
        ind_omt = mer.lut[qind_omt]
        merslice = vals[:, :, :, ind_omt]
    elif varname == 'omp':
        ind_omp = mer.lut[qind_omp]
        merslice = vals[:, :, :, ind_omp]
    elif varname == 'omr_prime':
        ind_omr = mer.lut[qind_omr]
        omr_slice = vals[:, :, :, ind_omr]
        if az is None:
            print("CAUTION: for %s, getting azav from merslice" %varname)
            omr_az = np.mean(omr_slice, axis=0)
        else:
            omr_az = az.vals[:, :, az.lut[qind_omr], j]
        merslice = omr_slice - omr_az.reshape((1, nt, nr))
    elif varname == 'omt_prime':
        ind_omt = mer.lut[qind_omt]
        omt_slice = vals[:, :, :, ind_omt]
        if az is None:
            print("CAUTION: for %s, getting azav from merslice" %varname)
            omt_az = np.mean(omt_slice, axis=0)
        else:
            omt_az = az.vals[:, :, az.lut[qind_omt], j]
        merslice = omt_slice - omt_az
    elif varname == 'omp_prime':
        ind_omp = mer.lut[qind_omp]
        omp_slice = vals[:, :, :, ind_omp]
        if az is None:
            print("CAUTION: for %s, getting azav from merslice" %varname)
            omp_az = np.mean(omp_slice, axis=0)
        else:
            omp_az = az.vals[:, :, az.lut[qind_omp], j]/100.
        merslice = omp_slice - omp_az
            
    # Cylindrical vorticities
    elif varname == 'oml':
        ind_omr, ind_omt = mer.lut[qind_omr], mer.lut[qind_omt]
        omr_slice = vals[:, :, :, ind_omr]
        omt_slice = vals[:, :, :, ind_omt]
        merslice = omr_slice*sint + omt_slice*cost       
    elif varname == 'omz':
        ind_omr, ind_omt = mer.lut[qind_omr], mer.lut[qind_omt]
        omr_slice = vals[:, :, :, ind_omr]
        omt_slice = vals[:, :, :, ind_omt]
        merslice = omr_slice*cost - omt_slice*sint
    elif varname == 'omz_prime':
        ind_omr, ind_omt = mer.lut[qind_omr], mer.lut[qind_omt]
        omr_slice = vals[:, :, :, ind_omr]
        omt_slice = vals[:, :, :, ind_omt]
        omz_slice = omr_slice*cost - omt_slice*sint
        if az is None:
            print("CAUTION: for %s, getting azav from merslice" %varname)
            omz_az = np.mean(omz_slice, axis=0)
        else:
            omr_az = az.vals[:, :, az.lut[qind_omr], j]
            omt_az = az.vals[:, :, az.lut[qind_omt], j]
            omz_az = omr_az*sint - omt_az*cost
        merslice = omz_slice - omz_az
            
    # Enstrophy (total and fluctuating)
    elif varname == 'omsq':
        ind_omr, ind_omt, ind_omp = mer.lut[qind_omr], mer.lut[qind_omt],\
                mer.lut[qind_omp]
        merslice = vals[:, :, :, ind_omr]**2 + vals[:, :, :, ind_omt]**2 +\
                vals[:, :, :, ind_omp]**2 
    elif varname == 'omsq_prime':
        ind_omr, ind_omt, ind_omp = mer.lut[qind_omr], mer.lut[qind_omt],\
                mer.lut[qind_omp]
        omr_slice = vals[:, :, :, ind_omr]
        omt_slice = vals[:, :, :, ind_omt]
        omp_slice = vals[:, :, :, ind_omp]
        if az is None:
            print("CAUTION: for %s, getting azav from merslice" %varname)
            omr_az = np.mean(omr_slice, axis=0)
            omt_az = np.mean(omt_slice, axis=0)
            omp_az = np.mean(omp_slice, axis=0)
        else:
            omr_az = az.vals[:, :, az.lut[qind_omr], j]
            omt_az = az.vals[:, :, az.lut[qind_omt], j]
            omp_az = az.vals[:, :, az.lut[qind_omp], j]
        omr_prime = omr_slice - omr_az
        omt_prime = omt_slice - omt_az
        omp_prime = omp_slice - omp_az
        merslice = omr_prime**2 + omt_prime**2 + omp_prime**2
        
    # Thermodynamic variables: deviations from reference state
    elif varname == 's':
        ind_s = mer.lut[qind_s]
        merslice = vals[:, :, :, ind_s]
    elif varname == 'p':
        ind_p = mer.lut[qind_p]
        merslice = vals[:, :, :, ind_p]
    
    # Thermodynamic variables: deviations from zonal mean
    elif varname == 's_prime':
        ind_s = mer.lut[qind_s]
        s_slice = vals[:, :, :, ind_s]
        if az is None:
            print("CAUTION: for %s, getting azav from merslice" %varname)
            s_az = np.mean(s_slice, axis=0)
        else:
            s_az = az.vals[:, :, az.lut[qind_s], j]
        merslice = s_slice - s_az
    elif varname == 'p_prime':
        ind_p = mer.lut[qind_p]
        p_slice = vals[:, :, :, ind_p]
        if az is None:
            print("CAUTION: for %s, getting azav from merslice" %varname)
            p_az = np.mean(p_slice, axis=0)
        else:
            p_az = az.vals[:, :, az.lut[qind_p], j]
        merslice = p_slice - p_az
            
    # Thermodynamic variables: deviations from spherical mean
    elif varname == 's_prime_sph':
        ind_s = mer.lut[qind_s]
        s_slice = vals[:, :, :, ind_s]
        if sh is None:
            print("CAUTION: for %s, getting azav from merslice" %varname)
            gi = GridInfo(dirname + '/grid_info')
            tw = gi.tweights.reshape((nt, 1))
            s_sh = np.mean(s_slice, axis=0) # first, average over lon.
            s_sh = np.sum(s_sh*tw, axis=0) # then over lat
        else:
            s_sh = sh.vals[:, 0, sh.lut[qind_s], j]
        merslice = s_slice - s_sh.reshape((1, 1, nr))
    elif varname == 'p_prime_sph':
        ind_p = mer.lut[qind_p]
        p_slice = vals[:, :, :, ind_p]
        if sh is None:
            print("CAUTION: for %s, getting azav from merslice" %varname)
            gi = GridInfo(dirname + '/grid_info')
            tw = gi.tweights.reshape((nt, 1))
            p_sh = np.mean(p_slice, axis=0) # first, average over lon.
            p_sh = np.sum(p_sh*tw, axis=0) # then over lat
        else:
            p_sh = sh.vals[:, 0, sh.lut[qind_p], j]
        merslice = p_slice - p_sh.reshape((1, 1, nr))

    # Spherical magnetic fields
    elif varname == 'br':
        ind_br = mer.lut[qind_br]
        merslice = vals[:, :, :, ind_br]
    elif varname == 'bt':
        ind_bt = mer.lut[qind_bt]
        merslice = vals[:, :, :, ind_bt]
    elif varname == 'bp':
        ind_bp = mer.lut[qind_bp]
        merslice = vals[:, :, :, ind_bp]
    elif varname == 'br_prime':
        ind_br = mer.lut[qind_br]
        br_slice = vals[:, :, :, ind_br]
        if az is None:
            print("CAUTION: for %s, getting azav from merslice" %varname)
            br_az = np.mean(br_slice, axis=0)
        else:
            br_az = az.vals[:, :, az.lut[qind_br], j]
        merslice = br_slice - br_az.reshape((1, nt, nr))
    elif varname == 'bt_prime':
        ind_bt = mer.lut[qind_bt]
        bt_slice = vals[:, :, :, ind_bt]
        if az is None:
            print("CAUTION: for %s, getting azav from merslice" %varname)
            bt_az = np.mean(bt_slice, axis=0)
        else:
            bt_az = az.vals[:, :, az.lut[qind_bt], j]
        merslice = bt_slice - bt_az
    elif varname == 'bp_prime':
        ind_bp = mer.lut[qind_bp]
        bp_slice = vals[:, :, :, ind_bp]
        if az is None:
            print("CAUTION: for %s, getting azav from merslice" %varname)
            bp_az = np.mean(bp_slice, axis=0)
        else:
            bp_az = az.vals[:, :, az.lut[qind_bp], j]
        merslice = bp_slice - bp_az
            
    # Cylindrical-coordinate belocities
    elif varname == 'bl':
        ind_br, ind_bt = mer.lut[qind_br], mer.lut[qind_bt]
        br_slice = vals[:, :, :, ind_br]
        bt_slice = vals[:, :, :, ind_bt]
        merslice = br_slice*sint + bt_slice*cost       
    elif varname == 'bz':
        ind_br, ind_bt = mer.lut[qind_br], mer.lut[qind_bt]
        br_slice = vals[:, :, :, ind_br]
        bt_slice = vals[:, :, :, ind_bt]
        merslice = br_slice*cost - bt_slice*sint   
    elif varname == 'bl_prime':
        ind_br, ind_bt = mer.lut[qind_br], mer.lut[qind_bt]        
        br_slice = vals[:, :, :, ind_br]
        bt_slice = vals[:, :, :, ind_bt]
        bl_slice = br_slice*sint + bt_slice*cost
        if az is None:
            print("CAUTION: for %s, getting azav from merslice" %varname)
            bl_az = np.mean(bl_slice, axis=0)
        else:
            br_az = az.vals[:, :, az.lut[qind_br], j]
            bt_az = az.vals[:, :, az.lut[qind_bt], j]
            bl_az = br_az*sint + bt_az*cost
        merslice = bl_slice - bl_az
    elif varname == 'bz_prime':
        ind_br, ind_bt = mer.lut[qind_br], mer.lut[qind_bt]        
        br_slice = vals[:, :, :, ind_br]
        bt_slice = vals[:, :, :, ind_bt]
        bz_slice = br_slice*cost - bt_slice*sint  
        if az is None:
            print("CAUTION: for %s, getting azav from merslice" %varname)
            bz_az = np.mean(bz_slice, axis=0)
        else:
            br_az = az.vals[:, :, az.lut[qind_br], j]
            bt_az = az.vals[:, :, az.lut[qind_bt], j]
            bz_az = br_az*sint - bt_az*cost
        merslice = bz_slice - bz_az
    else:
        print("get_merslice(): unknown variable name %s" %varname)
        print("exiting")
        sys.exit()
    return merslice
