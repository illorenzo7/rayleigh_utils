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
from common import thermo_R, thermo_gamma, c_P
from get_parameter import get_parameter
from get_eq import get_eq

def get_merslice(mer, varname, dirname=None, az=None, sh=None, old=False,\
        j=0):
    # Given a Meridional_Slices object (nphi, ntheta, nr, nq, nt), 
    # return the field (nphi, ntheta, nr) associated with [varname] 
    # (take the jth time slice (default j=0))
    # for "primed" variables, the script will use provided AZ_Avgs object az
    # for "sphericallly primed" variables, the script will use provided
    # Shell_Avgs object sh
    # integration weights (from GridInfo) to subtract off various 
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
    vals = a.vals[:, :, :, :, j]


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
               
        qind_omr = var_indices_old['br']
        qind_omt = var_indices_old['bt']
        qind_omp = var_indices_old['bp']
    
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
        vr_az = az.vals[:, :, az.lut[qind_vr], j]/100.
        merslice = vr_slice - vr_az.reshape((1, nt, nr))
    elif varname == 'vt_prime':
        ind_vt = mer.lut[qind_vt]
        vt_slice = vals[:, :, :, ind_vt]/100.
        vt_az = az.vals[:, :, az.lut[qind_vt], j]/100.
        merslice = vt_slice - vt_az
    elif varname == 'vp_prime':
        ind_vp = mer.lut[qind_vp]
        vp_slice = vals[:, :, :, ind_vp]/100.
        merslice = prime(vp_slice)
            
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
        sint = (a.sintheta).reshape((1, a.ntheta, 1))
        cost = (a.costheta).reshape((1, a.ntheta, 1))        
        ind_vr, ind_vt = mer.lut[qind_vr], mer.lut[qind_vt]        
        vr_slice = vals[:, :, :, ind_vr]/100.
        vt_slice = vals[:, :, :, ind_vt]/100.    
        vl_slice = vr_slice*sint + vt_slice*cost
        merslice = prime(vl_slice)
    elif varname == 'vz_prime':
        ind_vr, ind_vt = mer.lut[qind_vr], mer.lut[qind_vt]        
        vr_slice = vals[:, :, :, ind_vr]/100.
        vt_slice = vals[:, :, :, ind_vt]/100.    
        vz_slice = vr_slice*cost - vt_slice*sint  
        merslice = prime(vz_slice)
            
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
        merslice = prime(omr_slice)
    elif varname == 'omt_prime':
        ind_omt = mer.lut[qind_omt]
        omt_slice = vals[:, :, :, ind_omt]
        merslice = prime(omt_slice)
    elif varname == 'omp_prime':
        ind_omp = mer.lut[qind_omp]
        omp_slice = vals[:, :, :, ind_omp]
        merslice = prime(omp_slice)
            
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
    elif varname == 'oml_prime':
        sint = (a.sintheta).reshape((1, a.ntheta, 1))
        cost = (a.costheta).reshape((1, a.ntheta, 1))        
        ind_omr, ind_omt = mer.lut[qind_omr], mer.lut[qind_omt]
        omr_slice = vals[:, :, :, ind_omr]
        omt_slice = vals[:, :, :, ind_omt]
        oml_slice = omr_slice*sint + omt_slice*cost
        merslice = prime(oml_slice)
    elif varname == 'omz_prime':
        ind_omr, ind_omt = mer.lut[qind_omr], mer.lut[qind_omt]
        omr_slice = vals[:, :, :, ind_omr]
        omt_slice = vals[:, :, :, ind_omt]
        omz_slice = omr_slice*cost - omt_slice*sint
        merslice = prime(omz_slice)
            
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
        omr_prime_slice = prime(omr_slice)
        omt_prime_slice = prime(omt_slice)
        omp_prime_slice = prime(omp_slice)
        merslice = omr_prime_slice**2 + omt_prime_slice**2 +\
                omp_prime_slice**2
        
    # Thermodynamic variables: deviations from reference state
    elif varname == 's':
        ind_s = mer.lut[qind_s]
        merslice = vals[:, :, :, ind_s]
    elif varname == 'p':
        ind_p = mer.lut[qind_p]
        merslice = vals[:, :, :, ind_p]
    elif varname == 'rho':
        ind_s, ind_p = mer.lut[qind_s], mer.lut[qind_p]
        s_slice, p_slice = vals[:, :, :, ind_s], vals[:, :, :, ind_p]
        merslice = ref_rho*(p_slice/ref_P/thermo_gamma - s_slice/c_P) 
    elif varname == 't':
        ind_s, ind_p = mer.lut[qind_s], mer.lut[qind_p]
        s_slice, p_slice = vals[:, :, :, ind_s], vals[:, :, :, ind_p]
        merslice = ref_T*(p_slice/ref_P*(1. - 1./thermo_gamma) +\
                s_slice/c_P)
    
    # Thermodynamic variables: deviations from zonal mean
    elif varname == 's_prime':
        ind_s = mer.lut[qind_s]
        s_slice = vals[:, :, :, ind_s]
        merslice = prime(s_slice)
    elif varname == 'p_prime':
        ind_p = mer.lut[qind_p]
        p_slice = vals[:, :, :, ind_p]
        merslice = prime(p_slice)
    elif varname == 'rho_prime':
        ind_s, ind_p = mer.lut[qind_s], mer.lut[qind_p]
        s_slice, p_slice = vals[:, :, :, ind_s], vals[:, :, :, ind_p]
        rho_slice = ref_rho*(p_slice/ref_P/thermo_gamma - s_slice/c_P)
        merslice = prime(rho_slice)
    elif varname == 't_prime':
        ind_s, ind_p = mer.lut[qind_s], mer.lut[qind_p]
        s_slice, p_slice = vals[:, :, :, ind_s], vals[:, :, :, ind_p]
        t_slice = ref_T*(p_slice/ref_P*(1. - 1./thermo_gamma) +\
                s_slice/c_P)
        merslice = prime(t_slice)
            
    # Thermodynamic variables: deviations from spherical mean
    elif varname == 's_prime_sph':
        ind_s = mer.lut[qind_s]
        s_slice = vals[:, :, :, ind_s]
        merslice = prime_sph(s_slice, tw)
    elif varname == 'p_prime_sph':
        ind_p = mer.lut[qind_p]
        p_slice = vals[:, :, :, ind_p]
        merslice = prime_sph(p_slice, tw)
    elif varname == 'rho_prime_sph':
        ind_s, ind_p = mer.lut[qind_s], mer.lut[qind_p]
        s_slice, p_slice = vals[:, :, :, ind_s], vals[:, :, :, ind_p]
        rho_slice = ref_rho*(p_slice/ref_P/thermo_gamma - s_slice/c_P)
        merslice = prime_sph(rho_slice, tw)
    elif varname == 't_prime_sph':
        ind_s, ind_p = mer.lut[qind_s], mer.lut[qind_p]
        s_slice, p_slice = vals[:, :, :, ind_s], vals[:, :, :, ind_p]
        t_slice = ref_T*(p_slice/ref_P*(1. - 1./thermo_gamma) + s_slice/c_P)
        merslice = prime_sph(t_slice, tw)
                                              
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
        merslice = prime(br_slice)
    elif varname == 'bt_prime':
        ind_bt = mer.lut[qind_bt]
        bt_slice = vals[:, :, :, ind_bt]
        merslice = prime(bt_slice)
    elif varname == 'bp_prime':
        ind_bp = mer.lut[qind_bp]
        bp_slice = vals[:, :, :, ind_bp]
        merslice = prime(bp_slice)
            
    # Cylindrical magnetic fields
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
        merslice = prime(bl_slice)
    elif varname == 'bz_prime':
        ind_br, ind_bt = mer.lut[qind_br], mer.lut[qind_bt]
        br_slice = vals[:, :, :, ind_br]
        bt_slice = vals[:, :, :, ind_bt]
        bz_slice = br_slice*cost - bt_slice*sint
        merslice = prime(bz_slice)
    
    # Spherical fluctuating velocity products 
    # (multiply by rho first to get units of pressure)
    elif varname == 'vsq':
        ind_vr, ind_vt, ind_vp = mer.lut[qind_vr], mer.lut[qind_vt],\
                mer.lut[qind_vp]
        vr_slice = vals[:, :, :, ind_vr]
        vt_slice = vals[:, :, :, ind_vt]
        vp_slice = vals[:, :, :, ind_vp]      
        vr_prime_slice = prime(vr_slice)
        vt_prime_slice = prime(vt_slice)
        vp_prime_slice = prime(vp_slice)
        merslice = ref_rho*(vr_prime_slice**2 + vt_prime_slice**2 +\
                vp_prime_slice**2)
    elif varname == 'vrsq':
        ind_vr = mer.lut[qind_vr]
        vr_slice = vals[:, :, :, ind_vr]  
        vr_prime_slice = vr_slice - np.mean(vr_slice, axis=0)
        merslice = ref_rho*vr_prime_slice**2
    elif varname == 'vtsq':
        ind_vt = mer.lut[qind_vt]
        vt_slice = vals[:, :, :, ind_vt]  
        vt_prime_slice = prime(vt_slice)
        merslice = ref_rho*vt_prime_slice**2
    elif varname == 'vpsq':
        ind_vp = mer.lut[qind_vp]
        vp_slice = vals[:, :, :, ind_vp]  
        vp_prime_slice = prime(vp_slice)
        merslice = ref_rho*vp_prime_slice**2        
    elif varname == 'vhsq':
        ind_vt, ind_vp = mer.lut[qind_vt], mer.lut[qind_vp]
        vt_slice = vals[:, :, :, ind_vt]
        vp_slice = vals[:, :, :, ind_vp]      
        vt_prime_slice = prime(vt_slice)
        vp_prime_slice = prime(vp_slice)
        merslice = ref_rho*(vt_prime_slice**2 + vp_prime_slice**2)     
    elif varname == 'vrvp':
        ind_vr, ind_vp = mer.lut[qind_vr], mer.lut[qind_vp]
        vr_slice = vals[:, :, :, ind_vr]
        vp_slice = vals[:, :, :, ind_vp]      
        vr_prime_slice = prime(vr_slice)
        vp_prime_slice = prime(vp_slice)
        merslice = ref_rho*vr_prime_slice*vp_prime_slice      
    elif varname == 'vrvt':
        ind_vr, ind_vt = mer.lut[qind_vr], mer.lut[qind_vt]
        vr_slice = vals[:, :, :, ind_vr]
        vt_slice = vals[:, :, :, ind_vt]     
        vr_prime_slice = prime(vr_slice)
        vt_prime_slice = prime(vt_slice)
        merslice = ref_rho*vr_prime_slice*vt_prime_slice
    elif varname == 'vtvp':
        ind_vt, ind_vp = mer.lut[qind_vt], mer.lut[qind_vp]
        vt_slice = vals[:, :, :, ind_vt]
        vp_slice = vals[:, :, :, ind_vp]       
        vt_prime_slice = prime(vt_slice)
        vp_prime_slice = prime(vp_slice)
        merslice = ref_rho*vt_prime_slice*vp_prime_slice 
        
    # Cylindrical fluctuating velocity products velocity products  
    # Multiplied by rho to get units of pressure
    elif varname == 'vlsq':
        ind_vr, ind_vt = mer.lut[qind_vr], mer.lut[qind_vt]        
        vr_slice = vals[:, :, :, ind_vr]
        vt_slice = vals[:, :, :, ind_vt]   
        vl_slice = vr_slice*sint + vt_slice*cost
        vl_prime_slice = prime(vl_slice)
        merslice = ref_rho*vl_prime_slice**2
    elif varname == 'vzsq':
        ind_vr, ind_vt = mer.lut[qind_vr], mer.lut[qind_vt]        
        vr_slice = vals[:, :, :, ind_vr]
        vt_slice = vals[:, :, :, ind_vt]    
        vz_slice = vr_slice*cost - vt_slice*sint
        vz_prime_slice = prime(vz_slice)
        merslice = ref_rho*vz_prime_slice**2 
    elif varname == 'vmsq':
        ind_vr, ind_vt = mer.lut[qind_vr], mer.lut[qind_vt]        
        vr_slice = vals[:, :, :, ind_vr]
        vt_slice = vals[:, :, :, ind_vt]    
        vr_prime_slice = prime(vr_slice)
        vt_prime_slice = prime(vt_slice)
        merslice = ref_rho*(vr_prime_slice**2 + vt_prime_slice**2)        
    elif varname == 'vlvp':
        ind_vr, ind_vt, ind_vp = mer.lut[qind_vr], mer.lut[qind_vt],\
                mer.lut[qind_vp]
        vr_slice = vals[:, :, :, ind_vr]
        vt_slice = vals[:, :, :, ind_vt]
        vp_slice = vals[:, :, :, ind_vp]      
        vr_prime_slice = prime(vr_slice)
        vt_prime_slice = prime(vt_slice)
        vp_prime_slice = prime(vp_slice)
        vl_prime_slice = vr_prime_slice*sint + vt_prime_slice*cost    
        merslice = ref_rho*vl_prime_slice*vp_prime_slice
    elif varname == 'vlvz':
        ind_vr, ind_vt, ind_vp = mer.lut[qind_vr], mer.lut[qind_vt],\
                mer.lut[qind_vp]
        vr_slice = vals[:, :, :, ind_vr]
        vt_slice = vals[:, :, :, ind_vt] 
        vp_slice = vals[:, :, :, ind_vp]       
        vr_prime_slice = prime(vr_slice)
        vt_prime_slice = prime(vt_slice)
        vp_prime_slice = prime(vp_slice)
        vl_prime_slice = vr_prime_slice*sint + vt_prime_slice*cost   
        vz_prime_slice = vr_prime_slice*cost - vt_prime_slice*sint  
        merslice = ref_rho*vl_prime_slice*vz_prime_slice
    elif varname == 'vpvz':
        ind_vr, ind_vt, ind_vp = mer.lut[qind_vr], mer.lut[qind_vt],\
                mer.lut[qind_vp]
        vr_slice = vals[:, :, :, ind_vr]
        vt_slice = vals[:, :, :, ind_vt]  
        vp_slice = vals[:, :, :, ind_vp]       
        vr_prime_slice = prime(vr_slice)
        vt_prime_slice = prime(vt_slice)
        vp_prime_slice = prime(vp_slice)
        vz_prime_slice = vr_prime_slice*cost - vt_prime_slice*sint  
        merslice = ref_rho*vp_prime_slice*vz_prime_slice 

    # Spherical fluctuating magnetic field products 
    # (multiply by 1/4*pi to get units of magnetic pressure)
    elif varname == 'bsq':
        ind_br, ind_bt, ind_bp = mer.lut[qind_br], mer.lut[qind_bt],\
                mer.lut[qind_bp]
        br_slice = vals[:, :, :, ind_br]
        bt_slice = vals[:, :, :, ind_bt]
        bp_slice = vals[:, :, :, ind_bp]      
        br_prime_slice = prime(br_slice)
        bt_prime_slice = prime(bt_slice)
        bp_prime_slice = prime(bp_slice)
        merslice = (br_prime_slice**2 + bt_prime_slice**2 +\
                bp_prime_slice**2)/(4.*np.pi)
    elif varname == 'brsq':
        ind_br = mer.lut[qind_br]
        br_slice = vals[:, :, :, ind_br]  
        br_prime_slice = prime(br_slice)
        merslice = br_prime_slice**2/(4.*np.pi)
    elif varname == 'btsq':
        ind_bt = mer.lut[qind_bt]
        bt_slice = vals[:, :, :, ind_bt]  
        bt_prime_slice = prime(bt_slice)
        merslice = bt_prime_slice**2/(4.*np.pi)
    elif varname == 'bpsq':
        ind_bp = mer.lut[qind_bp]
        bp_slice = vals[:, :, :, ind_bp]  
        bp_prime_slice = prime(bp_slice)
        merslice = bp_prime_slice**2/(4.*np.pi)
    elif varname == 'bhsq':
        ind_bt, ind_bp = mer.lut[qind_bt], mer.lut[qind_bp]
        bt_slice = vals[:, :, :, ind_bt]
        bp_slice = vals[:, :, :, ind_bp]      
        bt_prime_slice = prime(bt_slice)
        bp_prime_slice = prime(bp_slice)
        merslice = (bt_prime_slice**2 + bp_prime_slice**2)/(4.*np.pi)
    elif varname == 'brbp':
        ind_br, ind_bp = mer.lut[qind_br], mer.lut[qind_bp]
        br_slice = vals[:, :, :, ind_br]
        bp_slice = vals[:, :, :, ind_bp]      
        br_prime_slice = prime(br_slice)
        bp_prime_slice = prime(bp_slice)
        merslice = br_prime_slice*bp_prime_slice/(4.*np.pi)
    elif varname == 'brbt':
        ind_br, ind_bt = mer.lut[qind_br], mer.lut[qind_bt]
        br_slice = vals[:, :, :, ind_br]
        bt_slice = vals[:, :, :, ind_bt]     
        br_prime_slice = prime(br_slice)
        bt_prime_slice = prime(bt_slice)
        merslice = br_prime_slice*bt_prime_slice/(4.*np.pi)
    elif varname == 'btbp':
        ind_bt, ind_bp = mer.lut[qind_bt], mer.lut[qind_bp]
        bt_slice = vals[:, :, :, ind_bt]
        bp_slice = vals[:, :, :, ind_bp]       
        bt_prime_slice = prime(bt_slice)
        bp_prime_slice = prime(bp_slice)
        merslice = bt_prime_slice*bp_prime_slice/(4.*np.pi)
        
    # Cylindrical fluctuating magnetic field products 
    # Multiplied by 1/4*pi to get units of magnetic pressure
    elif varname == 'blsq':
        sint = (a.sintheta).reshape((1, a.ntheta, 1))
        cost = (a.costheta).reshape((1, a.ntheta, 1))        
        ind_br, ind_bt = mer.lut[qind_br], mer.lut[qind_bt]        
        br_slice = vals[:, :, :, ind_br]
        bt_slice = vals[:, :, :, ind_bt]   
        bl_slice = br_slice*sint + bt_slice*cost
        bl_prime_slice = prime(bl_slice)
        merslice = bl_prime_slice**2/(4.*np.pi)
    elif varname == 'bzsq':
        ind_br, ind_bt = mer.lut[qind_br], mer.lut[qind_bt]        
        br_slice = vals[:, :, :, ind_br]
        bt_slice = vals[:, :, :, ind_bt]    
        bz_slice = br_slice*cost - bt_slice*sint
        bz_prime_slice = prime(bz_slice)
        merslice = bz_prime_slice**2/(4.*np.pi)
    elif varname == 'bmsq':
        ind_br, ind_bt = mer.lut[qind_br], mer.lut[qind_bt]        
        br_slice = vals[:, :, :, ind_br]
        bt_slice = vals[:, :, :, ind_bt]    
        br_prime_slice = prime(br_slice)
        bt_prime_slice = prime(bt_slice)
        merslice = (br_prime_slice**2 + bt_prime_slice**2)/(4.*np.pi)
    elif varname == 'blbp':
        ind_br, ind_bt, ind_bp = mer.lut[qind_br], mer.lut[qind_bt],\
                mer.lut[qind_bp]
        br_slice = vals[:, :, :, ind_br]
        bt_slice = vals[:, :, :, ind_bt]
        bp_slice = vals[:, :, :, ind_bp]      
        br_prime_slice = prime(br_slice)
        bt_prime_slice = prime(bt_slice)
        bp_prime_slice = prime(bp_slice)
        bl_prime_slice = br_prime_slice*sint + bt_prime_slice*cost    
        merslice = bl_prime_slice*bp_prime_slice/(4.*np.pi)
    elif varname == 'blbz':
        ind_br, ind_bt, ind_bp = mer.lut[qind_br], mer.lut[qind_bt],\
                mer.lut[qind_bp]
        br_slice = vals[:, :, :, ind_br]
        bt_slice = vals[:, :, :, ind_bt] 
        bp_slice = vals[:, :, :, ind_bp]       
        br_prime_slice = prime(br_slice)
        bt_prime_slice = prime(bt_slice)
        bp_prime_slice = prime(bp_slice)
        bl_prime_slice = br_prime_slice*sint + bt_prime_slice*cost   
        bz_prime_slice = br_prime_slice*cost - bt_prime_slice*sint  
        merslice = bl_prime_slice*bz_prime_slice/(4.*np.pi)
    elif varname == 'bpbz':
        ind_br, ind_bt, ind_bp = mer.lut[qind_br], mer.lut[qind_bt],\
                mer.lut[qind_bp]
        br_slice = vals[:, :, :, ind_br]
        bt_slice = vals[:, :, :, ind_bt]  
        bp_slice = vals[:, :, :, ind_bp]       
        br_prime_slice = prime(br_slice)
        bt_prime_slice = prime(bt_slice)
        bp_prime_slice = prime(bp_slice)
        bz_prime_slice = br_prime_slice*cost - bt_prime_slice*sint  
        merslice = bp_prime_slice*bz_prime_slice/(4.*np.pi)
        
    # Correlations between vertical flow (vr) and temperature/entropy
    elif varname == 'vrs':       
        ind_vr, ind_s = mer.lut[qind_vr], mer.lut[qind_s]
        vr_slice = vals[:, :, :, ind_vr]
        s_slice = vals[:, :, :, ind_s]      
        vr_prime_slice = prime(vr_slice)
        s_prime_slice = prime(s_slice)
        merslice = vr_prime_slice*s_prime_slice      
    elif varname == 'vrt':       
        ind_vr, ind_s, ind_p = mer.lut[qind_vr], mer.lut[qind_s], mer.lut[qind_p]
        vr_slice = vals[:, :, :, ind_vr]
        s_slice = vals[:, :, :, ind_s] 
        p_slice = vals[:, :, :, ind_p]         
        vr_prime_slice = prime(vr_slice)
        s_prime_slice = prime(s_slice)
        p_prime_slice = prime(p_slice)
        t_prime_slice = ref_T*(p_prime_slice/ref_P*(1. - 1./thermo_gamma)\
                + s_prime_slice/c_P)
        merslice = vr_prime_slice*t_prime_slice          
    elif varname == 'vtt':       
        ind_vt, ind_s, ind_p = mer.lut[qind_vt], mer.lut[qind_s], mer.lut[qind_p]
        vt_slice = vals[:, :, :, ind_vt]
        s_slice = vals[:, :, :, ind_s] 
        p_slice = vals[:, :, :, ind_p]         
        vt_prime_slice = prime(vt_slice)
        s_prime_slice = prime(s_slice)
        p_prime_slice = prime(p_slice)
        t_prime_slice = ref_T*(p_prime_slice/ref_P*(1. - 1./thermo_gamma)\
                + s_prime_slice/c_P)          
        merslice = vt_prime_slice*t_prime_slice          
    else:
        print("get_merslice(): unknown variable name %s" %varname)
        print("exiting")
        sys.exit()
    return merslice
