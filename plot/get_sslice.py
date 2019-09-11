# Author: Loren Matilsky
# Created: 01/03/2019
#
# Extremely long and boring script to find fundamental fluid quantities
# or derivative quantities from a shell slice. Takes a shellslice [a] and
# varname in [vr, vt, vp, vl, vz, ...]
# returns the slice for the variable as an array of shape (nphi, ntheta, nr)
#
import numpy as np
import sys, os
sys.path.append(os.environ['rapp'])
from varprops import var_indices, var_indices_old
from common import get_widest_range_file, thermo_R
from get_parameter import get_parameter
from rayleigh_diagnostics import ReferenceState
from reference_tools import equation_coefficients

def get_sslice(a, varname, dirname=None, old=False):
    # Given a shell_slice (nphi, ntheta, nr, nq, nt), return the shell slice
    # (nphi, ntheta, nr) associated with [varname] (use it = 0)
    # for "primed" variables, the script will need shell and/or azimuthal averages
    # to be pre-computed, or will throw an error, so provide "dirname" where these
    # averages are located
    
    vals = a.vals[:, :, :, :, 0]
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
    if (varname == 'vr'):
        ind_vr = a.lut[qind_vr]
        sslice = vals[:, :, :, ind_vr]/100. # measure velocities in m/s
    elif (varname == 'vt'):
        ind_vt = a.lut[qind_vt]
        sslice = vals[:, :, :, ind_vt]/100.
    elif (varname == 'vp'):
        ind_vp = a.lut[qind_vp]
        sslice = vals[:, :, :, ind_vp]/100.
    elif (varname == 'vr_prime'):
        ind_vr = a.lut[qind_vr]
        vr_slice = vals[:, :, :, ind_vr]/100.
        try: 
            azav_file = get_widest_range_file(dirname + '/data/', 'AZ_Avgs')
            di = np.load(dirname + '/data/' + azav_file).item()
            av_vals, lut = di['vals'], di['lut']
            vr_av = av_vals[:, a.inds, lut[qind_vr]]/100.
            sslice = vr_slice - vr_av
        except:
            sslice = vr_slice - np.mean(vr_slice, axis=0)     
    elif (varname == 'vt_prime'):
        ind_vt = a.lut[qind_vt]
        vt_slice = vals[:, :, :, ind_vt]/100.
        try: 
            azav_file = get_widest_range_file(dirname + '/data/', 'AZ_Avgs')
            di = np.load(dirname + '/data/' + azav_file).item()
            av_vals, lut = di['vals'], di['lut']
            vt_av = av_vals[:, a.inds, lut[qind_vt]]/100.
            sslice = vt_slice - vt_av
        except:
            sslice = vt_slice - np.mean(vt_slice, axis=0)             
    elif (varname == 'vp_prime'):
        ind_vp = a.lut[qind_vp]
        vp_slice = vals[:, :, :, ind_vp]/100.
        try: 
            azav_file = get_widest_range_file(dirname + '/data/', 'AZ_Avgs')
            di = np.load(dirname + '/data/' + azav_file).item()
            av_vals, lut = di['vals'], di['lut']
            vp_av = av_vals[:, a.inds, lut[qind_vp]]/100.
            sslice = vp_slice - vp_av
        except:
            sslice = vp_slice - np.mean(vp_slice, axis=0)
            
    # Cylindrical velocities
    elif (varname == 'vl'):
        sint = (a.sintheta).reshape((1, a.ntheta, 1))
        cost = (a.costheta).reshape((1, a.ntheta, 1))        
        ind_vr, ind_vt = a.lut[qind_vr], a.lut[qind_vt]
        vr_slice = vals[:, :, :, ind_vr]/100.
        vt_slice = vals[:, :, :, ind_vt]/100.
        sslice = vr_slice*sint + vt_slice*cost       
    elif (varname == 'vz'):
        sint = (a.sintheta).reshape((1, a.ntheta, 1))
        cost = (a.costheta).reshape((1, a.ntheta, 1))        
        ind_vr, ind_vt = a.lut[qind_vr], a.lut[qind_vt]
        vr_slice = vals[:, :, :, ind_vr]/100.
        vt_slice = vals[:, :, :, ind_vt]/100.
        sslice = vr_slice*cost - vt_slice*sint   
    elif (varname == 'vl_prime'):
        sint = (a.sintheta).reshape((1, a.ntheta, 1))
        cost = (a.costheta).reshape((1, a.ntheta, 1))        
        ind_vr, ind_vt = a.lut[qind_vr], a.lut[qind_vt]        
        vr_slice = vals[:, :, :, ind_vr]/100.
        vt_slice = vals[:, :, :, ind_vt]/100.    
        vl_slice = vr_slice*sint + vt_slice*cost
        try: 
            azav_file = get_widest_range_file(dirname + '/data/', 'AZ_Avgs')
            di = np.load(dirname + '/data/' + azav_file).item()
            av_vals, lut = di['vals'], di['lut']
            vr_av = av_vals[:, a.inds, lut[qind_vr]]/100.
            vt_av = av_vals[:, a.inds, lut[qind_vt]]/100. 
            vl_av = vr_av*sint + vt_av*cost 
            sslice = vl_slice - vl_av
        except:
            sslice = vl_slice - np.mean(vl_slice, axis=0)
    elif (varname == 'vz_prime'):
        sint = (a.sintheta).reshape((1, a.ntheta, 1))
        cost = (a.costheta).reshape((1, a.ntheta, 1))        
        ind_vr, ind_vt = a.lut[qind_vr], a.lut[qind_vt]        
        vr_slice = vals[:, :, :, ind_vr]/100.
        vt_slice = vals[:, :, :, ind_vt]/100.    
        vz_slice = vr_slice*cost - vt_slice*sint  
        try: 
            azav_file = get_widest_range_file(dirname + '/data/', 'AZ_Avgs')
            di = np.load(dirname + '/data/' + azav_file).item()
            av_vals, lut = di['vals'], di['lut']
            vr_av = av_vals[:, a.inds, lut[qind_vr]]/100.
            vt_av = av_vals[:, a.inds, lut[qind_vt]]/100. 
            vz_av = vl_av*cost - vt_av*sint  
            sslice = vz_slice - vz_av
        except:
            sslice = vz_slice - np.mean(vz_slice, axis=0)            
            
    # Spherical vorticities
    elif (varname == 'omr'):
        ind_omr = a.lut[qind_omr]
        sslice = vals[:, :, :, ind_omr]
    elif (varname == 'omt'):
        ind_omt = a.lut[qind_omt]
        sslice = vals[:, :, :, ind_omt]
    elif (varname == 'omp'):
        ind_omp = a.lut[qind_omp]
        sslice = vals[:, :, :, ind_omp]
    elif (varname == 'omr_prime'):
        ind_omr = a.lut[qind_omr]
        omr_slice = vals[:, :, :, ind_omr]
        try: 
            azav_file = get_widest_range_file(dirname + '/data/', 'AZ_Avgs')
            di = np.load(dirname + '/data/' + azav_file).item()
            av_vals, lut = di['vals'], di['lut']
            omr_av = av_vals[:, a.inds, lut[qind_omr]]/100.
            sslice = omr_slice - omr_av
        except:
            sslice = omr_slice - np.mean(omr_slice, axis=0)
    elif (varname == 'omt_prime'):
        ind_omt = a.lut[qind_omt]
        omt_slice = vals[:, :, :, ind_omt]
        try: 
            azav_file = get_widest_range_file(dirname + '/data/', 'AZ_Avgs')
            di = np.load(dirname + '/data/' + azav_file).item()
            av_vals, lut = di['vals'], di['lut']
            omt_av = av_vals[:, a.inds, lut[qind_omt]]/100.
            sslice = omt_slice - omt_av
        except:
            sslice = omt_slice - np.mean(omt_slice, axis=0)    
    elif (varname == 'omp_prime'):
        ind_omp = a.lut[qind_omp]
        omp_slice = vals[:, :, :, ind_omp]
        try: 
            azav_file = get_widest_range_file(dirname + '/data/', 'AZ_Avgs')
            di = np.load(dirname + '/data/' + azav_file).item()
            av_vals, lut = di['vals'], di['lut']
            omp_av = av_vals[:, a.inds, lut[qind_omp]]/100.
            sslice = omp_slice - omp_av
        except:
            sslice = omp_slice - np.mean(omp_slice, axis=0)  
            
    # Cylindrical vorticities
    elif (varname == 'oml'):
        sint = (a.sintheta).reshape((1, a.ntheta, 1))
        cost = (a.costheta).reshape((1, a.ntheta, 1))        
        ind_omr, ind_omt = a.lut[qind_omr], a.lut[qind_omt]
        omr_slice = vals[:, :, :, ind_omr]
        omt_slice = vals[:, :, :, ind_omt]
        sslice = omr_slice*sint + omt_slice*cost       
    elif (varname == 'omz'):
        sint = (a.sintheta).reshape((1, a.ntheta, 1))
        cost = (a.costheta).reshape((1, a.ntheta, 1))        
        ind_omr, ind_omt = a.lut[qind_omr], a.lut[qind_omt]
        omr_slice = vals[:, :, :, ind_omr]
        omt_slice = vals[:, :, :, ind_omt]
        sslice = omr_slice*cost - omt_slice*sint
    elif (varname == 'oml_prime'):
        sint = (a.sintheta).reshape((1, a.ntheta, 1))
        cost = (a.costheta).reshape((1, a.ntheta, 1))        
        ind_omr, ind_omt = a.lut[qind_omr], a.lut[qind_omt]
        omr_slice = vals[:, :, :, ind_omr]
        omt_slice = vals[:, :, :, ind_omt]
        oml_slice = omr_slice*sint + omt_slice*cost
        try: 
            azav_file = get_widest_range_file(dirname + '/data/', 'AZ_Avgs')
            di = np.load(dirname + '/data/' + azav_file).item()
            av_vals, lut = di['vals'], di['lut']
            omr_av = av_vals[:, a.inds, lut[qind_omr]]
            omt_av = av_vals[:, a.inds, lut[qind_omt]]
            oml_av = omr_av*sint + omt_av*cost
            sslice = oml_slice - oml_av            
        except:
            sslice = oml_slice - np.mean(oml_slice, axis=0)           
    elif (varname == 'omz_prime'):
        sint = (a.sintheta).reshape((1, a.ntheta, 1))
        cost = (a.costheta).reshape((1, a.ntheta, 1))        
        ind_omr, ind_omt = a.lut[qind_omr], a.lut[qind_omt]
        omr_slice = vals[:, :, :, ind_omr]
        omt_slice = vals[:, :, :, ind_omt]
        omz_slice = omr_slice*cost - omt_slice*sint
        try: 
            azav_file = get_widest_range_file(dirname + '/data/', 'AZ_Avgs')
            di = np.load(dirname + '/data/' + azav_file).item()
            av_vals, lut = di['vals'], di['lut']
            omr_av = av_vals[:, a.inds, lut[qind_omr]]
            omt_av = av_vals[:, a.inds, lut[qind_omt]]
            omz_av = omr_av*cost - omt_av*sint
            sslice = omz_slice - omz_av            
        except:
            sslice = omz_slice - np.mean(omz_slice, axis=0) 
            
    # Enstrophy (total and fluctuating)
    elif (varname == 'omsq'):
        ind_omr, ind_omt, ind_omp = a.lut[qind_omr], a.lut[qind_omt], a.lut[qind_omp]
        sslice = vals[:, :, :, ind_omr]**2 + vals[:, :, :, ind_omt]**2 + vals[:, :, :, ind_omp]**2 
    elif (varname == 'omsq_fluc'):
        ind_omr, ind_omt, ind_omp = a.lut[qind_omr], a.lut[qind_omt], a.lut[qind_omp]
        omr_slice = vals[:, :, :, ind_omr]
        omt_slice = vals[:, :, :, ind_omt]
        omp_slice = vals[:, :, :, ind_omp]
        try: 
            azav_file = get_widest_range_file(dirname + '/data/', 'AZ_Avgs')
            di = np.load(dirname + '/data/' + azav_file).item()
            av_vals, lut = di['vals'], di['lut']
            omr_av = av_vals[:, a.inds, lut[qind_omr]]
            omt_av = av_vals[:, a.inds, lut[qind_omt]]
            omp_av = av_vals[:, a.inds, lut[qind_omp]]
            omr_prime_slice = omr_slice - omr_av
            omt_prime_slice = omt_slice - omt_av    
            omp_prime_slice = omp_slice - omp_av                
        except:
            omr_prime_slice = omr_slice - np.mean(omr_slice, axis=0) 
            omt_prime_slice = omt_slice - np.mean(omt_slice, axis=0) 
            omp_prime_slice = omp_slice - np.mean(omp_slice, axis=0)
        sslice = omr_prime_slice**2 + omt_prime_slice**2 + omp_prime_slice**2
        
    # Thermodynamic variables: deviations from reference state
    elif (varname == 's'):
        ind_s = a.lut[qind_s]
        sslice = vals[:, :, :, ind_s]
    elif (varname == 'p'):
        ind_p = a.lut[qind_p]
        sslice = vals[:, :, :, ind_p]
    elif (varname == 'rho'):
        # must have the reference state for these derivative variables
        ref = ReferenceState(dirname + '/reference', '')
        cp = get_parameter(dirname, 'pressure_specific_heat')
        ref_p = (ref.pressure)[a.inds]
        ref_rho = (ref.density)[a.inds]
        ind_s, ind_p = a.lut[qind_s], a.lut[qind_p]
        s_slice, p_slice = vals[:, :, :, ind_s], vals[:, :, :, ind_p]
        try:
            poly_n = get_paramater(dirname, 'poly_n') # must check "poly_n" later
            poly_gamma = 1. + 1./poly_n
            sslice = ref_rho*(p_slice/ref_p/poly_gamma - s_slice/cp)
        except:
            poly_gamma = 5./3.
            sslice = ref_rho*(p_slice/ref_p/poly_gamma - s_slice/cp)            
    elif (varname == 't'):
        # must have the reference state for these derivative variables
        # assume ideal gas law
        ref = ReferenceState(dirname + '/reference', '')
        cp = get_parameter(dirname, 'pressure_specific_heat')
        ref_p = (ref.pressure)[a.inds]
        ref_temp = (ref.temperature)[a.inds]
        ind_s, ind_p = a.lut[qind_s], a.lut[qind_p]
        s_slice, p_slice = vals[:, :, :, ind_s], vals[:, :, :, ind_p]
        try:
            poly_n = get_paramater(dirname, 'poly_n') # must check "poly_n" later
            poly_gamma = 1. + 1./poly_n
            sslice = ref_temp*(p_slice/ref_p*(1. - 1./poly_gamma) + s_slice/cp)   
        except:
            poly_gamma = 5./3.
            sslice = ref_temp*(p_slice/ref_p*(1. - 1./poly_gamma) + s_slice/cp)                 
    
    # Thermodynamic variables: deviations from zonal mean
    elif (varname == 's_prime'):
        ind_s = a.lut[qind_s]
        s_slice = vals[:, :, :, ind_s]
        try: 
            azav_file = get_widest_range_file(dirname + '/data/', 'AZ_Avgs')
            di = np.load(dirname + '/data/' + azav_file).item()
            av_vals, lut = di['vals'], di['lut']
            s_av = av_vals[:, a.inds, lut[qind_s]]
            sslice = s_slice - s_av
        except:
            sslice = s_slice - np.mean(s_slice, axis=0)
    elif (varname == 'p_prime'):
        ind_p = a.lut[qind_p]
        p_slice = vals[:, :, :, ind_p]
        try: 
            azav_file = get_widest_range_file(dirname + '/data/', 'AZ_Avgs')
            di = np.load(dirname + '/data/' + azav_file).item()
            av_vals, lut = di['vals'], di['lut']
            p_av = av_vals[:, a.inds, lut[qind_p]]
            sslice = p_slice - p_av
        except:
            sslice = p_slice - np.mean(p_slice, axis=0)            
    elif (varname == 'rho_prime'):
        # Must have reference state for rho, t
        ref = ReferenceState(dirname + '/reference', '')
        cp = get_parameter(dirname, 'pressure_specific_heat')
        ref_p = (ref.pressure)[a.inds]
        ref_rho = (ref.density)[a.inds]
        ind_s, ind_p = a.lut[qind_s], a.lut[qind_p]
        s_slice, p_slice = vals[:, :, :, ind_s], vals[:, :, :, ind_p]
        try:
            poly_n = get_paramater(dirname, 'poly_n') # must check "poly_n" later
            poly_gamma = 1. + 1./poly_n
            rho_slice = ref_rho*(p_slice/ref_p/poly_gamma - s_slice/cp)
        except:
            poly_gamma = 5./3.
            rho_slice = ref_rho*(p_slice/ref_p/poly_gamma - s_slice/cp)
        try: 
            azav_file = get_widest_range_file(dirname + '/data/', 'AZ_Avgs')
            di = np.load(dirname + '/data/' + azav_file).item()
            av_vals, lut = di['vals'], di['lut']
            s_av = av_vals[:, a.inds, lut[qind_s]]
            p_av = av_vals[:, a.inds, lut[qind_p]]
            rho_av = ref_rho*(p_av/ref_p/poly_gamma - s_av/cp)
            sslice = rho_slice - rho_av            
        except:
            sslice = rho_slice - np.mean(rho_slice, axis=0)           
    elif (varname == 't_prime'):
        # Must have reference state for rho, t
        ref = ReferenceState(dirname + '/reference', '')
        cp = get_parameter(dirname, 'pressure_specific_heat')
        ref_p = (ref.pressure)[a.inds]
        ref_temp = (ref.temperature)[a.inds]
        ind_s, ind_p = a.lut[qind_s], a.lut[qind_p]
        s_slice, p_slice = vals[:, :, :, ind_s], vals[:, :, :, ind_p]
        try:
            poly_n = get_paramater(dirname, 'poly_n') # must check "poly_n" later
            poly_gamma = 1. + 1./poly_n
            t_slice = ref_temp*(p_slice/ref_p*(1. - 1./poly_gamma) + s_slice/cp)
        except:
            poly_gamma = 5./3.
            t_slice = ref_temp*(p_slice/ref_p*(1. - 1./poly_gamma) + s_slice/cp)
        try: 
            azav_file = get_widest_range_file(dirname + '/data/', 'AZ_Avgs')
            di = np.load(dirname + '/data/' + azav_file).item()
            av_vals, lut = di['vals'], di['lut']
            s_av = av_vals[:, a.inds, lut[qind_s]]
            p_av = av_vals[:, a.inds, lut[qind_p]]
            t_av = ref_temp*(p_av/ref_p*(1. - 1./poly_gamma) + s_av/cp)
            sslice = t_slice - t_av            
        except:
            sslice = t_slice - np.mean(t_slice, axis=0)         
            
    # Thermodynamic variables: deviations from spherical mean
    elif (varname == 's_prime_sph'):
        ind_s = a.lut[qind_s]
        s_slice = vals[:, :, :, ind_s]
        try: 
            shav_file = get_widest_range_file(dirname + '/data/', 'Shell_Avgs')
            di = np.load(dirname + '/data/' + shav_file).item()
            av_vals, lut = di['vals'], di['lut']
            s_av = av_vals[a.inds, lut[qind_s]]
            sslice = s_slice - s_av
        except:
            sslice = s_slice - np.mean(s_slice) # Should probably fix this sloppy average
    elif (varname == 'p_prime_sph'):
        ind_p = a.lut[qind_p]
        p_slice = vals[:, :, :, ind_p]
        try: 
            shav_file = get_widest_range_file(dirname + '/data/', 'Shell_Avgs')
            di = np.load(dirname + '/data/' + shav_file).item()
            av_vals, lut = di['vals'], di['lut']
            p_av = av_vals[a.inds, lut[qind_p]]
            sslice = p_slice - p_av
        except:
            sslice = p_slice - np.mean(p_slice)  # Should probably fix this sloppy average          
    elif (varname == 'rho_prime_sph'):
        # must have the reference state for these derivative variables
        ref = ReferenceState(dirname + '/reference', '')
        cp = get_parameter(dirname, 'pressure_specific_heat')
        ref_p = (ref.pressure)[a.inds]
        ref_rho = (ref.density)[a.inds]
        ind_s, ind_p = a.lut[qind_s], a.lut[qind_p]
        s_slice, p_slice = vals[:, :, :, ind_s], vals[:, :, :, ind_p]
        try:
            poly_n = get_paramater(dirname, 'poly_n') # must check "poly_n" later
            poly_gamma = 1. + 1./poly_n
            rho_slice = ref_rho*(p_slice/ref_p/poly_gamma - s_slice/cp)
        except:
            poly_gamma = 5./3.
            rho_slice = ref_rho*(p_slice/ref_p/poly_gamma - s_slice/cp)
        try: 
            shav_file = get_widest_range_file(dirname + '/data/', 'Shell_Avgs')
            di = np.load(dirname + '/data/' + shav_file).item()
            av_vals, lut = di['vals'], di['lut']
            s_av = av_vals[a.inds, lut[qind_s]]
            p_av = av_vals[a.inds, lut[qind_p]]
            rho_av = ref_rho*(p_av/ref_p/poly_gamma - s_av/cp)
            sslice = rho_slice - rho_av            
        except:
            sslice = rho_slice - np.mean(rho_slice) # Should probably fix this sloppy average          
    elif (varname == 't_prime_sph'):
        # must have the reference state for these derivative variables
        # assume ideal gas law                                     
        ref = ReferenceState(dirname + '/reference', '')
        cp = get_parameter(dirname, 'pressure_specific_heat')
        ref_p = (ref.pressure)[a.inds]
        ref_temp = (ref.temperature)[a.inds]
        ind_s, ind_p = a.lut[qind_s], a.lut[qind_p]
        s_slice, p_slice = vals[:, :, :, ind_s], vals[:, :, :, ind_p]
        try:
            poly_n = get_paramater(dirname, 'poly_n') # must check "poly_n" later
            poly_gamma = 1. + 1./poly_n
            t_slice = ref_temp*(p_sliceref_p*(1. - 1./poly_gamma) + s_slice/cp)
        except:
            poly_gamma = 5./3.
            t_slice = ref_temp*(p_slice/ref_p*(1. - 1./poly_gamma) + s_slice/cp) 
        try: 
            shav_file = get_widest_range_file(dirname + '/data/', 'Shell_Avgs')
            di = np.load(dirname + '/data/' + shav_file).item()
            av_vals, lut = di['vals'], di['lut']
            s_av = av_vals[a.inds, lut[qind_s]]
            p_av = av_vals[a.inds, lut[qind_p]]
            t_av = ref_temp*(p_av/ref_p*(1. - 1./poly_gamma) + s_av/cp)
            sslice = t_slice - t_av            
        except:
            sslice = t_slice - np.mean(t_slice) # Should probably fix this sloppy average
                                              
    # Spherical magnetic fields
    elif (varname == 'br'):
        ind_br = a.lut[qind_br]
        sslice = vals[:, :, :, ind_br]
    elif (varname == 'bt'):
        ind_bt = a.lut[qind_bt]
        sslice = vals[:, :, :, ind_bt]
    elif (varname == 'bp'):
        ind_bp = a.lut[qind_bp]
        sslice = vals[:, :, :, ind_bp]
    elif (varname == 'br_prime'):
        ind_br = a.lut[qind_br]
        br_slice = vals[:, :, :, ind_br]
        print(np.max(br_slice))
        try: 
            azav_file = get_widest_range_file(dirname + '/data/', 'AZ_Avgs')
            di = np.load(dirname + '/data/' + azav_file).item()
            av_vals, lut = di['vals'], di['lut']
            br_av = av_vals[:, a.inds, lut[qind_br]]
            sslice = br_slice - br_av
        except:
            sslice = br_slice - np.mean(br_slice, axis=0)
    elif (varname == 'bt_prime'):
        ind_bt = a.lut[qind_bt]
        bt_slice = vals[:, :, :, ind_bt]
        try: 
            azav_file = get_widest_range_file(dirname + '/data/', 'AZ_Avgs')
            di = np.load(dirname + '/data/' + azav_file).item()
            av_vals, lut = di['vals'], di['lut']
            bt_av = av_vals[:, a.inds, lut[qind_bt]]
            sslice = bt_slice - bt_av
        except:
            sslice = bt_slice - np.mean(bt_slice, axis=0)      
    elif (varname == 'bp_prime'):
        ind_bp = a.lut[qind_bp]
        bp_slice = vals[:, :, :, ind_bp]
        try: 
            azav_file = get_widest_range_file(dirname + '/data/', 'AZ_Avgs')
            di = np.load(dirname + '/data/' + azav_file).item()
            av_vals, lut = di['vals'], di['lut']
            bp_av = av_vals[:, a.inds, lut[qind_bp]]
            sslice = bp_slice - bp_av
        except:
            sslice = bp_slice - np.mean(bp_slice, axis=0)
            
    # Cylindrical magnetic fields
    elif (varname == 'bl'):
        ind_br, ind_bt = a.lut[qind_br], a.lut[qind_bt]
        br_slice = vals[:, :, :, ind_br]
        bt_slice = vals[:, :, :, ind_bt]
        sint = (a.sintheta).reshape((1, a.ntheta, 1))
        cost = (a.costheta).reshape((1, a.ntheta, 1))
        sslice = br_slice*sint + bt_slice*cost       
    elif (varname == 'bz'):
        ind_br, ind_bt = a.lut[qind_br], a.lut[qind_bt]
        br_slice = vals[:, :, :, ind_br]
        bt_slice = vals[:, :, :, ind_bt]
        sint = (a.sintheta).reshape((1, a.ntheta, 1))
        cost = (a.costheta).reshape((1, a.ntheta, 1))
        sslice = br_slice*cost - bt_slice*sint           
    elif (varname == 'bl_prime'):
        ind_br, ind_bt = a.lut[qind_br], a.lut[qind_bt]
        br_slice = vals[:, :, :, ind_br]
        bt_slice = vals[:, :, :, ind_bt]
        sint = (a.sintheta).reshape((1, a.ntheta, 1))
        cost = (a.costheta).reshape((1, a.ntheta, 1))
        bl_slice = br_slice*sint + bt_slice*cost
        try: 
            azav_file = get_widest_range_file(dirname + '/data/', 'AZ_Avgs')
            di = np.load(dirname + '/data/' + azav_file).item()
            av_vals, lut = di['vals'], di['lut']
            br_av = av_vals[:, a.inds, lut[qind_br]]
            bt_av = av_vals[:, a.inds, lut[qind_bt]]
            bl_av = br_av*sint + bt_av*cost
            sslice = bl_slice - bl_av            
        except:
            sslice = bl_slice - np.mean(bl_slice, axis=0)      
    elif (varname == 'bz_prime'):
        ind_br, ind_bt = a.lut[qind_br], a.lut[qind_bt]
        br_slice = vals[:, :, :, ind_br]
        bt_slice = vals[:, :, :, ind_bt]
        sint = (a.sintheta).reshape((1, a.ntheta, 1))
        cost = (a.costheta).reshape((1, a.ntheta, 1))
        bz_slice = br_slice*cost - bt_slice*sint
        try: 
            azav_file = get_widest_range_file(dirname + '/data/', 'AZ_Avgs')
            di = np.load(dirname + '/data/' + azav_file).item()
            av_vals, lut = di['vals'], di['lut']
            br_av = av_vals[:, a.inds, lut[qind_br]]
            bt_av = av_vals[:, a.inds, lut[qind_bt]]
            bz_av = br_av*cost - bt_av*sint
            sslice = bz_slice - bz_av            
        except:
            sslice = bz_slice - np.mean(bz_slice, axis=0)         
    
    # Spherical fluctuating velocity products 
    # (multiply by rho first to get units of pressure)
    elif (varname == 'vsq'):
        ref = ReferenceState(dirname + '/reference', '')
        prefactor = ref.density[a.inds]
        ind_vr, ind_vt, ind_vp = a.lut[qind_vr], a.lut[qind_vt], a.lut[qind_vp]
        vr_slice = vals[:, :, :, ind_vr]
        vt_slice = vals[:, :, :, ind_vt]
        vp_slice = vals[:, :, :, ind_vp]      
        try: 
            azav_file = get_widest_range_file(dirname + '/data/', 'AZ_Avgs')
            di = np.load(dirname + '/data/' + azav_file).item()
            av_vals, lut = di['vals'], di['lut']
            vr_av = av_vals[:, a.inds, lut[qind_vr]]
            vt_av = av_vals[:, a.inds, lut[qind_vt]]
            vp_av = av_vals[:, a.inds, lut[qind_vp]]
            vr_prime_slice = vr_slice - vr_av
            vt_prime_slice = vt_slice - vt_av
            vp_prime_slice = vp_slice - vp_av
        except:
            vr_prime_slice = vr_slice - np.mean(vr_slice, axis=0)
            vt_prime_slice = vt_slice - np.mean(vt_slice, axis=0)
            vp_prime_slice = vp_slice - np.mean(vp_slice, axis=0)
        sslice = prefactor*(vr_prime_slice**2 + vt_prime_slice**2 + vp_prime_slice**2)
    elif (varname == 'vrsq'):
        ref = ReferenceState(dirname + '/reference', '')
        prefactor = ref.density[a.inds]
        ind_vr = a.lut[qind_vr]
        vr_slice = vals[:, :, :, ind_vr]  
        try: 
            azav_file = get_widest_range_file(dirname + '/data/', 'AZ_Avgs')
            di = np.load(dirname + '/data/' + azav_file).item()
            av_vals, lut = di['vals'], di['lut']
            vr_av = av_vals[:, a.inds, lut[qind_vr]]
            vr_prime_slice = vr_slice - vr_av
        except:
            vr_prime_slice = vr_slice - np.mean(vr_slice, axis=0)
        sslice = prefactor*vr_prime_slice**2
    elif (varname == 'vtsq'):
        ref = ReferenceState(dirname + '/reference', '')
        prefactor = ref.density[a.inds]
        ind_vt = a.lut[qind_vt]
        vt_slice = vals[:, :, :, ind_vt]  
        try: 
            azav_file = get_widest_range_file(dirname + '/data/', 'AZ_Avgs')
            di = np.load(dirname + '/data/' + azav_file).item()
            av_vals, lut = di['vals'], di['lut']
            vt_av = av_vals[:, a.inds, lut[qind_vt]]
            vt_prime_slice = vt_slice - vt_av
        except:
            vt_prime_slice = vt_slice - np.mean(vt_slice, axis=0)
        sslice = prefactor*vt_prime_slice**2
    elif (varname == 'vpsq'):
        ref = ReferenceState(dirname + '/reference', '')
        prefactor = ref.density[a.inds]
        ind_vp = a.lut[qind_vp]
        vp_slice = vals[:, :, :, ind_vp]  
        try: 
            azav_file = get_widest_range_file(dirname + '/data/', 'AZ_Avgs')
            di = np.load(dirname + '/data/' + azav_file).item()
            av_vals, lut = di['vals'], di['lut']
            vp_av = av_vals[:, a.inds, lut[qind_vp]]
            vp_prime_slice = vp_slice - vp_av
        except:
            vp_prime_slice = vp_slice - np.mean(vp_slice, axis=0)
        sslice = prefactor*vp_prime_slice**2        
    elif (varname == 'vhsq'):
        ref = ReferenceState(dirname + '/reference', '')
        prefactor = ref.density[a.inds]        
        ind_vt, ind_vp = a.lut[qind_vt], a.lut[qind_vp]
        vt_slice = vals[:, :, :, ind_vt]
        vp_slice = vals[:, :, :, ind_vp]      
        try: 
            azav_file = get_widest_range_file(dirname + '/data/', 'AZ_Avgs')
            di = np.load(dirname + '/data/' + azav_file).item()
            av_vals, lut = di['vals'], di['lut']
            vt_av = av_vals[:, a.inds, lut[qind_vt]]
            vp_av = av_vals[:, a.inds, lut[qind_vp]]
            vt_prime_slice = vt_slice - vt_av
            vp_prime_slice = vp_slice - vp_av
        except:
            vt_prime_slice = vt_slice - np.mean(vt_slice, axis=0)
            vp_prime_slice = vp_slice - np.mean(vp_slice, axis=0)
        sslice = prefactor*(vt_prime_slice**2 + vp_prime_slice**2)     
    elif (varname == 'vrvp'):
        ref = ReferenceState(dirname + '/reference', '')
        prefactor = ref.density[a.inds]        
        ind_vr, ind_vp = a.lut[qind_vr], a.lut[qind_vp]
        vr_slice = vals[:, :, :, ind_vr]
        vp_slice = vals[:, :, :, ind_vp]      
        try: 
            azav_file = get_widest_range_file(dirname + '/data/', 'AZ_Avgs')
            di = np.load(dirname + '/data/' + azav_file).item()
            av_vals, lut = di['vals'], di['lut']
            vr_av = av_vals[:, a.inds, lut[qind_vr]]
            vp_av = av_vals[:, a.inds, lut[qind_vp]]
            vr_prime_slice = vr_slice - vr_av
            vp_prime_slice = vp_slice - vp_av
        except:
            vr_prime_slice = vr_slice - np.mean(vr_slice, axis=0)
            vp_prime_slice = vp_slice - np.mean(vp_slice, axis=0)    
        sslice = prefactor*vr_prime_slice*vp_prime_slice      
    elif (varname == 'vrvt'):
        ref = ReferenceState(dirname + '/reference', '')
        prefactor = ref.density[a.inds]        
        ind_vr, ind_vt = a.lut[qind_vr], a.lut[qind_vt]
        vr_slice = vals[:, :, :, ind_vr]
        vt_slice = vals[:, :, :, ind_vt]     
        try: 
            azav_file = get_widest_range_file(dirname + '/data/', 'AZ_Avgs')
            di = np.load(dirname + '/data/' + azav_file).item()
            av_vals, lut = di['vals'], di['lut']
            vr_av = av_vals[:, a.inds, lut[qind_vr]]
            vt_av = av_vals[:, a.inds, lut[qind_vt]]
            vr_prime_slice = vr_slice - vr_av
            vt_prime_slice = vt_slice - vt_av
        except:
            vr_prime_slice = vr_slice - np.mean(vr_slice, axis=0)
            vt_prime_slice = vt_slice - np.mean(vt_slice, axis=0)
        sslice = prefactor*vr_prime_slice*vt_prime_slice
    elif (varname == 'vtvp'):
        ref = ReferenceState(dirname + '/reference', '')
        prefactor = ref.density[a.inds]        
        ind_vt, ind_vp = a.lut[qind_vt], a.lut[qind_vp]
        vt_slice = vals[:, :, :, ind_vt]
        vp_slice = vals[:, :, :, ind_vp]       
        try: 
            azav_file = get_widest_range_file(dirname + '/data/', 'AZ_Avgs')
            di = np.load(dirname + '/data/' + azav_file).item()
            av_vals, lut = di['vals'], di['lut']
            vt_av = av_vals[:, a.inds, lut[qind_vt]]
            vp_av = av_vals[:, a.inds, lut[qind_vp]]
            vt_prime_slice = vt_slice - vt_av
            vp_prime_slice = vp_slice - vp_av
        except:
            vt_prime_slice = vt_slice - np.mean(vt_slice, axis=0)
            vp_prime_slice = vp_slice - np.mean(vp_slice, axis=0)    
        sslice = prefactor*vt_prime_slice*vp_prime_slice 
        
    # Cylindrical fluctuating velocity products velocity products  
    # Multiplied by rho to get units of pressure
    elif (varname == 'vlsq'):
        ref = ReferenceState(dirname + '/reference', '')
        prefactor = ref.density[a.inds]          
        sint = (a.sintheta).reshape((1, a.ntheta, 1))
        cost = (a.costheta).reshape((1, a.ntheta, 1))        
        ind_vr, ind_vt = a.lut[qind_vr], a.lut[qind_vt]        
        vr_slice = vals[:, :, :, ind_vr]
        vt_slice = vals[:, :, :, ind_vt]   
        vl_slice = vr_slice*sint + vt_slice*cost
        try: 
            azav_file = get_widest_range_file(dirname + '/data/', 'AZ_Avgs')
            di = np.load(dirname + '/data/' + azav_file).item()
            av_vals, lut = di['vals'], di['lut']
            vr_av = av_vals[:, a.inds, lut[qind_vr]]
            vt_av = av_vals[:, a.inds, lut[qind_vt]]
            vl_av = vr_av*sint + vt_av*cost 
            vl_prime_slice = vl_slice - vl_av
        except:
            vl_prime_slice = vl_slice - np.mean(vl_slice, axis=0) 
        sslice = prefactor*vl_prime_slice**2
    elif (varname == 'vzsq'):
        ref = ReferenceState(dirname + '/reference', '')
        prefactor = ref.density[a.inds]          
        sint = (a.sintheta).reshape((1, a.ntheta, 1))
        cost = (a.costheta).reshape((1, a.ntheta, 1))        
        ind_vr, ind_vt = a.lut[qind_vr], a.lut[qind_vt]        
        vr_slice = vals[:, :, :, ind_vr]
        vt_slice = vals[:, :, :, ind_vt]    
        vz_slice = vr_slice*cost - vt_slice*sint
        try: 
            azav_file = get_widest_range_file(dirname + '/data/', 'AZ_Avgs')
            di = np.load(dirname + '/data/' + azav_file).item()
            av_vals, lut = di['vals'], di['lut']
            vr_av = av_vals[:, a.inds, lut[qind_vr]]
            vt_av = av_vals[:, a.inds, lut[qind_vt]]
            vz_av = vr_av*cost - vt_av*sint 
            vz_prime_slice = vz_slice - vz_av
        except:
            vz_prime_slice = vz_slice - np.mean(vz_slice, axis=0) 
        sslice = prefactor*vz_prime_slice**2 
    elif (varname == 'vpolsq'):
        ref = ReferenceState(dirname + '/reference', '')
        prefactor = ref.density[a.inds]                 
        ind_vr, ind_vt = a.lut[qind_vr], a.lut[qind_vt]        
        vr_slice = vals[:, :, :, ind_vr]
        vt_slice = vals[:, :, :, ind_vt]    
        try: 
            azav_file = get_widest_range_file(dirname + '/data/', 'AZ_Avgs')
            di = np.load(dirname + '/data/' + azav_file).item()
            av_vals, lut = di['vals'], di['lut']
            vr_av = av_vals[:, a.inds, lut[qind_vr]]
            vt_av = av_vals[:, a.inds, lut[qind_vt]]
            vr_prime_slice = vr_slice - vr_av
            vt_prime_slice = vt_slice - vt_av
        except:
            vr_prime_slice = vr_slice - np.mean(vr_slice, axis=0) 
            vt_prime_slice = vt_slice - np.mean(vt_slice, axis=0)             
        sslice = prefactor*(vr_prime_slice**2 + vt_prime_slice**2)        
    elif (varname == 'vlvp'):
        ref = ReferenceState(dirname + '/reference', '')
        prefactor = ref.density[a.inds]        
        sint = (a.sintheta).reshape((1, a.ntheta, 1))
        cost = (a.costheta).reshape((1, a.ntheta, 1))
        ind_vr, ind_vt, ind_vp = a.lut[qind_vr], a.lut[qind_vt], a.lut[qind_vp]
        vr_slice = vals[:, :, :, ind_vr]
        vt_slice = vals[:, :, :, ind_vt]
        vp_slice = vals[:, :, :, ind_vp]      
        try: 
            azav_file = get_widest_range_file(dirname + '/data/', 'AZ_Avgs')
            di = np.load(dirname + '/data/' + azav_file).item()
            av_vals, lut = di['vals'], di['lut']
            vr_av = av_vals[:, a.inds, lut[qind_vr]]
            vt_av = av_vals[:, a.inds, lut[qind_vt]]
            vp_av = av_vals[:, a.inds, lut[qind_vp]]
            vr_prime_slice = vr_slice - vr_av
            vt_prime_slice = vt_slice - vt_av
            vp_prime_slice = vp_slice - vp_av
        except:
            vr_prime_slice = vr_slice - np.mean(vr_slice, axis=0)
            vt_prime_slice = vt_slice - np.mean(vt_slice, axis=0)
            vp_prime_slice = vp_slice - np.mean(vp_slice, axis=0)          
        vl_prime_slice = vr_prime_slice*sint + vt_prime_slice*cost    
        sslice = prefactor*vl_prime_slice*vp_prime_slice
    elif (varname == 'vlvz'):
        ref = ReferenceState(dirname + '/reference', '')
        prefactor = ref.density[a.inds]        
        sint = (a.sintheta).reshape((1, a.ntheta, 1))
        cost = (a.costheta).reshape((1, a.ntheta, 1))
        ind_vr, ind_vt, ind_vp = a.lut[qind_vr], a.lut[qind_vt], a.lut[qind_vp]
        vr_slice = vals[:, :, :, ind_vr]
        vt_slice = vals[:, :, :, ind_vt] 
        vp_slice = vals[:, :, :, ind_vp]       
        try: 
            azav_file = get_widest_range_file(dirname + '/data/', 'AZ_Avgs')
            di = np.load(dirname + '/data/' + azav_file).item()
            av_vals, lut = di['vals'], di['lut']
            vr_av = av_vals[:, a.inds, lut[qind_vr]]
            vt_av = av_vals[:, a.inds, lut[qind_vt]]
            vp_av = av_vals[:, a.inds, lut[qind_vp]]
            vr_prime_slice = vr_slice - vr_av
            vt_prime_slice = vt_slice - vt_av
            vp_prime_slice = vp_slice - vp_av
        except:
            vr_prime_slice = vr_slice - np.mean(vr_slice, axis=0)
            vt_prime_slice = vt_slice - np.mean(vt_slice, axis=0)
            vp_prime_slice = vp_slice - np.mean(vp_slice, axis=0)          
        vl_prime_slice = vr_prime_slice*sint + vt_prime_slice*cost   
        vz_prime_slice = vr_prime_slice*cost - vt_prime_slice*sint  
        sslice = prefactor*vl_prime_slice*vz_prime_slice
    elif (varname == 'vpvz'):
        ref = ReferenceState(dirname + '/reference', '')
        prefactor = ref.density[a.inds]        
        sint = (a.sintheta).reshape((1, a.ntheta, 1))
        cost = (a.costheta).reshape((1, a.ntheta, 1))
        ind_vr, ind_vt, ind_vp = a.lut[qind_vr], a.lut[qind_vt], a.lut[qind_vp]
        vr_slice = vals[:, :, :, ind_vr]
        vt_slice = vals[:, :, :, ind_vt]  
        vp_slice = vals[:, :, :, ind_vp]       
        try: 
            azav_file = get_widest_range_file(dirname + '/data/', 'AZ_Avgs')
            di = np.load(dirname + '/data/' + azav_file).item()
            av_vals, lut = di['vals'], di['lut']
            vr_av = av_vals[:, a.inds, lut[qind_vr]]
            vt_av = av_vals[:, a.inds, lut[qind_vt]]
            vp_av = av_vals[:, a.inds, lut[qind_vp]]
            vr_prime_slice = vr_slice - vr_av
            vt_prime_slice = vt_slice - vt_av
            vp_prime_slice = vp_slice - vp_av
        except:
            vr_prime_slice = vr_slice - np.mean(vr_slice, axis=0)
            vt_prime_slice = vt_slice - np.mean(vt_slice, axis=0)
            vp_prime_slice = vp_slice - np.mean(vp_slice, axis=0)          
        vz_prime_slice = vr_prime_slice*cost - vt_prime_slice*sint  
        sslice = prefactor*vp_prime_slice*vz_prime_slice 

    # Spherical fluctuating magnetic field products 
    # (multiply by 1/4*pi to get units of magnetic pressure)
    elif (varname == 'bsq'):
        prefactor = 1/(4*np.pi)
        ind_br, ind_bt, ind_bp = a.lut[qind_br], a.lut[qind_bt], a.lut[qind_bp]
        br_slice = vals[:, :, :, ind_br]
        bt_slice = vals[:, :, :, ind_bt]
        bp_slice = vals[:, :, :, ind_bp]      
        try: 
            azav_file = get_widest_range_file(dirname + '/data/', 'AZ_Avgs')
            di = np.load(dirname + '/data/' + azav_file).item()
            av_vals, lut = di['vals'], di['lut']
            br_av = av_vals[:, a.inds, lut[qind_br]]
            bt_av = av_vals[:, a.inds, lut[qind_bt]]
            bp_av = av_vals[:, a.inds, lut[qind_bp]]
            br_prime_slice = br_slice - br_av
            bt_prime_slice = bt_slice - bt_av
            bp_prime_slice = bp_slice - bp_av
        except:
            br_prime_slice = br_slice - np.mean(br_slice, axis=0)
            bt_prime_slice = bt_slice - np.mean(bt_slice, axis=0)
            bp_prime_slice = bp_slice - np.mean(bp_slice, axis=0)
        sslice = prefactor*(br_prime_slice**2 + bt_prime_slice**2 + bp_prime_slice**2)
    elif (varname == 'brsq'):
        prefactor = 1/(4*np.pi)
        ind_br = a.lut[qind_br]
        br_slice = vals[:, :, :, ind_br]  
        try: 
            azav_file = get_widest_range_file(dirname + '/data/', 'AZ_Avgs')
            di = np.load(dirname + '/data/' + azav_file).item()
            vals, lut = di['vals'], di['lut']
            br_av = vals[:, a.inds, lut[qind_br]]
            br_prime_slice = br_slice - br_av
        except:
            br_prime_slice = br_slice - np.mean(br_slice, axis=0)
        sslice = prefactor*br_prime_slice**2
    elif (varname == 'btsq'):
        prefactor = 1/(4*np.pi)
        ind_bt = a.lut[qind_bt]
        bt_slice = vals[:, :, :, ind_bt]  
        try: 
            azav_file = get_widest_range_file(dirname + '/data/', 'AZ_Avgs')
            di = np.load(dirname + '/data/' + azav_file).item()
            av_vals, lut = di['vals'], di['lut']
            bt_av = av_vals[:, a.inds, lut[qind_bt]]
            bt_prime_slice = bt_slice - bt_av
        except:
            bt_prime_slice = bt_slice - np.mean(bt_slice, axis=0)
        sslice = prefactor*bt_prime_slice**2
    elif (varname == 'bpsq'):
        prefactor = 1/(4*np.pi)
        ind_bp = a.lut[qind_bp]
        bp_slice = vals[:, :, :, ind_bp]  
        try: 
            azav_file = get_widest_range_file(dirname + '/data/', 'AZ_Avgs')
            di = np.load(dirname + '/data/' + azav_file).item()
            av_vals, lut = di['vals'], di['lut']
            bp_av = av_vals[:, a.inds, lut[qind_bp]]
            bp_prime_slice = bp_slice - bp_av
        except:
            bp_prime_slice = bp_slice - np.mean(bp_slice, axis=0)
        sslice = prefactor*bp_prime_slice**2        
    elif (varname == 'bhsq'):
        prefactor = 1/(4*np.pi)        
        ind_bt, ind_bp = a.lut[qind_bt], a.lut[qind_bp]
        bt_slice = vals[:, :, :, ind_bt]
        bp_slice = vals[:, :, :, ind_bp]      
        try: 
            azav_file = get_widest_range_file(dirname + '/data/', 'AZ_Avgs')
            di = np.load(dirname + '/data/' + azav_file).item()
            av_vals, lut = di['vals'], di['lut']
            bt_av = av_vals[:, a.inds, lut[qind_bt]]
            bp_av = av_vals[:, a.inds, lut[qind_bp]]
            bt_prime_slice = bt_slice - bt_av
            bp_prime_slice = bp_slice - bp_av
        except:
            bt_prime_slice = bt_slice - np.mean(bt_slice, axis=0)
            bp_prime_slice = bp_slice - np.mean(bp_slice, axis=0)
        sslice = prefactor*(bt_prime_slice**2 + bp_prime_slice**2)     
    elif (varname == 'brbp'):
        prefactor = 1/(4*np.pi)        
        ind_br, ind_bp = a.lut[qind_br], a.lut[qind_bp]
        br_slice = vals[:, :, :, ind_br]
        bp_slice = vals[:, :, :, ind_bp]      
        try: 
            azav_file = get_widest_range_file(dirname + '/data/', 'AZ_Avgs')
            di = np.load(dirname + '/data/' + azav_file).item()
            av_vals, lut = di['vals'], di['lut']
            br_av = av_vals[:, a.inds, lut[qind_br]]
            bp_av = av_vals[:, a.inds, lut[qind_bp]]
            br_prime_slice = br_slice - br_av
            bp_prime_slice = bp_slice - bp_av
        except:
            br_prime_slice = br_slice - np.mean(br_slice, axis=0)
            bp_prime_slice = bp_slice - np.mean(bp_slice, axis=0)    
        sslice = prefactor*br_prime_slice*bp_prime_slice      
    elif (varname == 'brbt'):
        prefactor = 1/(4*np.pi)        
        ind_br, ind_bt = a.lut[qind_br], a.lut[qind_bt]
        br_slice = vals[:, :, :, ind_br]
        bt_slice = vals[:, :, :, ind_bt]     
        try: 
            azav_file = get_widest_range_file(dirname + '/data/', 'AZ_Avgs')
            di = np.load(dirname + '/data/' + azav_file).item()
            av_vals, lut = di['vals'], di['lut']
            br_av = av_vals[:, a.inds, lut[qind_br]]
            bt_av = av_vals[:, a.inds, lut[qind_bt]]
            br_prime_slice = br_slice - br_av
            bt_prime_slice = bt_slice - bt_av
        except:
            br_prime_slice = br_slice - np.mean(br_slice, axis=0)
            bt_prime_slice = bt_slice - np.mean(bt_slice, axis=0)
        sslice = prefactor*br_prime_slice*bt_prime_slice
    elif (varname == 'btbp'):
        prefactor = 1/(4*np.pi)        
        ind_bt, ind_bp = a.lut[qind_bt], a.lut[qind_bp]
        bt_slice = vals[:, :, :, ind_bt]
        bp_slice = vals[:, :, :, ind_bp]       
        try: 
            azav_file = get_widest_range_file(dirname + '/data/', 'AZ_Avgs')
            di = np.load(dirname + '/data/' + azav_file).item()
            av_vals, lut = di['vals'], di['lut']
            bt_av = av_vals[:, a.inds, lut[qind_bt]]
            bp_av = av_vals[:, a.inds, lut[qind_bp]]
            bt_prime_slice = bt_slice - bt_av
            bp_prime_slice = bp_slice - bp_av
        except:
            bt_prime_slice = bt_slice - np.mean(bt_slice, axis=0)
            bp_prime_slice = bp_slice - np.mean(bp_slice, axis=0)    
        sslice = prefactor*bt_prime_slice*bp_prime_slice 
        
    # Cylindrical fluctuating magnetic field products velocity products  
    # Multiplied by 1/4*pi to get units of magnetic pressure
    elif (varname == 'blsq'):
        prefactor = 1/(4*np.pi)          
        sint = (a.sintheta).reshape((1, a.ntheta, 1))
        cost = (a.costheta).reshape((1, a.ntheta, 1))        
        ind_br, ind_bt = a.lut[qind_br], a.lut[qind_bt]        
        br_slice = vals[:, :, :, ind_br]
        bt_slice = vals[:, :, :, ind_bt]   
        bl_slice = br_slice*sint + bt_slice*cost
        try: 
            azav_file = get_widest_range_file(dirname + '/data/', 'AZ_Avgs')
            di = np.load(dirname + '/data/' + azav_file).item()
            av_vals, lut = di['vals'], di['lut']
            br_av = av_vals[:, a.inds, lut[qind_br]]
            bt_av = av_vals[:, a.inds, lut[qind_bt]]
            bl_av = br_av*sint + bt_av*cost 
            bl_prime_slice = bl_slice - bl_av
        except:
            bl_prime_slice = bl_slice - np.mean(bl_slice, axis=0) 
        sslice = prefactor*bl_prime_slice**2
    elif (varname == 'bzsq'):
        prefactor = 1/(4*np.pi)          
        sint = (a.sintheta).reshape((1, a.ntheta, 1))
        cost = (a.costheta).reshape((1, a.ntheta, 1))        
        ind_br, ind_bt = a.lut[qind_br], a.lut[qind_bt]        
        br_slice = vals[:, :, :, ind_br]
        bt_slice = vals[:, :, :, ind_bt]    
        bz_slice = br_slice*cost - bt_slice*sint
        try: 
            azav_file = get_widest_range_file(dirname + '/data/', 'AZ_Avgs')
            di = np.load(dirname + '/data/' + azav_file).item()
            av_vals, lut = di['vals'], di['lut']
            br_av = av_vals[:, a.inds, lut[qind_br]]
            bt_av = av_vals[:, a.inds, lut[qind_bt]]
            bz_av = br_av*cost - bt_av*sint 
            bz_prime_slice = bz_slice - bz_av
        except:
            bz_prime_slice = bz_slice - np.mean(bz_slice, axis=0) 
        sslice = prefactor*bz_prime_slice**2 
    elif (varname == 'bpolsq'):
        prefactor = 1/(4*np.pi)                 
        ind_br, ind_bt = a.lut[qind_br], a.lut[qind_bt]        
        br_slice = vals[:, :, :, ind_br]
        bt_slice = vals[:, :, :, ind_bt]    
        try: 
            azav_file = get_widest_range_file(dirname + '/data/', 'AZ_Avgs')
            di = np.load(dirname + '/data/' + azav_file).item()
            av_vals, lut = di['vals'], di['lut']
            br_av = av_vals[:, a.inds, lut[qind_br]]
            bt_av = av_vals[:, a.inds, lut[qind_bt]]
            br_prime_slice = br_slice - br_av
            bt_prime_slice = bt_slice - bt_av
        except:
            br_prime_slice = br_slice - np.mean(br_slice, axis=0) 
            bt_prime_slice = bt_slice - np.mean(bt_slice, axis=0)             
        sslice = prefactor*(br_prime_slice**2 + bt_prime_slice**2)     
    elif (varname == 'blbp'):
        prefactor = 1/(4*np.pi)        
        sint = (a.sintheta).reshape((1, a.ntheta, 1))
        cost = (a.costheta).reshape((1, a.ntheta, 1))
        ind_br, ind_bt, ind_bp = a.lut[qind_br], a.lut[qind_bt], a.lut[qind_bp]
        br_slice = vals[:, :, :, ind_br]
        bt_slice = vals[:, :, :, ind_bt]
        bp_slice = vals[:, :, :, ind_bp]      
        try: 
            azav_file = get_widest_range_file(dirname + '/data/', 'AZ_Avgs')
            di = np.load(dirname + '/data/' + azav_file).item()
            av_vals, lut = di['vals'], di['lut']
            br_av = av_vals[:, a.inds, lut[qind_br]]
            bt_av = av_vals[:, a.inds, lut[qind_bt]]
            bp_av = av_vals[:, a.inds, lut[qind_bp]]
            br_prime_slice = br_slice - br_av
            bt_prime_slice = bt_slice - bt_av
            bp_prime_slice = bp_slice - bp_av
        except:
            br_prime_slice = br_slice - np.mean(br_slice, axis=0)
            bt_prime_slice = bt_slice - np.mean(bt_slice, axis=0)
            bp_prime_slice = bp_slice - np.mean(bp_slice, axis=0)          
        bl_prime_slice = br_prime_slice*sint + bt_prime_slice*cost    
        sslice = prefactor*bl_prime_slice*bp_prime_slice
    elif (varname == 'blbz'):
        prefactor = 1/(4*np.pi)        
        sint = (a.sintheta).reshape((1, a.ntheta, 1))
        cost = (a.costheta).reshape((1, a.ntheta, 1))
        ind_br, ind_bt, ind_bp = a.lut[qind_br], a.lut[qind_bt], a.lut[qind_bp]
        br_slice = vals[:, :, :, ind_br]
        bt_slice = vals[:, :, :, ind_bt] 
        bp_slice = vals[:, :, :, ind_bp]       
        try: 
            azav_file = get_widest_range_file(dirname + '/data/', 'AZ_Avgs')
            di = np.load(dirname + '/data/' + azav_file).item()
            av_vals, lut = di['vals'], di['lut']
            br_av = av_vals[:, a.inds, lut[qind_br]]
            bt_av = av_vals[:, a.inds, lut[qind_bt]]
            bp_av = av_vals[:, a.inds, lut[qind_bp]]
            br_prime_slice = br_slice - br_av
            bt_prime_slice = bt_slice - bt_av
            bp_prime_slice = bp_slice - bp_av
        except:
            br_prime_slice = br_slice - np.mean(br_slice, axis=0)
            bt_prime_slice = bt_slice - np.mean(bt_slice, axis=0)
            bp_prime_slice = bp_slice - np.mean(bp_slice, axis=0)          
        bl_prime_slice = br_prime_slice*sint + bt_prime_slice*cost   
        bz_prime_slice = br_prime_slice*cost - bt_prime_slice*sint  
        sslice = prefactor*bl_prime_slice*bz_prime_slice
    elif (varname == 'bpbz'):
        prefactor = 1/(4*np.pi)        
        sint = (a.sintheta).reshape((1, a.ntheta, 1))
        cost = (a.costheta).reshape((1, a.ntheta, 1))
        ind_br, ind_bt, ind_bp = a.lut[qind_br], a.lut[qind_bt], a.lut[qind_bp]
        br_slice = vals[:, :, :, ind_br]
        bt_slice = vals[:, :, :, ind_bt]  
        bp_slice = vals[:, :, :, ind_bp]       
        try: 
            azav_file = get_widest_range_file(dirname + '/data/', 'AZ_Avgs')
            di = np.load(dirname + '/data/' + azav_file).item()
            av_vals, lut = di['vals'], di['lut']
            br_av = av_vals[:, a.inds, lut[qind_br]]
            bt_av = av_vals[:, a.inds, lut[qind_bt]]
            bp_av = av_vals[:, a.inds, lut[qind_bp]]
            br_prime_slice = br_slice - br_av
            bt_prime_slice = bt_slice - bt_av
            bp_prime_slice = bp_slice - bp_av
        except:
            br_prime_slice = br_slice - np.mean(br_slice, axis=0)
            bt_prime_slice = bt_slice - np.mean(bt_slice, axis=0)
            bp_prime_slice = bp_slice - np.mean(bp_slice, axis=0)          
        bz_prime_slice = br_prime_slice*cost - bt_prime_slice*sint  
        sslice = prefactor*bp_prime_slice*bz_prime_slice 
        
    # Correlations between vertical flow (vr) and temperature/entropy
    elif (varname == 'vrs'):       
        ind_vr, ind_s = a.lut[qind_vr], a.lut[qind_s]
        vr_slice = vals[:, :, :, ind_vr]
        s_slice = vals[:, :, :, ind_s]      
        try: 
            azav_file = get_widest_range_file(dirname + '/data/', 'AZ_Avgs')
            di = np.load(dirname + '/data/' + azav_file).item()
            av_vals, lut = di['vals'], di['lut']
            vr_av = av_vals[:, a.inds, lut[qind_vr]]
            s_av = av_vals[:, a.inds, lut[qind_s]]
            vr_prime_slice = vr_slice - vr_av
            s_prime_slice = s_slice - s_av
        except:
            vr_prime_slice = vr_slice - np.mean(vr_slice, axis=0)
            s_prime_slice = s_slice - np.mean(s_slice, axis=0)    
        sslice = vr_prime_slice*s_prime_slice      
    elif (varname == 'vrt'):       
        ind_vr, ind_s, ind_p = a.lut[qind_vr], a.lut[qind_s], a.lut[qind_p]
        vr_slice = vals[:, :, :, ind_vr]
        s_slice = vals[:, :, :, ind_s] 
        p_slice = vals[:, :, :, ind_p]         
        try: 
            azav_file = get_widest_range_file(dirname + '/data/', 'AZ_Avgs')
            di = np.load(dirname + '/data/' + azav_file).item()
            av_vals, lut = di['vals'], di['lut']
            vr_av = av_vals[:, a.inds, lut[qind_vr]]
            s_av = av_vals[:, a.inds, lut[qind_s]]
            p_av = av_vals[:, a.inds, lut[qind_p]]            
            vr_prime_slice = vr_slice - vr_av
            s_prime_slice = s_slice - s_av
            p_prime_slice = p_slice - p_av            
        except:
            vr_prime_slice = vr_slice - np.mean(vr_slice, axis=0)
            s_prime_slice = s_slice - np.mean(s_slice, axis=0)  
            p_prime_slice = p_slice - np.mean(p_slice, axis=0)               
        ref = ReferenceState(dirname + '/reference', '')
        cp = get_parameter(dirname, 'pressure_specific_heat')
        ref_p = (ref.pressure)[a.inds]
        ref_temp = (ref.temperature)[a.inds]
        try:
            poly_n = get_paramater(dirname, 'poly_n') # must check "poly_n" later
            poly_gamma = 1. + 1./poly_n  
        except:
            poly_gamma = 5./3.
        t_prime_slice = ref_temp*(p_prime_slice/ref_p*(1. - 1./poly_gamma) + s_prime_slice/cp)          
        sslice = vr_prime_slice*t_prime_slice          
    elif (varname == 'vtt'):       
        ind_vt, ind_s, ind_p = a.lut[qind_vt], a.lut[qind_s], a.lut[qind_p]
        vt_slice = vals[:, :, :, ind_vt]
        s_slice = vals[:, :, :, ind_s] 
        p_slice = vals[:, :, :, ind_p]         
        try: 
            azav_file = get_widest_range_file(dirname + '/data/', 'AZ_Avgs')
            di = np.load(dirname + '/data/' + azav_file).item()
            av_vals, lut = di['vals'], di['lut']
            vt_av = av_vals[:, a.inds, lut[qind_vt]]
            s_av = av_vals[:, a.inds, lut[qind_s]]
            p_av = av_vals[:, a.inds, lut[qind_p]]            
            vt_prime_slice = vt_slice - vt_av
            s_prime_slice = s_slice - s_av
            p_prime_slice = p_slice - p_av            
        except:
            vt_prime_slice = vt_slice - np.mean(vt_slice, axis=0)
            s_prime_slice = s_slice - np.mean(s_slice, axis=0)  
            p_prime_slice = p_slice - np.mean(p_slice, axis=0)               
        cp = get_parameter(dirname, 'pressure_specific_heat')
        try:
            ref = ReferenceState(dirname + '/reference', '')
            ref_p = (ref.pressure)[a.inds]
            ref_temp = (ref.temperature)[a.inds]
        except:
            eq = equation_coefficients()
            eq.read(dirname + '/equation_coefficients')
            ref_temp = (eq.temperature)[a.inds]
            ref_dens = (eq.density)[a.inds]
            ref_p = ref_dens*thermo_R*ref_temp
        try:
            poly_n = get_paramater(dirname, 'poly_n') # must check "poly_n" later
            poly_gamma = 1. + 1./poly_n  
        except:
            poly_gamma = 5./3.
        t_prime_slice = ref_temp*(p_prime_slice/ref_p*(1. - 1./poly_gamma) + s_prime_slice/cp)          
        sslice = vt_prime_slice*t_prime_slice          
    return sslice
