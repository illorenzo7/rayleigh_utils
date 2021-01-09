# Author: Loren Matilsky
# Created: 12/19/2018
# Stores properties of various Rayliegh output variables;
# units, quantity codes, LaTex variable names, etc.
# Add to these lists as need be.
# Not all the variables have quantity codes but are derivative quantities of 
# other fluid variables

texlabels = {\
    'vr'        :       r'$v_r$',
    'vt'        :       r'$v_\theta$',
    'vp'        :       r'$v_\phi$',
    'vr_prime'  :       r'$v_r - \langle v_r\rangle$',
    'vt_prime'  :       r'$v_\theta - \langle v_\theta\rangle$',
    'vp_prime'  :       r'$v_\phi - \langle v_\phi\rangle$',
    'vl'        :       r'$v_\lambda$',
    'vz'        :       r'$v_z$',   
    'vl_prime'  :       r'$v_\lambda - \langle v_\lambda\rangle$',
    'vz_prime'  :       r'$v_z - \langle v_z\rangle$',
             
    'br'        :       r'$B_r$',
    'bt'        :       r'$B_\theta$',
    'bp'        :       r'$B_\phi$',
    'br_prime'  :       r'$B_r - \langle B_r\rangle$',
    'bt_prime'  :       r'$B_\theta - \langle B_\theta\rangle$',
    'bp_prime'  :       r'$B_\phi - \langle B_\phi\rangle$',             
    'bl'        :       r'$B_\lambda$',
    'bz'        :       r'$B_z$',   
    'bl_prime'  :       r'$B_\lambda - \langle B_\lambda\rangle$',
    'bz_prime'  :       r'$B_z - \langle B_z\rangle$',
             
    'omr'       :       r'$\omega_r$',
    'omt'       :       r'$\omega_\theta$',
    'omp'       :       r'$\omega_\phi$',
    'omr_prime' :       r'$\omega_r - \langle \omega_r\rangle$',  
    'omt_prime' :       r'$\omega_\theta - \langle \omega_\theta\rangle$', 
    'omp_prime' :       r'$\omega_\phi - \langle \omega_\phi\rangle$',
    'oml'       :       r'$\omega_\lambda$',
    'omz'       :       r'$\omega_z$',
    'oml_prime' :       r'$\omega_\lambda - \langle \omega_\lambda\rangle$',
    'omz_prime' :       r'$\omega_z - \langle \omega_z\rangle$',
             
    'vsq'       :       r'$\overline{\rho}(v^\prime)^2$',
    'vrsq'      :       r'$\overline{\rho}(v^\prime_r)^2$',
    'vtsq'      :       r'$\overline{\rho}(v^\prime_\theta)^2$',
    'vpsq'      :       r'$\overline{\rho}(v^\prime_\phi)^2$',             
    'vhsq'      :       r'$\overline{\rho}[(v^\prime_\theta)^2\ +\ (v^\prime_\phi)^2]$',
    'vlsq'      :       r'$\overline{\rho}(v^\prime_\lambda)^2$',    
    'vzsq'      :       r'$\overline{\rho}(v^\prime_z)^2$', 
    'vmsq'    :       r'$\overline{\rho}[(v^\prime_r)^2 + (v^\prime_\theta)^2]$', 
             
    'omsq'      :       r'$\omega^2$',
    'omsq_prime' :       r'$(\omega^\prime)^2$',
             
    'bsq'       :       r'$(1/4\pi)(B^\prime)^2$',
    'brsq'      :       r'$(1/4\pi)(B^\prime_r)^2$',
    'btsq'      :       r'$(1/4\pi)(B^\prime_\theta)^2$',
    'bpsq'      :       r'$(1/4\pi)(B^\prime_\phi)^2$',               
    'bhsq'      :       r'$(1/4\pi)[(B^\prime_\theta)^2\ +\ (B^\prime_\phi)^2]$',
             
    'blsq'      :       r'$(1/4\pi)(B^\prime_\lambda)^2$',    
    'bzsq'      :       r'$(1/4\pi)(B^\prime_z)^2$', 
    'bmsq'    :       r'$(1/4\pi)[(B^\prime_r)^2 + (B^\prime_\theta)^2]$',             
             
    's'         :       r'$S$',
    'p'         :       r'$P$',  
    'rho'       :       r'$\rho$', 
    't'         :       r'$T$',  
             
    's_prime'   :       r'$S - \langle S \rangle$',
    'p_prime'   :       r'$P - \langle P \rangle$',             
    'rho_prime' :       r'$\rho - \langle \rho \rangle$',
    't_prime'   :       r'$T - \langle T \rangle$',  
             
    's_prime_sph'   :       r'$S - \langle S \rangle_s$',
    'p_prime_sph'   :       r'$P - \langle P \rangle_s$',             
    'rho_prime_sph' :       r'$\rho - \langle \rho \rangle_s$',
    't_prime_sph'   :       r'$T - \langle T \rangle_s$',             
             
    'vrt_sph'   :       r'$\overline{\rho}c_p(v_r - \langle v_r\rangle_{\rm{sph}})(T - \langle T\rangle_{\rm{sph}})$',
    'vrs_sph'   :       r'$\overline{\rho}\overline{T}(v_r - \langle v_r\rangle_{\rm{sph}})(S - \langle S\rangle_{\rm{sph}})$',
    'vrp_sph'   :       r'$(v_r - \langle v_r\rangle_{\rm{sph}})(P - \langle P\rangle_{\rm{sph}})$',
    'vtt'       :       r'$v_\theta^\prime T^\prime$',
             
    'vrvp'      :       r'$\overline{\rho}v_r^\prime v_\phi^\prime$',
    'vrvt'      :       r'$\overline{\rho}v_r^\prime v_\theta^\prime$',
    'vtvp'      :       r'$\overline{\rho}v_\theta^\prime v_\phi^\prime$',
    'vlvp'      :       r'$\overline{\rho}v_\lambda^\prime v_\phi^\prime$',
    'vpvz'      :       r'$\overline{\rho}v_\phi^\prime v_z^\prime$',
    'vlvz'      :       r'$\overline{\rho}v_\lambda^\prime v_z^\prime$',
             
    'brbp'      :       r'$(1/4\pi)B_r^\prime B_\phi^\prime$',
    'brbt'      :       r'$(1/4\pi)B_r^\prime B_\theta^\prime$',
    'btbp'      :       r'$(1/4\pi)B_\theta^\prime B_\phi^\prime$',
    'blbp'      :       r'$(1/4\pi)B_\lambda^\prime B_\phi^\prime$',
    'bpbz'      :       r'$(1/4\pi)B_\phi^\prime B_z^\prime$',
    'blbz'      :       r'$(1/4\pi)B_\lambda^\prime B_z^\prime$',
    'dvtdr'     :       r'$\partial v_\theta^\prime/\partial r$',
    'brdvtdr'   :       r'$B_r^\prime \partial v_\theta^\prime/\partial r$'
    }

texunits = {\
    'vr'        :       r'$\rm{m}\ \rm{s}^{-1}$',
    'vt'        :       r'$\rm{m}\ \rm{s}^{-1}$',
    'vp'        :       r'$\rm{m}\ \rm{s}^{-1}$',
    'vr_prime'  :       r'$\rm{m}\ \rm{s}^{-1}$',
    'vt_prime'  :       r'$\rm{m}\ \rm{s}^{-1}$',           
    'vp_prime'  :       r'$\rm{m}\ \rm{s}^{-1}$',
    'vl'        :       r'$\rm{m}\ \rm{s}^{-1}$',
    'vz'        :       r'$\rm{m}\ \rm{s}^{-1}$',
    'vl_prime'  :       r'$\rm{m}\ \rm{s}^{-1}$',
    'vz_prime'  :       r'$\rm{m}\ \rm{s}^{-1}$',
            
    'br'        :       'G',
    'bt'        :       'G',
    'bp'        :       'G',
    'br_prime'  :       'G',   
    'bt_prime'  :       'G',               
    'bp_prime'  :       'G',            
    'bl'        :       'G',
    'bz'        :       'G',
    'bl_prime'  :       'G',               
    'bz_prime'  :       'G',   
            
    'omr'       :       r'$\rm{rad}\ \rm{s}^{-1}$',
    'omt'       :       r'$\rm{rad}\ \rm{s}^{-1}$',
    'omp'       :       r'$\rm{rad}\ \rm{s}^{-1}$',
    'omr_prime' :       r'$\rm{rad}\ \rm{s}^{-1}$',  
    'omt_prime' :       r'$\rm{rad}\ \rm{s}^{-1}$',  
    'omp_prime' :       r'$\rm{rad}\ \rm{s}^{-1}$',              
    'oml'       :       r'$\rm{rad}\ \rm{s}^{-1}$',
    'omz'       :       r'$\rm{rad}\ \rm{s}^{-1}$',
    'oml_prime' :       r'$\rm{rad}\ \rm{s}^{-1}$',              
    'omz_prime' :       r'$\rm{rad}\ \rm{s}^{-1}$',  
            
    'vsq'       :       r'$\rm{dyn}\ \rm{cm}^{-2}$',
    'vrsq'      :       r'$\rm{dyn}\ \rm{cm}^{-2}$',
    'vtsq'      :       r'$\rm{dyn}\ \rm{cm}^{-2}$',
    'vpsq'      :       r'$\rm{dyn}\ \rm{cm}^{-2}$',            
    'vhsq'      :       r'$\rm{dyn}\ \rm{cm}^{-2}$',
             
    'vlsq'      :       r'$\rm{dyn}\ \rm{cm}^{-2}$',    
    'vzsq'      :       r'$\rm{dyn}\ \rm{cm}^{-2}$', 
    'vmsq'    :       r'$\rm{dyn}\ \rm{cm}^{-2}$', 
            
    'omsq'      :       r'$\rm{rad}^2\ \rm{s}^{-2}$',
    'omsq_prime' :       r'$\rm{rad}^2\ \rm{s}^{-2}$',            

    'bsq'       :       r'$\rm{dyn}\ \rm{cm}^{-2}$',
    'brsq'      :       r'$\rm{dyn}\ \rm{cm}^{-2}$',
    'btsq'      :       r'$\rm{dyn}\ \rm{cm}^{-2}$',
    'bpsq'      :       r'$\rm{dyn}\ \rm{cm}^{-2}$',            
    'bhsq'      :       r'$\rm{dyn}\ \rm{cm}^{-2}$',
             
    'blsq'      :       r'$\rm{dyn}\ \rm{cm}^{-2}$',    
    'bzsq'      :       r'$\rm{dyn}\ \rm{cm}^{-2}$', 
    'bmsq'    :       r'$\rm{dyn}\ \rm{cm}^{-2}$',    
            
            
    's'         :       r'$\rm{erg}\ \rm{K}^{-1}\ \rm{g}^{-1}$',
    's_prime'   :       r'$\rm{erg}\ \rm{K}^{-1}\ \rm{g}^{-1}$',
    's_prime_sph'   :       r'$\rm{erg}\ \rm{K}^{-1}\ \rm{g}^{-1}$',            
    'p'         :       r'$\rm{dyn}\ \rm{cm}^{-2}$',
    'p_prime'   :       r'$\rm{dyn}\ \rm{cm}^{-2}$',            
    'p_prime_sph'   :       r'$\rm{dyn}\ \rm{cm}^{-2}$',  
    'rho'       :       r'$\rm{g}\ \rm{cm}^{-3}$',
    'rho_prime' :       r'$\rm{g}\ \rm{cm}^{-3}$',            
    'rho_prime_sph' :   r'$\rm{g}\ \rm{cm}^{-3}$', 
    't'       :         r'$\rm{K}$',
    't_prime' :         r'$\rm{K}$',          
    't_prime_sph' :     r'$\rm{K}$',            
    
    'vrt_sph'   :       r'$\rm{erg\ cm^{-2}\ s^{-1}}$', 
    'vrs_sph'   :       r'$\rm{erg\ cm^{-2}\ s^{-1}}$', 
    'vrp_sph'   :       r'$\rm{erg\ cm^{-2}\ s^{-1}}$', 
    'vtt'       :       r'$\rm{cm\ s^{-1}\ K}$', 
           
    'vrvp'      :       r'$\rm{dyn}\ \rm{cm}^{-2}$', 
    'vrvt'      :       r'$\rm{dyn}\ \rm{cm}^{-2}$', 
    'vtvp'      :       r'$\rm{dyn}\ \rm{cm}^{-2}$', 
    
    'vlvp'      :       r'$\rm{dyn}\ \rm{cm}^{-2}$', 
    'vpvz'      :       r'$\rm{dyn}\ \rm{cm}^{-2}$', 
    'vlvz'      :       r'$\rm{dyn}\ \rm{cm}^{-2}$', 
           
    'brbp'      :       r'$\rm{dyn}\ \rm{cm}^{-2}$', 
    'brbt'      :       r'$\rm{dyn}\ \rm{cm}^{-2}$', 
    'btbp'      :       r'$\rm{dyn}\ \rm{cm}^{-2}$', 
    
    'blbp'      :       r'$\rm{dyn}\ \rm{cm}^{-2}$', 
    'bpbz'      :       r'$\rm{dyn}\ \rm{cm}^{-2}$', 
    'blbz'      :       r'$\rm{dyn}\ \rm{cm}^{-2}$',
    'brdvtdr'   :       r'$\rm{G\ cm^{-1}}$',
    'dvtdr'     :       r'$\rm{cm^{-1}}$' }

var_indices = {\
    'vr'    :       1, 
    'vt'    :       2,
    'vp'    :       3,              
    'omr'   :       301,
    'omt'   :       302,
    'omp'   :       303,
    's'     :       501,
    'p'     :       502,
    'br'    :       801,
    'bt'    :       802,
    'bp'    :       803,
    'dvtdr' :       11}

var_indices_old = {\
    'vr'    :       1, 
    'vt'    :       2,
    'vp'    :       3,
    's'     :       64,
    'p'     :       65,              
    'omr'   :       103,
    'omt'   :       104,
    'omp'   :       105,
    'br'    :       301,  # not sure if these are right, must check (12/24/2018)
    'bt'    :       302,
    'bp'    :       303,}
