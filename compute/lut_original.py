"""
Module to hold the lookup table for Rayleigh outputs.
There is support for the current lookup table as well
as the older version, the "alternate" lookup table.

Python dictionaries are used where the key is the
index and the value is the Fortran variable name
taken from the source code

For example the radial velocity and the radial
component of {B dot grad v} are stored as:

    lut[1] = 'v_r'
    lut[506] = 'induction_shear_r'

and for the alternate lookup table, total kinetic
energy and entropy/temperature would be:

    alt_lut[4] = 'temperature'
    alt_lut[7] = 'kinetic_energy'

R. Orvedahl 7-27-2016
"""

from collections import OrderedDict

lut = OrderedDict()      # current lookup table
lut_tex = OrderedDict()
alt_lut = OrderedDict()  # older lookup table

lut[-10]     = 'version 0'
alt_lut[-10] = 'version 0'

##################################
# current values as of 7/15/16
#   index/values are from the source code
##################################
#----------------
# velocity
#----------------
voffset = 0
# field components
lut[voffset+1] = 'v_r';  lut[voffset+2] = 'v_theta';  lut[voffset+3] = 'v_phi'   # full
lut[voffset+4] = 'vp_r'; lut[voffset+5] = 'vp_theta'; lut[voffset+6] = 'vp_phi'  # fluctuating
lut[voffset+7] = 'vm_r'; lut[voffset+8] = 'vm_theta'; lut[voffset+9] = 'vm_phi'  # mean

# radial derivatives
lut[voffset+10] = 'dv_r_dr';  lut[voffset+11] = 'dv_theta_dr';  lut[voffset+12] = 'dv_phi_dr'   # full
lut[voffset+13] = 'dvp_r_dr'; lut[voffset+14] = 'dvp_theta_dr'; lut[voffset+15] = 'dvp_phi_dr'  # fluctuating
lut[voffset+16] = 'dvm_r_dr'; lut[voffset+17] = 'dvm_theta_dr'; lut[voffset+18] = 'dvm_phi_dr'  # mean

# theta derivatives
lut[voffset+19] = 'dv_r_dt';  lut[voffset+20] = 'dv_theta_dt';  lut[voffset+21] = 'dv_phi_dt'   # full
lut[voffset+22] = 'dvp_r_dt'; lut[voffset+23] = 'dvp_theta_dt'; lut[voffset+24] = 'dvp_phi_dt'  # fluctuating
lut[voffset+25] = 'dvm_r_dt'; lut[voffset+26] = 'dvm_theta_dt'; lut[voffset+27] = 'dvm_phi_dt'  # mean

# phi derivatives
lut[voffset+28] = 'dv_r_dp';  lut[voffset+29] = 'dv_theta_dp';  lut[voffset+30] = 'dv_phi_dp'   # full
lut[voffset+31] = 'dvp_r_dp'; lut[voffset+32] = 'dvp_theta_dp'; lut[voffset+33] = 'dvp_phi_dp'  # fluctuating
lut[voffset+34] = 'dvm_r_dp'; lut[voffset+35] = 'dvm_theta_dp'; lut[voffset+36] = 'dvm_phi_dp'  # mean

# (1/r)*theta derivatives
lut[voffset+37] = 'dv_r_dtr';  lut[voffset+38] = 'dv_theta_dtr';  lut[voffset+39] = 'dv_phi_dtr'   # full
lut[voffset+40] = 'dvp_r_dtr'; lut[voffset+41] = 'dvp_theta_dtr'; lut[voffset+42] = 'dvp_phi_dtr'  # fluctuating
lut[voffset+43] = 'dvm_r_dtr'; lut[voffset+44] = 'dvm_theta_dtr'; lut[voffset+45] = 'dvm_phi_dtr'  # mean

# (1/{r sintheta})*phi derivatives
lut[voffset+46] = 'dv_r_dprs';  lut[voffset+47] = 'dv_theta_dprs';  lut[voffset+48] = 'dv_phi_dprs'   # full
lut[voffset+49] = 'dvp_r_dprs'; lut[voffset+50] = 'dvp_theta_dprs'; lut[voffset+51] = 'dvp_phi_dprs'  # fluctuating
lut[voffset+52] = 'dvm_r_dprs'; lut[voffset+53] = 'dvm_theta_dprs'; lut[voffset+54] = 'dvm_phi_dprs'  # mean

#----------------
# mass flux
#----------------
lut[voffset+55] = 'rhov_r';  lut[voffset+56] = 'rhov_theta';  lut[voffset+57] = 'rhov_phi'   # full
lut[voffset+58] = 'rhovp_r'; lut[voffset+59] = 'rhovp_theta'; lut[voffset+60] = 'rhovp_phi'  # fluctuating
lut[voffset+61] = 'rhovm_r'; lut[voffset+62] = 'rhovm_theta'; lut[voffset+63] = 'rhovm_phi'  # mean

#----------------
# thermodynamics
#----------------
pt_off = voffset+63 # pt_off = 63
lut[pt_off+1] = 'entropy';   lut[pt_off+2] = 'pressure'    # full
lut[pt_off+3] = 'entropy_p'; lut[pt_off+4] = 'pressure_p'  # fluctuating
lut[pt_off+5] = 'entropy_m'; lut[pt_off+6] = 'pressure_m'  # mean

# radial derivs
lut[pt_off+7]  = 'entropy_dr';   lut[pt_off+8]  = 'pressure_dr'    # full
lut[pt_off+9]  = 'entropy_p_dr'; lut[pt_off+10] = 'pressure_p_dr'  # fluctuating
lut[pt_off+11] = 'entropy_m_dr'; lut[pt_off+12] = 'pressure_m_dr'  # mean

# theta derivs
lut[pt_off+13] = 'entropy_dtheta';   lut[pt_off+14] = 'pressure_dtheta'    # full
lut[pt_off+15] = 'entropy_p_dtheta'; lut[pt_off+16] = 'pressure_p_dtheta'  # fluctuating
lut[pt_off+17] = 'entropy_m_dtheta'; lut[pt_off+18] = 'pressure_m_dtheta'  # mean

# phi derivs
lut[pt_off+19] = 'entropy_dphi';   lut[pt_off+20] = 'pressure_dphi'    # full
lut[pt_off+21] = 'entropy_p_dphi'; lut[pt_off+22] = 'pressure_p_dphi'  # fluctuating
lut[pt_off+23] = 'entropy_m_dphi'; lut[pt_off+24] = 'pressure_m_dphi'  # mean

# (1/r)*theta derivs
lut[pt_off+25] = 'entropy_dtr';   lut[pt_off+26] = 'pressure_dtr'    # full
lut[pt_off+27] = 'entropy_p_dtr'; lut[pt_off+28] = 'pressure_p_dtr'  # fluctuating
lut[pt_off+29] = 'entropy_m_dtr'; lut[pt_off+30] = 'pressure_m_dtr'  # mean

# (1/{r sintheta})*phi derivs
lut[pt_off+31] = 'entropy_dprs';   lut[pt_off+32] = 'pressure_dprs'    # full
lut[pt_off+33] = 'entropy_p_dprs'; lut[pt_off+34] = 'pressure_p_dprs'  # fluctuating
lut[pt_off+35] = 'entropy_m_dprs'; lut[pt_off+36] = 'pressure_m_dprs'  # mean

# rho_bar * d_by_dr(P/rho_bar)
lut[pt_off+37] = 'rhopressure_dr'   # full
lut[pt_off+38] = 'rhopressurep_dr'  # fluc
lut[pt_off+39] = 'rhopressurem_dr'  # mean

#----------------
# vorticity
#----------------
vort_off = pt_off + 39 # vort_off = 102
lut[vort_off+1] = 'vort_r';  lut[vort_off+2] = 'vort_theta';  lut[vort_off+3] = 'vort_phi'   # full
lut[vort_off+4] = 'vortp_r'; lut[vort_off+5] = 'vortp_theta'; lut[vort_off+6] = 'vortp_phi'  # fluctuating
lut[vort_off+7] = 'vortm_r'; lut[vort_off+8] = 'vortm_theta'; lut[vort_off+9] = 'vortm_phi'  # mean

# enstrophy
lut[vort_off+10] = 'enstrophy'     # enstrophy
lut[vort_off+11] = 'enstrophy_pm'  # (fluctuating-mean)
lut[vort_off+12] = 'enstrophy_mm'  # (mean-mean)
lut[vort_off+13] = 'enstrophy_pp'  # (fluct-fluct)

#----------------
# fluxes
#----------------
eoffset = vort_off + 13 # eoffest = 115
lut[eoffset+1] = 'ecrossb_r'             # [ExB]_r (un-normalized Poynting flux)
lut[eoffset+2] = 'ke_flux_radial'        # vr*KE
lut[eoffset+3] = 'thermale_flux_radial'  # thermal energy radial flux: vr*rho_bar*T_bar*S OR vr*T
lut[eoffset+4] = 'enth_flux_radial'      # vr*cp*rho_bar*T
lut[eoffset+5] = 'visc_flux_r'           # -[Div dot D]_r
lut[eoffset+6] = 'vol_heat_flux'         # "Flux" associated with volumetric heating
lut[eoffset+7] = 'cond_flux_r'           # thermal conductive flux
lut[eoffset+8] = 'vol_heating'           # volumetric heating function
lut[eoffset+9] = 'visc_heating'          # viscous heating

#----------------
# KE
#----------------
keoffset = eoffset + 9 # keoffest = 124
lut[keoffset+1] = 'kinetic_energy'  # 1/2 rho_bar v^2
lut[keoffset+2] = 'radial_ke'       # 1/2 rho_bar {v_r}^2
lut[keoffset+3] = 'theta_ke'        # 1/2 rho_bar {v_theta}^2
lut[keoffset+4] = 'phi_ke'          # 1/2 rho_bar {v_phi}^2

lut[keoffset+5] = 'mkinetic_energy'  # 1/2 rho_bar <v>^2
lut[keoffset+6] = 'radial_mke'       # 1/2 rho_bar <v_r>^2
lut[keoffset+7] = 'theta_mke'        # 1/2 rho_bar <v_theta>^2
lut[keoffset+8] = 'phi_mke'          # 1/2 rho_bar <v_phi>^2

lut[keoffset+9] = 'pkinetic_energy'  # 1/2 rho_bar {v'}^2
lut[keoffset+10] = 'radial_pke'      # 1/2 rho_bar {v_r'}^2
lut[keoffset+11] = 'theta_pke'       # 1/2 rho_bar {v_theta'}^2
lut[keoffset+12] = 'phi_pke'         # 1/2 rho_bar {v_phi'}^2

lut[keoffset+13] = 'vsq'         # v^2
lut[keoffset+14] = 'radial_vsq'  # {v_r}^2
lut[keoffset+15] = 'theta_vsq'   # {v_theta}^2
lut[keoffset+16] = 'phi_vsq'     # {v_phi}^2

lut[keoffset+17] = 'mvsq'         # <v^2
lut[keoffset+18] = 'radial_mvsq'  # <v_r>^2
lut[keoffset+19] = 'theta_mvsq'   # <v_theta>^2
lut[keoffset+20] = 'phi_mvsq'     # <v_phi>^2

lut[keoffset+21] = 'pvsq'         # {v'}^2
lut[keoffset+22] = 'radial_pvsq'  # {v_r'}^2
lut[keoffset+23] = 'theta_pvsq'   # {v_theta'}^2
lut[keoffset+24] = 'phi_pvsq'     # {v_phi'}^2

#----------------
# angular momentum
#----------------
amoff = keoffset + 24 # amoff = 148
lut[amoff+1] = 'amom_fluct_r'      # rho_bar * r *sintheta * {v_r'} * {v_phi'}
lut[amoff+2] = 'amom_fluct_theta'  # rho_bar * r *sintheta * {v_theta'} * {v_phi'}
lut[amoff+3] = 'amom_dr_r'         # rho_bar * r *sintheta * <v_r> * <v_phi>}
lut[amoff+4] = 'amom_dr_theta'     # rho_bar * r *sintheta * <v_theta> * <v_phi>}
lut[amoff+5] = 'amom_mean_r'       # rho_bar * r^2 *sintheta^2 * <v_r> * Omega
lut[amoff+6] = 'amom_mean_theta'   # rho_bar * r^2 *sintheta^2 * <v_theta> * Omega

#----------------
# advection
#----------------
vgv = amoff + 6 # vgv = 154
lut[vgv+1]  = 'v_grad_v_r';   lut[vgv+2]  = 'v_grad_v_theta';   lut[vgv+3]  = 'v_grad_v_phi'    # v dot grad v
lut[vgv+4]  = 'vp_grad_vm_r'; lut[vgv+5]  = 'vp_grad_vm_theta'; lut[vgv+6]  = 'vp_grad_vm_phi'  # v' dot grad <v>
lut[vgv+7]  = 'vm_grad_vp_r'; lut[vgv+8]  = 'vm_grad_vp_theta'; lut[vgv+9]  = 'vm_grad_vp_phi'  # <v> dot grad v'
lut[vgv+10] = 'vp_grad_vp_r'; lut[vgv+11] = 'vp_grad_vp_theta'; lut[vgv+12] = 'vp_grad_vp_phi'  # v' dot grad v'
lut[vgv+13] = 'vm_grad_vm_r'; lut[vgv+14] = 'vm_grad_vm_theta'; lut[vgv+15] = 'vm_grad_vm_phi'  # <v> dot grad <v>

#----------------
# linear forces
#----------------
force_offset = vgv + 15 # force_offset = 169
lut[force_offset+1] = 'buoyancy_force';    lut[force_offset+2] = 'buoyancy_pforce';       lut[force_offset+3] = 'buoyancy_mforce'
lut[force_offset+4] = 'coriolis_force_r';  lut[force_offset+5] = 'coriolis_force_theta';  lut[force_offset+6] = 'coriolis_force_phi'
lut[force_offset+7] = 'coriolis_pforce_r'; lut[force_offset+8] = 'coriolis_pforce_theta'; lut[force_offset+9] = 'coriolis_pforce_phi'
lut[force_offset+10]= 'coriolis_mforce_r'; lut[force_offset+11]= 'coriolis_mforce_theta'; lut[force_offset+12]= 'coriolis_mforce_phi'

# viscous forces
lut[force_offset+13]= 'viscous_force_r';   lut[force_offset+14]= 'viscous_force_theta';   lut[force_offset+15]= 'viscous_force_phi'
lut[force_offset+16]= 'viscous_pforce_r';  lut[force_offset+17]= 'viscous_pforce_theta';  lut[force_offset+18]= 'viscous_pforce_phi'
lut[force_offset+19]= 'viscous_mforce_r';  lut[force_offset+20]= 'viscous_mforce_theta';  lut[force_offset+21]= 'viscous_mforce_phi'

# checks
dcheck = force_offset + 21 # = 190
lut[dcheck+1] = 'diagnostic1'
lut[dcheck+2] = 'diagnostic2'

#----------------
# custom hydro
#----------------
custom_hydro = 300
lut[custom_hydro+1] = 'custom-number-1'

#----------------
# B fields
#----------------
boffset = 400
# field components
lut[boffset+1] = 'b_r';  lut[boffset+2] = 'b_theta';  lut[boffset+3] = 'b_phi'   # full
lut[boffset+4] = 'bp_r'; lut[boffset+5] = 'bp_theta'; lut[boffset+6] = 'bp_phi'  # fluctuating
lut[boffset+7] = 'bm_r'; lut[boffset+8] = 'bm_theta'; lut[boffset+9] = 'bm_phi'  # mean

# radial derivatives
lut[boffset+10] = 'db_r_dr';  lut[boffset+11] = 'db_theta_dr';  lut[boffset+12] = 'db_phi_dr'   # full
lut[boffset+13] = 'dbp_r_dr'; lut[boffset+14] = 'dbp_theta_dr'; lut[boffset+15] = 'dbp_phi_dr'  # fluctuating
lut[boffset+16] = 'dbm_r_dr'; lut[boffset+17] = 'dbm_theta_dr'; lut[boffset+18] = 'dbm_phi_dr'  # mean

# theta derivatives
lut[boffset+19] = 'db_r_dt';  lut[boffset+20] = 'db_theta_dt';  lut[boffset+21] = 'db_phi_dt'   # full
lut[boffset+22] = 'dbp_r_dt'; lut[boffset+23] = 'dbp_theta_dt'; lut[boffset+24] = 'dbp_phi_dt'  # fluctuating
lut[boffset+25] = 'dbm_r_dt'; lut[boffset+26] = 'dbm_theta_dt'; lut[boffset+27] = 'dbm_phi_dt'  # mean

# phi derivatives
lut[boffset+28] = 'db_r_dp';  lut[boffset+29] = 'db_theta_dp';  lut[boffset+30] = 'db_phi_dp'   # full
lut[boffset+31] = 'dbp_r_dp'; lut[boffset+32] = 'dbp_theta_dp'; lut[boffset+33] = 'dbp_phi_dp'  # fluctuating
lut[boffset+34] = 'dbm_r_dp'; lut[boffset+35] = 'dbm_theta_dp'; lut[boffset+36] = 'dbm_phi_dp'  # mean

# (1/r)*theta derivatives
lut[boffset+37] = 'db_r_dtr';  lut[boffset+38] = 'db_theta_dtr';  lut[boffset+39] = 'db_phi_dtr'   # full
lut[boffset+40] = 'dbp_r_dtr'; lut[boffset+41] = 'dbp_theta_dtr'; lut[boffset+42] = 'dbp_phi_dtr'  # fluctuating
lut[boffset+43] = 'dbm_r_dtr'; lut[boffset+44] = 'dbm_theta_dtr'; lut[boffset+45] = 'dbm_phi_dtr'  # mean

# (1/{r sintheta})*phi derivatives
lut[boffset+46] = 'db_r_dprs';  lut[boffset+47] = 'db_theta_dprs';  lut[boffset+48] = 'db_phi_dprs'   # full
lut[boffset+49] = 'dbp_r_dprs'; lut[boffset+50] = 'dbp_theta_dprs'; lut[boffset+51] = 'dbp_phi_dprs'  # fluctuating
lut[boffset+52] = 'dbm_r_dprs'; lut[boffset+53] = 'dbm_theta_dprs'; lut[boffset+54] = 'dbm_phi_dprs'  # mean

#----------------
# current density
#----------------
joffset = boffset + 54 # = 454
lut[joffset+1] = 'j_r';     lut[joffset+2] = 'jp_r';     lut[joffset+3] = 'jm_r'
lut[joffset+4] = 'j_theta'; lut[joffset+5] = 'jp_theta'; lut[joffset+6] = 'jm_theta'
lut[joffset+7] = 'j_phi';   lut[joffset+8] = 'jp_phi';   lut[joffset+9] = 'jm_phi'
lut[joffset+10] = 'j_r_sq'           # (j_r)^2
lut[joffset+11] = 'jp_r_sq'          # (jp_r)^2
lut[joffset+12] = 'j_theta_sq'       # (j_theat)^2
lut[joffset+13] = 'jp_theta_sq'      # (jp_theat)^2
lut[joffset+14] = 'j_phi_sq'         # (j_phi)^2
lut[joffset+15] = 'jp_phi_sq'        # (jp_phi)^2
lut[joffset+16] = 'j_sq'             # j dot j
lut[joffset+17] = 'jp_sq'            # j' dot j'
lut[joffset+18] = 'ohmic_heat'       # eta{  j dot j }
lut[joffset+19] = 'ohmic_heat_pp'    # eta{  j' dot j'}
lut[joffset+20] = 'ohmic_heat_mm'    # eta{ <j> dot <j>}

#----------------
# magnetic energies
#----------------
meoffset = joffset + 20 # = 474
lut[meoffset+1] = 'magnetic_energy'  # B^2
lut[meoffset+2] = 'radial_me'        # {B_r}^2
lut[meoffset+3] = 'theta_me'         # {B_theta}^2
lut[meoffset+4] = 'phi_me'           # {B_phi}^2

lut[meoffset+5] = 'mmagnetic_energy' # <B>^2
lut[meoffset+6] = 'radial_mme'       # <B_r>^2
lut[meoffset+7] = 'theta_mme'        # <B_theta>^2
lut[meoffset+8] = 'phi_mme'          # <B_phi>^2

lut[meoffset+9] = 'pmagnetic_energy' # {B'}^2
lut[meoffset+10] = 'radial_pme'      # {B_r'}^2
lut[meoffset+11] = 'theta_pme'       # {B_theta'}^2
lut[meoffset+12] = 'phi_pme'         # {B_phi'}^2

#----------------
# lorentz forces
#----------------
loff = meoffset + 12 # = 486
lut[loff+1]  = 'j_cross_b_r';   lut[loff+2]  = 'j_cross_b_theta';   lut[loff+3]  = 'j_cross_b_phi'
lut[loff+4]  = 'jp_cross_bm_r'; lut[loff+5]  = 'jp_cross_bm_theta'; lut[loff+6]  = 'jp_cross_bm_phi'
lut[loff+7]  = 'jm_cross_bp_r'; lut[loff+8]  = 'jm_cross_bp_theta'; lut[loff+9]  = 'jm_cross_bp_phi'
lut[loff+10] = 'jm_cross_bm_r'; lut[loff+11] = 'jm_cross_bm_theta'; lut[loff+12] = 'jm_cross_bm_phi'
lut[loff+13] = 'jp_cross_bp_r'; lut[loff+14] = 'jp_cross_bp_theta'; lut[loff+15] = 'jp_cross_bp_phi'

# maxwell stresses
lut[loff+16] = 'maxwell_stress_r'      # -rsinheta {B_r'}{B_phi'}*Lorentz_Coeff
lut[loff+17] = 'maxwell_stress_theta'  # -rsinheta {B_theta'}{B_phi'}*Lorentz_Coeff
lut[loff+18] = 'maxwell_torque_r'      # -rsintheta <B_r><B_phi>*Lorentz_Coeff
lut[loff+19] = 'maxwell_torque_theta'  # -rsintheta <B_theta><B_phi>*Lorentz_Coeff

#----------------
# induction terms
#----------------
indoff = loff + 19 # = 505
lut[indoff+1]  = 'induction_shear_r'      # {B dot grad v}
lut[indoff+2]  = 'induction_comp_r'       # -{div dot v}B
lut[indoff+3]  = 'induction_advec_r'      # -{v dot grad B}
lut[indoff+4]  = 'induction_r'            # del cross {v x B}
lut[indoff+5]  = 'induction_diff_r'       # -del x (eta {del x B})
lut[indoff+6]  = 'induction_shear_theta'  # theta components
lut[indoff+7]  = 'induction_comp_theta'
lut[indoff+8]  = 'induction_advec_theta'
lut[indoff+9]  = 'induction_theta'
lut[indoff+10] = 'induction_diff_theta'
lut[indoff+11] = 'induction_shear_phi'    # phi components
lut[indoff+12] = 'induction_comp_phi'
lut[indoff+13] = 'induction_advec_phi'
lut[indoff+14] = 'induction_phi'
lut[indoff+15] = 'induction_diff_phi'

# <v> x <B> terms
lut[indoff+16] = 'induction_shear_vmbm_r'      # {<B> dot grad <v>}
lut[indoff+17] = 'induction_comp_vmbm_r'       # -{div dot <v>}<B>
lut[indoff+18] = 'induction_advec_vmbm_r'      # -{<v> dot grad <B>}
lut[indoff+19] = 'induction_vmbm_r'              # del cross {<v> x <B>}
lut[indoff+20] = 'induction_diff_bm_r'         # -del x (eta {del x <B>})
lut[indoff+21] = 'induction_shear_vmbm_theta'  # theta components
lut[indoff+22] = 'induction_comp_vmbm_theta'
lut[indoff+23] = 'induction_advec_vmbm_theta'
lut[indoff+24] = 'induction_vmbm_theta'
lut[indoff+25] = 'induction_diff_bm_theta'
lut[indoff+26] = 'induction_shear_vmbm_phi'    # phi components
lut[indoff+27] = 'induction_comp_vmbm_phi'
lut[indoff+28] = 'induction_advec_vmbm_phi'
lut[indoff+29] = 'induction_vmbm_phi'
lut[indoff+30] = 'induction_diff_bm_phi'

# <v> x B'
lut[indoff+31] = 'induction_shear_vmbp_r'      # {B' dot grad <v>}
lut[indoff+32] = 'induction_comp_vmbp_r'       # -{div dot <v>}B'
lut[indoff+33] = 'induction_advec_vmbp_r'      # -{<v> dot grad B'}
lut[indoff+34] = 'induction_vmbp_r'            # del cross {<v> x B'}
lut[indoff+35] = 'induction_diff_bp_r'         # -del x (eta {del x B'})
lut[indoff+36] = 'induction_shear_vmbp_theta'  # theta components
lut[indoff+37] = 'induction_comp_vmbp_theta'
lut[indoff+38] = 'induction_advec_vmbp_theta'
lut[indoff+39] = 'induction_vmbp_theta'
lut[indoff+40] = 'induction_diff_bp_theta'
lut[indoff+41] = 'induction_shear_vmbp_phi'    # phi components
lut[indoff+42] = 'induction_comp_vmbp_phi'
lut[indoff+43] = 'induction_advec_vmbp_phi'
lut[indoff+44] = 'induction_vmbp_phi'
lut[indoff+45] = 'induction_diff_bp_phi'

# v' x <B>
lut[indoff+46] = 'induction_shear_vpbm_r'      # {<B> dot grad v'}
lut[indoff+47] = 'induction_comp_vpbm_r'       # -{div dot v'}<B>
lut[indoff+48] = 'induction_advec_vpbm_r'      # -{v' dot grad <B>}
lut[indoff+49] = 'induction_vpbm_r'            # del cross {v' x <B>}
lut[indoff+50] = 'induction_shear_vpbm_theta'  # theta components
lut[indoff+51] = 'induction_comp_vpbm_theta'
lut[indoff+52] = 'induction_advec_vpbm_theta'
lut[indoff+53] = 'induction_vpbm_theta'
lut[indoff+54] = 'induction_shear_vpbm_phi'    # phi components
lut[indoff+55] = 'induction_comp_vpbm_phi'
lut[indoff+56] = 'induction_advec_vpbm_phi'
lut[indoff+57] = 'induction_vpbm_phi'

# v' x B'
lut[indoff+58] = 'induction_shear_vpbp_r'      # {B' dot grad v'}
lut[indoff+59] = 'induction_comp_vpbp_r'       # -{div dot v'}B'
lut[indoff+60] = 'induction_advec_vpbp_r'      # -{v' dot grad B'}
lut[indoff+61] = 'induction_vpbp_r'            # del cross {v' x B'}
lut[indoff+62] = 'induction_shear_vpbp_theta'  # theta components
lut[indoff+63] = 'induction_comp_vpbp_theta'
lut[indoff+64] = 'induction_advec_vpbp_theta'
lut[indoff+65] = 'induction_vpbp_theta'
lut[indoff+66] = 'induction_shear_vpbp_phi'    # phi components
lut[indoff+67] = 'induction_comp_vpbp_phi'
lut[indoff+68] = 'induction_advec_vpbp_phi'
lut[indoff+69] = 'induction_vpbp_phi'

#----------------
# custom magnetic
#----------------
custom_mag = 700
lut[custom_mag+1] = 'cross_helicity'  # v dot B
lut[custom_mag+2] = 'turb_cross_helicity'

##################################
# old values from early version of Rayleigh
#     various definitions should in theory
#     be the same as those listed above
##################################
alt_lut[1] = 'v_r'
alt_lut[2] = 'v_theta'
alt_lut[3] = 'v_phi'
alt_lut[4] = 'temperature' # this is entropy
alt_lut[5] = 'pressure'
alt_lut[6] = 'v_sq'
alt_lut[7] = 'kinetic_energy'
alt_lut[8] = 'gradt_r'
alt_lut[9] = 'cond_flux_r'
alt_lut[10] = 'zonal_ke'
alt_lut[11] = 'merid_ke'
alt_lut[12] = 'vol_heating'
alt_lut[13] = 'rhov_r'
alt_lut[14] = 'rhov_theta'
alt_lut[15] = 'rhov_phi'
alt_lut[16] = 'thermale_flux_radial'
alt_lut[17] = 'radial_ke'
alt_lut[18] = 'ke_flux_radial'
alt_lut[19] = 'enth_flux_radial'
alt_lut[20] = 'buoyancy_work'
alt_lut[21] = 'vort_r'
alt_lut[22] = 'vort_theta'
alt_lut[23] = 'vort_phi'
alt_lut[24] = 'enstrophy'
alt_lut[25] = 'amom_fluct_r'
alt_lut[26] = 'amom_fluct_theta'
alt_lut[27] = 'amom_dr_r'
alt_lut[28] = 'amom_dr_theta'
alt_lut[29] = 'amom_mean_r'
alt_lut[30] = 'amom_mean_theta'
alt_lut[31] = 'visc_flux_r'

vgv = 40
alt_lut[vgv+1] = 'v_grad_v_r'; alt_lut[vgv+2] = 'v_grad_v_theta'; alt_lut[vgv+3] = 'v_grad_v_phi'
alt_lut[vgv+4] = 'vp_grad_vm_r'; alt_lut[vgv+5] = 'vp_grad_vm_theta'; alt_lut[vgv+6] = 'vp_grad_vm_phi'
alt_lut[vgv+7] = 'vm_grad_vp_r'; alt_lut[vgv+8] = 'vm_grad_vp_theta'; alt_lut[vgv+9] = 'vm_grad_vp_phi'
alt_lut[vgv+10] = 'vp_grad_vp_r'; alt_lut[vgv+11] = 'vp_grad_vp_theta'; alt_lut[vgv+12] = 'vp_grad_vp_phi'
alt_lut[vgv+13] = 'vm_grad_vm_r'; alt_lut[vgv+14] = 'vm_grad_vm_theta'; alt_lut[vgv+15] = 'vm_grad_vm_phi'

alt_lut[99] = 'diagnostic1'; alt_lut[100] = 'diagnostic2'
alt_lut[101] = 'vr2'; alt_lut[102] = 'vt2'; alt_lut[103] = 'vp2'
alt_lut[104] = 'vr3'; alt_lut[105] = 'vt3'; alt_lut[106] = 'vp3'

alt_lut[201] = 'b_r'; alt_lut[202] = 'b_theta'; alt_lut[203] = 'b_phi'
alt_lut[204] = 'j_r'; alt_lut[205] = 'j_theta'; alt_lut[206] = 'j_phi'
alt_lut[207] = 'b_sq'
alt_lut[208] = 'magnetic_energy'; alt_lut[209] = 'zonal_me'; alt_lut[210] = 'merid_me'
alt_lut[211] = 'b_r2'; alt_lut[212] = 'b_theta2'; alt_lut[213] = 'b_phi2'

alt_lut[220] = 'j_cross_b_r'; alt_lut[221] = 'j_cross_b_theta'; alt_lut[222] = 'j_cross_b_phi'
alt_lut[223] = 'jp_cross_bm_r'; alt_lut[224] = 'jp_cross_bm_theta'; alt_lut[225] = 'jp_cross_bm_phi'
alt_lut[226] = 'jm_cross_bp_r'; alt_lut[227] = 'jm_cross_bp_theta'; alt_lut[228] = 'jm_cross_bp_phi'
alt_lut[229] = 'jm_cross_bm_r'; alt_lut[230] = 'jm_cross_bm_theta'; alt_lut[231] = 'jm_cross_bm_phi'
alt_lut[232] = 'jp_cross_bp_r'; alt_lut[233] = 'jp_cross_bp_theta'; alt_lut[234] = 'jp_cross_bp_phi'

indoff = 250
alt_lut[indoff+1]  = 'induction_shear_r'      # {B dot grad v}
alt_lut[indoff+2]  = 'induction_comp_r'       # -{div dot v}B
alt_lut[indoff+3]  = 'induction_advec_r'      # -{v dot grad B}
alt_lut[indoff+4]  = 'induction_r'            # del cross {v x B}
alt_lut[indoff+5]  = 'induction_diff_r'       # -del x (eta {del x B})
alt_lut[indoff+6]  = 'induction_shear_theta'  # theta components
alt_lut[indoff+7]  = 'induction_comp_theta'
alt_lut[indoff+8]  = 'induction_advec_theta'
alt_lut[indoff+9]  = 'induction_theta'
alt_lut[indoff+10] = 'induction_diff_theta'
alt_lut[indoff+11] = 'induction_shear_phi'    # phi components
alt_lut[indoff+12] = 'induction_comp_phi'
alt_lut[indoff+13] = 'induction_advec_phi'
alt_lut[indoff+14] = 'induction_phi'
alt_lut[indoff+15] = 'induction_diff_phi'

alt_lut[indoff+16] = 'induction_shear_vmbm_r'      # {<B> dot grad <v>}
alt_lut[indoff+17] = 'induction_comp_vmbm_r'       # -{div dot <v>}<B>
alt_lut[indoff+18] = 'induction_advec_vmbm_r'      # -{<v> dot grad <B>}
alt_lut[indoff+19] = 'induction_vmbm_r'              # del cross {<v> x <B>}
alt_lut[indoff+20] = 'induction_diff_bm_r'         # -del x (eta {del x <B>})
alt_lut[indoff+21] = 'induction_shear_vmbm_theta'  # theta components
alt_lut[indoff+22] = 'induction_comp_vmbm_theta'
alt_lut[indoff+23] = 'induction_advec_vmbm_theta'
alt_lut[indoff+24] = 'induction_vmbm_theta'
alt_lut[indoff+25] = 'induction_diff_bm_theta'
alt_lut[indoff+26] = 'induction_shear_vmbm_phi'    # phi components
alt_lut[indoff+27] = 'induction_comp_vmbm_phi'
alt_lut[indoff+28] = 'induction_advec_vmbm_phi'
alt_lut[indoff+29] = 'induction_vmbm_phi'
alt_lut[indoff+30] = 'induction_diff_bm_phi'

alt_lut[indoff+31] = 'induction_shear_vmbp_r'      # {B' dot grad <v>}
alt_lut[indoff+32] = 'induction_comp_vmbp_r'       # -{div dot <v>}B'
alt_lut[indoff+33] = 'induction_advec_vmbp_r'      # -{<v> dot grad B'}
alt_lut[indoff+34] = 'induction_vmbp_r'            # del cross {<v> x B'}
alt_lut[indoff+35] = 'induction_diff_bp_r'         # -del x (eta {del x B'})
alt_lut[indoff+36] = 'induction_shear_vmbp_theta'  # theta components
alt_lut[indoff+37] = 'induction_comp_vmbp_theta'
alt_lut[indoff+38] = 'induction_advec_vmbp_theta'
alt_lut[indoff+39] = 'induction_vmbp_theta'
alt_lut[indoff+40] = 'induction_diff_bp_theta'
alt_lut[indoff+41] = 'induction_shear_vmbp_phi'    # phi components
alt_lut[indoff+42] = 'induction_comp_vmbp_phi'
alt_lut[indoff+43] = 'induction_advec_vmbp_phi'
alt_lut[indoff+44] = 'induction_vmbp_phi'
alt_lut[indoff+45] = 'induction_diff_bp_phi'

alt_lut[indoff+46] = 'induction_shear_vpbm_r'      # {<B> dot grad v'}
alt_lut[indoff+47] = 'induction_comp_vpbm_r'       # -{div dot v'}<B>
alt_lut[indoff+48] = 'induction_advec_vpbm_r'      # -{v' dot grad <B>}
alt_lut[indoff+49] = 'induction_vpbm_r'            # del cross {v' x <B>}
alt_lut[indoff+50] = 'induction_shear_vpbm_theta'  # theta components
alt_lut[indoff+51] = 'induction_comp_vpbm_theta'
alt_lut[indoff+52] = 'induction_advec_vpbm_theta'
alt_lut[indoff+53] = 'induction_vpbm_theta'
alt_lut[indoff+54] = 'induction_shear_vpbm_phi'    # phi components
alt_lut[indoff+55] = 'induction_comp_vpbm_phi'
alt_lut[indoff+56] = 'induction_advec_vpbm_phi'
alt_lut[indoff+57] = 'induction_vpbm_phi'

alt_lut[indoff+58] = 'induction_shear_vpbp_r'      # {B' dot grad v'}
alt_lut[indoff+59] = 'induction_comp_vpbp_r'       # -{div dot v'}B'
alt_lut[indoff+60] = 'induction_advec_vpbp_r'      # -{v' dot grad B'}
alt_lut[indoff+61] = 'induction_vpbp_r'            # del cross {v' x B'}
alt_lut[indoff+62] = 'induction_shear_vpbp_theta'  # theta components
alt_lut[indoff+63] = 'induction_comp_vpbp_theta'
alt_lut[indoff+64] = 'induction_advec_vpbp_theta'
alt_lut[indoff+65] = 'induction_vpbp_theta'
alt_lut[indoff+66] = 'induction_shear_vpbp_phi'    # phi components
alt_lut[indoff+67] = 'induction_comp_vpbp_phi'
alt_lut[indoff+68] = 'induction_advec_vpbp_phi'
alt_lut[indoff+69] = 'induction_vpbp_phi'

