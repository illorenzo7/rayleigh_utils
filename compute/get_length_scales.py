# Author: Loren Matilsky
# Created: 04/10/2020
# This script computes various lengthscales associated with
# Rayleigh run in directory [dirname]:
# SPECTRAL LENGTH SCALES (if Shell_Spectra file available)
# L_vh, L_vr, and L_v: pi*r / rms l-value for horizontal, vertical, and total convective flows
# L_Bm, L_Bp, and  L_B: pi*r / rms l-value for meridional, longitudinal, and total convective magnetic fields
# L_S: pi*r / rms l-value of convective entropy
# FLOW LENGTH SCALE (if Shell_Avgs file available)
# L_om: < (v')^2 / (om')^2 >_sph^(1/2)
# MIXING LENGTH SCALE
# L_rho: -(dlnhro/dr)^(-1)
# Returns length scales as a dictionary with arrays

import numpy as np
import sys, os
sys.path.append(os.environ['rapp'])
sys.path.append(os.environ['raco'])
from rayleigh_diagnostics import Shell_Avgs
from common import *

def get_length_scales(dirname):
    # Make empty dictionary for length_scale arrays
    di_out = dict([])

    # See if run is magnetic
    magnetism = get_parameter(dirname, 'magnetism')

    # First get mixing length scale (the one we know will be there)
    eq = get_eq(dirname)
    L_rho = -1./eq.dlnrho
    rr = eq.radius
    nr = len(rr)
    di_out['rr'] = rr
    di_out['nr'] = nr
    di_out['L_rho'] = L_rho

    # Get data directory
    datadir = dirname + '/data/'

    # Read in the Shell_Avgs data
    the_file = get_widest_range_file(datadir, 'Shell_Avgs')
    if not the_file == '':
        print ('length_scales(): Getting velocity/vorticity from ' +\
                the_file)
        di = get_dict(datadir + the_file)
        di_out['iter1'], di_out['iter2'] = get_iters_from_file(datadir +\
                the_file)
        vals = di['vals']
        lut = di['lut']
        try:
            # Read in enstrophy of convective flows
            vortsqr = vals[:, lut[317]] 
            vortsqt = vals[:, lut[318]]
            vortsqp = vals[:, lut[319]]
            vortsqh = vortsqt + vortsqp
            enstr = vortsqr + vortsqt + vortsqp
            # Read in velocity-squared of convective flows
            vsq = vals[:, lut[422]] + vals[:, lut[423]] + vals[:, lut[424]]
            # Compute length scale and put it in dictionary
            L_omr = (vsq/vortsqr)**0.5
            L_omh = (vsq/vortsqh)**0.5
            L_om = (vsq/enstr)**0.5
            di_out['L_omr'] = L_omr
            di_out['L_omh'] = L_omh
            di_out['L_om'] = L_om
        except:
            print ("length_scales(): one or more quantities needed for enstrophy or")
            print("velocity-squared were not output for Shell_Avgs data")
            print("failed to compute L_om")

        if magnetism:
            print ('Getting B fields/currents from ' + the_file)
            try:
                # Read in current of convective fields
                del_crossB2 = vals[:, lut[1015]] + vals[:, lut[1018]] +\
                        vals[:, lut[1021]]
                # Read in B-squared of convective flows
                B2 = 8.*np.pi*(vals[:, lut[1110]] + vals[:, lut[1111]] +\
                        vals[:, lut[1112]])
                # Compute length scale and put it in dictionary
                L_J = (B2/del_crossB2)**0.5
                di_out['L_J'] = L_J
            except:
                print ("one or more quantities needed for current or")
                print("B-squared were not output for Shell_Avgs data")
                print("failed to compute L_J, assigning it L_om ")
                di_out['L_J'] = L_om

    # Read in the Shell_Spectra data
    the_file = get_widest_range_file(datadir, 'Shell_Spectra')
    if not the_file == '':
        print ('length_scales(): Reading Shell_Spectra data from ' +\
                the_file)
        di = get_dict(datadir + the_file)
        lpower = di['lpower']
        rinds = di['rinds']
        nell = di['nell']
        lvals = di['lvals'].reshape((nell, 1))
        lut = di['lut']
        rr_spec = di['rvals']
        # get the convective power    
        vrsq_power = lpower[:, :, lut[1], 2]
        vtsq_power = lpower[:, :, lut[2], 2] 
        vpsq_power = lpower[:, :, lut[3], 2] 
        vhsq_power = vtsq_power + vpsq_power
        vsq_power = vrsq_power + vtsq_power + vpsq_power
        # Compute rms l-values
        l_rms_vr = np.sum(vrsq_power*(lvals + 1.), axis=0)/np.sum(vrsq_power, axis=0)
        l_rms_vh = np.sum(vhsq_power*(lvals + 1.), axis=0)/np.sum(vhsq_power, axis=0)
        l_rms_v = np.sum(vsq_power*(lvals + 1.), axis=0)/np.sum(vsq_power, axis=0)
        # Compute lengthscales and add to dictionary
        pir = np.pi*rr_spec
        L_vr = pir/l_rms_vr
        L_vh = pir/l_rms_vh
        L_v = pir/l_rms_v
        di_out['rr_spec'] = rr_spec
        di_out['nr_spec'] = len(rr_spec)
        di_out['L_vr'] = L_vr
        di_out['L_vh'] = L_vh
        di_out['L_v'] = L_v

        if magnetism:
            Brsq_power = lpower[:, :, lut[801], 2]
            Btsq_power = lpower[:, :, lut[802], 2] 
            Bpsq_power = lpower[:, :, lut[803], 2] 
            Bmsq_power = Brsq_power + Btsq_power
            Bsq_power = Brsq_power + Btsq_power + Bpsq_power
            # Compute rms l-values
            l_rms_Bp = np.sum(Bpsq_power*(lvals + 1.), axis=0)/np.sum(Bpsq_power, axis=0)
            l_rms_Bm = np.sum(Bmsq_power*(lvals + 1.), axis=0)/np.sum(Bmsq_power, axis=0)
            l_rms_B = np.sum(Bsq_power*(lvals + 1.), axis=0)/np.sum(Bsq_power, axis=0)

            # Compute lengthscales and add to dictionary
            L_Bp = pir/l_rms_Bp
            L_Bm = pir/l_rms_Bm
            L_B = pir/l_rms_B
            di_out['L_Bp'] = L_Bp
            di_out['L_Bm'] = L_Bm
            di_out['L_B'] = L_v
    # For consistency also compute the shell depth
    di_out['shell_depth'] = np.max(rr) - np.min(rr)

    # Return the dictionary 
    return di_out
