
    # Read in the Shell_Spectra data
    if the_file_spec is None:
        the_file_spec = get_widest_range_file(datadir, 'Shell_Spectra')

    print ('length_scales(): reading ' + the_file)
    di = get_dict(the_file_spec)
    lpower = di['lpower']
    nell = np.shape(lpower)[0]
    lvals = np.arange(nell)
    lvals = lvals.reshape((nell, 1))
    lut = di['lut']
    rr_spec = di['rvals']

    # get the full power
    vsqr = lpower[:, :, lut[1], 0]
    vsqt = lpower[:, :, lut[2], 0] 
    vsqp = lpower[:, :, lut[3], 0] 

    vsq = vsqr + vsqt + vsqp
    vsqhor = vsqt + vsqp
    vsqpol = vsqr + vsqt

    # get the mean (m=0) power    
    vsqr_mean = lpower[:, :, lut[1], 1]
    vsqt_mean = lpower[:, :, lut[2], 1] 
    vsqp_mean = lpower[:, :, lut[3], 1] 

    vsq_mean = vsqr_mean + vsqt_mean + vsqp_mean
    vsqhor_mean = vsqt_mean + vsqp_mean
    vsqpol_mean = vsqr_mean + vsqt_mean

    # get the fluc (convective) power    
    vsqr_fluc = lpower[:, :, lut[1], 2]
    vsqt_fluc = lpower[:, :, lut[2], 2] 
    vsqp_fluc = lpower[:, :, lut[3], 2] 

    vsq_fluc = vsqr_fluc + vsqt_fluc + vsqp_fluc
    vsqhor_fluc = vsqt_fluc + vsqp_fluc
    vsqpol_fluc = vsqr_fluc + vsqt_fluc

    # Compute rms l-values
    
    # full v
    l_rms_vr = np.sum(vsqr*lvals, axis=0)/np.sum(vsqr, axis=0)
    l_rms_vh = np.sum(vhsq_power*(lvals + 1.), axis=0)/np.sum(vhsq_power, axis=0)
    l_rms_v = np.sum(vsq_power*(lvals + 1.), axis=0)/np.sum(vsq_power, axis=0)
    # Compute lengthscales and add to dictionary
    pir = np.pi*rr_spec
    L_vr = pir/l_rms_vr
    L_vh = pir/l_rms_vh
    L_v = pir/l_rms_v
    di_out['rr_spec'] = rr_spec
    di_out['ir_spec'] = inds_from_vals(rr, di_out['rr_spec'])
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

