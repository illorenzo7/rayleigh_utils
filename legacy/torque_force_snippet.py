    
if forced: # compute torque associated with forcing function
    mean_vp = vals[:, :, lut[3]]
    tacho_r = get_parameter(dirname, 'tacho_r')
    print ("read tacho_r = %1.2e" %tacho_r)
    tacho_dr = get_parameter(dirname, 'tacho_dr')
    tacho_tau = get_parameter(dirname, 'tacho_tau')
    forcing = np.zeros((nt, nr))
    eq = get_eq(dirname)
    rho = eq.rho

    if os.path.exists(dirname + '/inner_vp'):
        print ("inner_vp file exists, so I assume you have a forcing function which\n quartically matches on to a CZ differential rotation")
        inner_vp = read_inner_vp(dirname + '/inner_vp')
        for it in range(nt):
            for ir in range(nr):
                if rr[ir] <= tacho_r - tacho_dr*rr[0]:
                    desired_vp = 0.0
                elif rr[ir] <= tacho_r:
                    desired_vp = inner_vp[it]*(1.0 - ( (rr[ir] - tacho_r)/(tacho_dr*rr[0]) )**2)**2
                else:
                    desired_vp = mean_vp[it, ir] # No forcing in CZ
                forcing[it, ir] = -rho[ir]*(mean_vp[it, ir] - desired_vp)/tacho_tau
    elif os.path.exists(dirname + '/eq_vp'):
        print ("eq_vp file exists, so I assume you have a forcing function which\n quartically matches on to a CZ differential rotation\n with viscous-torque-free buffer zone")
        tacho_r2 = get_parameter(dirname, 'tacho_r2')
        i_tacho_r = np.argmin(np.abs(rr - tacho_r))
        print ("read tacho_r2 = %1.2e" %tacho_r2)
        eq_vp = read_eq_vp(dirname + '/eq_vp', nt, nr)
        for it in range(nt):
            for ir in range(nr):
                if rr[ir] <= tacho_r - tacho_dr*rr[0]:
                    desired_vp = 0.0
                elif rr[ir] <= tacho_r:
                    desired_vp = eq_vp[it, i_tacho_r]*(1.0 - ( (rr[ir] - tacho_r)/(tacho_dr*rr[0]) )**2)**2
                elif rr[ir] <= tacho_r2:
                    desired_vp = eq_vp[it, ir]
                else:
                    desired_vp = mean_vp[it, ir] # No forcing in CZ
                forcing[it, ir] = -rho[ir]*(mean_vp[it, ir] - desired_vp)/tacho_tau
    else:
        forcing_coeff = -rho/tacho_tau*0.5*(1.0 - np.tanh((rr - tacho_r)/(tacho_dr*rr[0])))
        forcing = forcing_coeff.reshape((1, nr))*mean_vp

    # convert forcing function into a torque
    torque_forcing = di['xx']*forcing
    torque_tot += torque_forcing

torques = [torque_rs, torque_mc, torque_visc, torque_tot]
titles = [r'$\tau_{\rm{rs}}$', r'$\tau_{\rm{mc}}$', r'$\tau_{\rm{v}}$',\
          r'$\tau_{\rm{tot}}$']
units = r'$\rm{g}\ \rm{cm}^{-1}\ \rm{s}^{-2}$'

ind_insert = 3
if magnetism:
    torques.insert(ind_insert, torque_Maxwell_mean)
    torques.insert(ind_insert + 1, torque_Maxwell_rs)
    titles.insert(ind_insert, r'$\tau_{\rm{mm}}$') 
    titles.insert(ind_insert + 1, r'$\tau_{\rm{ms}}$')
    ind_insert += 2
if forced:
    torques.insert(ind_insert, torque_forcing)
    titles.insert(ind_insert, r'$\tau_{\rm{forcing}}$')
    ind_insert =+ 1
    
if nadd > 0:
    jstart = 0
    for i in range(nadd):
        ind = torques_to_add[jstart]
        torque_sum = np.copy(torques[ind])
        title_sum = titles[ind]
        for j in range(jstart + 1, jstart + nsubset[i]):
            ind = torques_to_add[j]
            torque_sum += torques[ind]
            title_sum += r'$\ +\ $' + titles[ind]
        jstart += nsubset[i]

        torques.insert(ind_insert, torque_sum)
        titles.insert(ind_insert, title_sum)
        ind_insert += 1


