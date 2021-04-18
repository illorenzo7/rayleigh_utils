# routine to get groupings of Rayleigh quantities and also deal with units
# "building block" units
import numpy as np

def make_unit(st, exp=1):
    the_unit = r'$\rm{%s}$' %st
    if exp != 1:
        the_unit += r'$^{%i}$' %exp
    return the_unit

def get_quantity_group(tag, magnetism):
    di_out = dict({})

    # set default qvals: velocity + vorticity, Pressure/ entropy
    # then possibly B, del x B
    if tag == 'default':
        qvals = [1, 2, 3, 301, 302, 303, 501, 502]
        titles = [r'$v_r$', r'$v_\theta$', r'$v_\phi$',\
                r'$\omega_r$', r'$\omega_\theta$', r'$\omega_\phi$',\
                r'$S$', r'$P$']

        v_unit = make_unit('cm') + make_unit('s', -1)
        om_unit = make_unit('s', -1)
        units = [v_unit]*3 + [om_unit]*3 +\
            [make_unit('erg') + make_unit('g', -1) + make_unit('K', -1),\
            make_unit('dyn') + make_unit('cm', -2)]

        magnetism = get_parameter(dirname, 'magnetism')
        if magnetism:
            qvals += [801, 802, 803, 1001, 1002, 1003]
            titles += [r'$B_r$', r'$B_\theta$', r'$B_\phi$',\
                r'$\mathcal{J}_r$', r'$\mathcal{J}_\theta$', r'$\mathcal{J}_\phi$']
            b_unit = make_unit('G')
            j_unit = make_unit('G') + make_unit('cm', -1)
            units += [b_unit]*3 + [j_unit]*3

    if tag == 'torque':
        qvals = [3, 1801, 1802, 1803, 1804, 1819]
        titles = [r'$v_\phi$', r'$-\tau_{\rm{rs}}$', r'$-\tau_{\rm{mc,v_\phi}}$',r'$\tau_{\rm{mc,\Omega_0}}$', r'$\tau_{\rm{v}}$',\
          r'$\mathcal{L}_z$']
        units = [make_unit('cm') + make_unit('s', -1)]
        torque_unit = make_unit('g') + make_unit('cm', -1) +\
                make_unit('s', -2)
        for i in range(4):
            units.append(torque_unit)
        units.append(make_unit('g') + make_unit('cm', -1) +\
                make_unit('s', -1))

        if magnetism:
            qvals.append(1805)
            qvals.append(1806)
            titles.append(r'$\tau_{\rm{mm}}$')
            titles.append(r'$\tau_{\rm{ms}}$')
            units.append(torque_unit)
            units.append(torque_unit)

    if tag == 'induction':
        qvals = [801, 802, 803]            
        for j in range(1, 31):
            qvals.append(1600 + j)
    di_out['qvals'] = np.array(qvals)
    di_out['titles'] = np.array(titles)
    di_out['units'] = np.array(units)

    return di_out
