# Author: Loren Matilsky
# Created: 12/19/2018
# Stores properties of various Rayliegh output variables;
# units, quantity codes, LaTex variable names, etc.
# also groups of quantities, like torques, induction, etc.
# Add to these lists as need be.
# Not all the variables have quantity codes but are derivative quantities
# of other fluid variables
import numpy as np
from common import array_of_strings

def make_unit(st, exp=1):
    the_unit = r'$\rm{%s}$' %st
    if exp != 1:
        the_unit += r'$^{%i}$' %exp
    return the_unit

# basic units
utype = dict({})
space = ' '
# field variables
utype['v'] = make_unit('cm') + make_unit('s', -1)
utype['om'] = make_unit('s', -1)
utype['b'] = make_unit('G')
utype['j'] = make_unit('G') + make_unit('cm', -1)
utype['s'] = make_unit('erg') + make_unit('g', -1) + make_unit('K', -1)
utype['p'] = make_unit('dyn') + make_unit('cm', -2)
utype['rho'] = make_unit('g') + make_unit('cm', -3)
utype['t'] = make_unit('K')

# derived units
utype['energy'] = make_unit('erg') + space + make_unit('cm', -3)
utype['eprod'] = utype['energy'] + space + make_unit('s', -1)
utype['power'] = make_unit('erg') + make_unit('s', -1)
utype['eflux'] = make_unit('erg') + make_unit('cm', -2) +\
        make_unit('s', -1)
utype['dsdt'] = make_unit('erg') + make_unit('g', -1) +\
        make_unit('K', -1) + make_unit('s', -1)
utype['dsdr'] = make_unit('erg') + make_unit('g', -1) +\
        make_unit('K', -1) + make_unit('cm', -1)
utype['amom'] = make_unit('g') + make_unit('cm', -1) +\
        make_unit('s', -1)
utype['aflux'] = make_unit('g') + make_unit('cm', -3) +\
        make_unit('s', -2)
utype['force'] = make_unit('g') + make_unit('cm', -2) +\
        make_unit('s', -2)
utype['torque'] = make_unit('g') + make_unit('cm', -1) +\
        make_unit('s', -2)
utype['induct'] = make_unit('G') + make_unit('s', -1)

# groups of quantities
def get_quantity_group(tag, magnetism):
    di_out = dict({})

    # set default qvals: velocity + vorticity, Pressure/ entropy
    # then possibly B, del x B
    if tag == 'default':
        qvals = [1, 2, 3, 301, 302, 303, 501, 502]
        titles = [r'$v_r$', r'$v_\theta$', r'$v_\phi$',\
                r'$\omega_r$', r'$\omega_\theta$', r'$\omega_\phi$',\
                r'$S$', r'$P$']

        units = [utype['v']]*3 + [utype['om']]*3 + [utype['s'], utype['p']]

        if magnetism:
            qvals += [801, 802, 803, 1001, 1002, 1003]
            titles += [r'$B_r$', r'$B_\theta$', r'$B_\phi$',\
                r'$\mathcal{J}_r$', r'$\mathcal{J}_\theta$', r'$\mathcal{J}_\phi$']
            units += [utype['b']]*3 + [utype['j']]*3

        ncol = 4 + 3*magnetism

    if tag == 'torque':
        qvals = [3, 1801, 1802, 1803, 1804, 1819]
        titles = [r'$v_\phi$', r'$-\tau_{\rm{rs}}$', r'$-\tau_{\rm{mc,v_\phi}}$',r'$\tau_{\rm{mc,\Omega_0}}$', r'$\tau_{\rm{v}}$',\
          r'$\mathcal{L}_z$']
        units = [utype['v']] + [utype['torque']]*4 + [utype['amom']]

        if magnetism:
            qvals += [1805, 1806]
            titles += [r'$\tau_{\rm{mm}}$', r'$\tau_{\rm{ms}}$']
            units += [utype['torque']]*2
        ncol = 3 + magnetism

    if tag == 'induct':
        qvals = [801, 802, 803]            
        for j in range(1, 31):
            qvals.append(1600 + j)
        titles = array_of_strings(qvals)
        units = [utype['induct']]
        ncol = 11

    if tag == 'v':
        qvals = [1, 2, 3]            
        titles = [texlabels['vr'], texlabels['vt'], texlabels['vp']]
        units = [utype['b']]
        ncol = 3

    if tag == 'b':
        qvals = [801, 802, 803]            
        titles = [texlabels['br'], texlabels['bt'], texlabels['bp']]
        units = [utype['b']]
        ncol = 3

    if tag == 'forcer': # linear forces, radial
        qvals = [1201, 1219, 1237, 1216, 1228]
        titles= [r'$-(\mathbf{f}_{\rm{adv}})_r$', r'$(\mathbf{f}_{\rm{cor}})_r$',  r'$(\mathbf{f}_{\rm{p}})_r$', r'$(\mathbf{f}_{\rm{buoy}})_r$', r'$(\mathbf{f}_{\rm{v}})_r$']
        units = [utype['force']]
        if magnetism:
            qvals += [1248]
            titles += [r'$(\mathbf{f}_{\rm{mag}})_r$']
        ncol = 3

    if tag == 'forcet': # linear forces, theta
        qvals = [1202, 1220, 1238, 1229]
        titles= [r'$-(\mathbf{f}_{\rm{adv}})_\theta$', r'$(\mathbf{f}_{\rm{cor}})_\theta$',  r'$(\mathbf{f}_{\rm{p}})_\theta$',  r'$(\mathbf{f}_{\rm{v}})_\theta$']
        units = [utype['force']]
        if magnetism:
            qvals += [1249]
            titles += [r'$(\mathbf{f}_{\rm{mag}})_\theta$']
        ncol = 3

    if tag == 'forcep': # linear forces, phi
        qvals = [1203, 1221, 1239, 1230]
        titles= [r'$-(\mathbf{f}_{\rm{adv}})_\phi$', r'$(\mathbf{f}_{\rm{cor}})_\phi$',  r'$(\mathbf{f}_{\rm{p}})_\phi$',  r'$(\mathbf{f}_{\rm{v}})_\phi$']
        units = [utype['force']]
        if magnetism:
            qvals += [1250]
            titles += [r'$(\mathbf{f}_{\rm{mag}})_\phi$']
        ncol = 3

    if tag == 'efr': # energy fluxes, r
        qvals = [1455, 1458, 1470, 1935, 1923]
        titles = [r'$(\mathbf{\mathcal{F}}_{\rm{enth}})_r$',\
          r'$(\mathbf{\mathcal{F}}_{\rm{enth,pp}})_r$',\
          r'$(\mathbf{\mathcal{F}}_{\rm{cond}})_r$',\
          r'$-(\mathbf{\mathcal{F}}_{\rm{visc}})_r$',\
          r'$(\mathbf{\mathcal{F}}_{\rm{KE}})_r$']
        units = [utype['eflux']]
        if magnetism:
            qvals += [2001]
            titles += [r'$(\mathbf{\mathcal{F}}_{\rm{Poynt}})_r$']
        ncol = 3

    if 'meprodnum' in tag:
        nq = 12 # (r, t, p) x (ind, shear, adv, comp)
        ext = tag[-3:]
        if ext == 'tot':
            iqstart = 0
        if ext == 'pmp':
            iqstart = nq
        if ext == 'ppm':
            iqstart = 2*nq
        if ext == 'mmm':
            iqstart = 3*nq
        if ext == 'mpp':
            iqstart = 4*nq
        if ext == 'ppp':
            iqstart = 5*nq
        qvals = np.arange(iqstart, iqstart + nq)
        titles = []
        for direc in ['r', 'th', 'ph']:
            app = ' (' + direc + ')'
            titles += ['induct' + app, 'shear' + app, 'advec' + app,\
                    'comp' + app]
        units = [utype['eprod']]
        ncol = 4

    if tag == 'eft': # energy fluxes, theta
        qvals = [1456, 1459, 1471, 1936, 1924]
        titles = [r'$(\mathbf{\mathcal{F}}_{\rm{enth}})_\theta$',\
          r'$(\mathbf{\mathcal{F}}_{\rm{enth,pp}})_\theta$',\
          r'$(\mathbf{\mathcal{F}}_{\rm{cond}})_\theta$',\
          r'$-(\mathbf{\mathcal{F}}_{\rm{visc}})_\theta$',\
          r'$(\mathbf{\mathcal{F}}_{\rm{KE}})_\theta$']
        units = [utype['eflux']]
        if magnetism:
            qvals += [2002]
            titles += [r'$(\mathbf{\mathcal{F}}_{\rm{Poynt}})_\theta$']
        ncol = 3

    if tag == 'efp': # energy fluxes, phi
        qvals = [1457, 1460, 1937, 1925]
        titles = [r'$(\mathbf{\mathcal{F}}_{\rm{enth}})_\phi$',\
          r'$(\mathbf{\mathcal{F}}_{\rm{enth,pp}})_\phi$',\
          r'$-(\mathbf{\mathcal{F}}_{\rm{visc}})_\phi$',\
          r'$(\mathbf{\mathcal{F}}_{\rm{KE}})_\phi$']
        units = [utype['eflux']]
        if magnetism:
            qvals += [2003]
            titles += [r'$(\mathbf{\mathcal{F}}_{\rm{Poynt}})_\phi$']
        ncol = 3

    if tag == 'indr': # induction, r
        qvals = [1601, 1602, 1603, 1604, 1605, 801]
        titles = [r'$[\left\langle\mathbf{B}\cdot\nabla\mathbf{v}\right\rangle]_r$',\
            r'$-\left\langleB_r(\nabla\cdot\mathbf{v})\right\rangle$',\
            r'$-[\left\langle\mathbf{v}\cdot\nabla\mathbf{B}\right\rangle]_r$',\
            r'$[\nabla\times(\left\langle\mathbf{v}\times\mathbf{B}\right\rangle)]_r$',\
            r'$-[\nabla\times(\eta\nabla\times\langle\mathbf{B}\rangle)]_r$',\
            r'$\left\langle B_r\right\rangle$']
        units = [utype['induct']]*5 + [utype['b']]
        ncol = 3

    if tag == 'indrmean': # energy fluxes, r, mean
        qvals = [1616, 1617, 1618, 1619, 1620, 801]
        titles = [r'$[\left\langle\mathbf{B}\right\rangle\cdot\nabla\left\langle\mathbf{v}\right\rangle]_r$',\
    r'$-\left\langleB_r\right\rangle(\nabla\cdot\left\langle\mathbf{v}\right\rangle)$',\
    r'$-[\left\langle\mathbf{v}\right\rangle\cdot\nabla\left\langle\mathbf{B}\right\rangle]_r$',\
    r'$[\nabla\times(\left\langle\mathbf{v}\right\rangle\times\left\langle\mathbf{B}\right\rangle)]_r$',\
    r'$-[\nabla\times(\eta\nabla\times\langle\mathbf{B}\rangle)]_r$',\
    r'$\left\langleB_r\right\rangle$']
        units = [utype['induct']]*5 + [utype['b']]
        ncol = 3

    if tag == 'indt': # energy fluxes, theta
        qvals = [1606, 1607, 1608, 1609, 1610, 802]
        titles = [r'$[\left\langle\mathbf{B}\cdot\nabla\mathbf{v}\right\rangle]_\theta$',\
            r'$-\left\langleB_\theta(\nabla\cdot\mathbf{v})\right\rangle$',\
            r'$-[\left\langle\mathbf{v}\cdot\nabla\mathbf{B}\right\rangle]_\theta$',\
            r'$[\nabla\times(\left\langle\mathbf{v}\times\mathbf{B}\right\rangle)]_\theta$',\
            r'$-[\nabla\times(\eta\nabla\times\langle\mathbf{B}\rangle)]_\theta$',\
            r'$\left\langle B_\theta\right\rangle$']
        units = [utype['induct']]*5 + [utype['b']]
        ncol = 3

    if tag == 'indtmean': # energy fluxes, theta, mean
        qvals = [1621, 1622, 1623, 1624, 1625, 802]
        titles = [r'$[\left\langle\mathbf{B}\right\rangle\cdot\nabla\left\langle\mathbf{v}\right\rangle]_\theta$',\
            r'$-\left\langleB_\theta\right\rangle(\nabla\cdot\left\langle\mathbf{v}\right\rangle)$',\
            r'$-[\left\langle\mathbf{v}\right\rangle\cdot\nabla\left\langle\mathbf{B}\right\rangle]_\theta$',\
            r'$[\nabla\times(\left\langle\mathbf{v}\right\rangle\times\left\langle\mathbf{B}\right\rangle)]_\theta$',\
            r'$-[\nabla\times(\eta\nabla\times\langle\mathbf{B}\rangle)]_\theta$',\
            r'$\left\langleB_\theta\right\rangle$']
        units = [utype['induct']]*5 + [utype['b']]
        ncol = 3

    if tag == 'indp': # energy fluxes, phi
        qvals = [1611, 1612, 1613, 1614, 1615, 803]

        titles = [r'$[\left\langle\mathbf{B}\cdot\nabla\mathbf{v}\right\rangle]_\phi$',\
        r'$-\left\langleB_\phi(\nabla\cdot\mathbf{v})\right\rangle$',\
        r'$-[\left\langle\mathbf{v}\cdot\nabla\mathbf{B}\right\rangle]_\phi$',\
        r'$[\nabla\times(\left\langle\mathbf{v}\times\mathbf{B}\right\rangle)]_\phi$',\
        r'$-[\nabla\times(\eta\nabla\times\langle\mathbf{B}\rangle)]_\phi$',\
        r'$\left\langle B_\phi\right\rangle$']
        units = [utype['induct']]*5 + [utype['b']]
        ncol = 3

    if tag == 'indpmean': # energy fluxes, phi, mean
        qvals = [1626, 1627, 1628, 1629, 1630, 803]
        titles = [r'$[\left\langle\mathbf{B}\right\rangle\cdot\nabla\left\langle\mathbf{v}\right\rangle]_\phi$',\
            r'$-\left\langleB_\phi\right\rangle(\nabla\cdot\left\langle\mathbf{v}\right\rangle)$',\
            r'$-[\left\langle\mathbf{v}\right\rangle\cdot\nabla\left\langle\mathbf{B}\right\rangle]_\phi$',\
            r'$[\nabla\times(\left\langle\mathbf{v}\right\rangle\times\left\langle\mathbf{B}\right\rangle)]_\phi$',\
            r'$-[\nabla\times(\eta\nabla\times\langle\mathbf{B}\rangle)]_\phi$',\
            r'$\left\langleB_\phi\right\rangle$']
        units = [utype['induct']]*5 + [utype['b']]
        ncol = 3

    if tag == 'ke':
        qvals = [402, 403, 404, 410, 411, 412]
        titles =\
            [r'$\overline{\rm{KE}}_{\rm{r}}$',\
            r'$\overline{\rm{KE}}_{\rm{\theta}}$',\
            r'$\overline{\rm{KE}}_{\rm{\phi}}$',\
            r'$\rm{KE}^\prime_r$',\
            r'$\rm{KE}^\prime_\theta$',\
            r'$\rm{KE}^\prime_\phi$']
        units = [utype['energy']]
        ncol = 3
    
    if tag == 'me':
        qvals = [1102, 1103, 1104, 1110, 1111, 1112]
        titles =\
            [r'$\overline{\rm{ME}}_{\rm{r}}$',\
            r'$\overline{\rm{ME}}_{\rm{\theta}}$',\
            r'$\overline{\rm{ME}}_{\rm{\phi}}$',\
            r'$\rm{ME}^\prime_r$',\
            r'$\rm{ME}^\prime_\theta$',\
            r'$\rm{ME}^\prime_\phi$']
        units = [utype['energy']]
        ncol = 3

    if tag[:-3] == 'meprod': # this is the exact stuff
        ext = tag[-3:]
        basetitles = ['induct', 'shear', 'advec', 'comp', 'diff']

        custom_offset = 2200
        set_offset = 15
        set_offset2 = 12
        set_offset3 = 9

        if ext == 'tot':
            ncol = 5
            iqstart = custom_offset + 1

        if ext == 'pmp':
            ncol = 4
            iqstart = custom_offset + 1 + set_offset

        if ext == 'ppm':
            ncol = 4
            iqstart = custom_offset + 1 + set_offset + set_offset2

        if ext == 'mmm':
            ncol = 5
            iqstart = custom_offset + 1 + set_offset + 2*set_offset2

        if ext == 'mpp':
            ncol = 3
            iqstart = custom_offset + 1 + 2*set_offset + 2*set_offset2

        if ext == 'ppp':
            ncol = 5
            iqstart = custom_offset + 1 + 2*set_offset + 2*set_offset2 + set_offset3

        # do r, theta, then phi, production terms, in that order
        qvals = []
        titles = []
        count = 0
        for direc in ['r', 'th', 'ph']:
            qvals_loc = iqstart + np.arange(0, ncol*3, 3) + count
            qvals += qvals_loc.tolist()
            titles_loc = []
            for j in range(len(basetitles)):
                titles_loc.append(basetitles[j] + ' (' + ext + ', ' +\
                        direc + ')')
            titles += titles_loc
            count += 1
        units = [utype['eprod']]
        ncol = 5

    if tag in ['indralt', 'indraltnum']:
        ind_off = 0
        qvals = np.arange(15) + ind_off
        titles = ['(d/dT)(vr*bt)', '-(d/dT)(vt*br)', '-(d/dP)(vp*br)', '(d/dP)(vr*bp)', 'curv.']
        units = [utype['induct']]

    if tag in ['indtalt', 'indtaltnum']:
        ind_off = 15
        qvals = np.arange(15) + ind_off
        titles = ['(d/dP)(vt*bp)', '-(d/dP)(vp*bt)', '-(d/dr)(vr*bt)', '(d/dr)(vt*br)', 'curv.']
        units = [utype['induct']]

    if tag in ['indpalt', 'indpaltnum']:
        ind_off = 30
        qvals = np.arange(15) + ind_off
        titles = ['(d/dr)(vp*br)', '-(d/dr)(vr*bp)', '-(d/dT)(vt*bp)', '(d/dT)(vp*bt)', 'curv.']
        units = [utype['induct']]

    if 'ind' in tag and 'alt' in tag:
        for ext in ['mm', 'pp']:
            more_titles = []
            for j in range(5):
                more_titles.append(titles[j] + '_' + ext)
            titles += more_titles
        ncol = 5 

    if tag == 'magtorquemm':
        titles = [r'$\tau_{\rm{mm,r}}$', r'$\frac{r\sin\theta}{4\pi}\left\langle B_r\right\rangle\left\langle\frac{\partial B_\phi}{\partial r}\right\rangle$', r'$\frac{r\sin\theta}{4\pi}\left\langle B_\phi\right\rangle\left\langle\frac{\partial B_r}{\partial r}\right\rangle$', r'$\frac{3\sin\theta}{4\pi}\langle B_r\rangle\langle B_\phi\rangle$',\
        r'$\tau_{\rm{mm,\theta}}$', r'$\frac{\sin\theta}{4\pi}\left\langle B_\theta\right\rangle \left\langle\frac{\partial B_\phi}{\partial \theta}\right\rangle$', r'$\frac{\sin\theta}{4\pi}\left\langle B_\phi\right\rangle\left\langle\frac{\partial B_\theta}{\partial \theta}\right\rangle$', r'$\frac{2\cos\theta}{4\pi}\langle B_\theta\rangle\langle B_\phi\rangle$']
        qvals = [2, 3, 6, 7, 8, 9, 10, 11]
        units = [utype['torque']]
        ncol = 4

    if tag == 'magtorquems':
        titles = [r'$\tau_{\rm{ms,r}}$', r'$\frac{r\sin\theta}{4\pi}\left\langle B_r^\prime\frac{\partial B_\phi^\prime}{\partial r}\right\rangle$', r'$\frac{r\sin\theta}{4\pi}\left\langle B_\phi^\prime\frac{\partial B_r^\prime}{\partial r}\right\rangle$', r'$\frac{3\sin\theta}{4\pi}\langle B_r^\prime B_\phi^\prime\rangle$',\
        r'$\tau_{\rm{ms,\theta}}$', r'$\frac{\sin\theta}{4\pi}\left\langle B_\theta^\prime\frac{\partial B_\phi^\prime}{\partial \theta}\right\rangle$', r'$\frac{\sin\theta}{4\pi}\left\langle B_\phi^\prime\frac{\partial B_\theta^\prime}{\partial \theta}\right\rangle$', r'$\frac{2\cos\theta}{4\pi}\langle B_\theta^\prime B_\phi^\prime\rangle$']
        qvals = [4, 5, 12, 13, 14, 15, 16, 17]
        units = [utype['torque']]
        ncol = 4

    if 'meprodmean' in tag:
        nq = 15 # (shear, adv, comp, ind, diff) x (r, th, ph)
        ext = tag[-3:]
        if ext == 'tot':
            iqstart = 0
        if ext == 'mmm':
            iqstart = 1*nq
        if ext == 'mpp':
            iqstart = 2*nq
        qvals = np.arange(iqstart, iqstart + nq)
        titles = []
        for direc in ['r', 'th', 'ph']:
            app = ' (' + direc + ')'
            titles += ['shear' + app, 'comp' + app, 'advec' + app,
                    'induct' + app, 'diff' + app]
        units = [utype['eprod']]
        ncol = 5

    di_out['qvals'] = np.array(qvals)
    di_out['titles'] = np.array(titles)
    if len(units) > 1:
        di_out['units'] = np.array(units)
    else:
        di_out['units'] = units[0]
    di_out['ncol'] = ncol

    return di_out

# basic tex labels
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
