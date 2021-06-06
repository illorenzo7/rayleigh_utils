# Author: Loren Matilsky
# Created: 12/19/2018
# Stores properties of various Rayliegh output variables;
# quantity codes, LaTex variable names, etc.
# also groups of quantities, like torques, induction, etc.
# Add to these lists as need be.
# Not all the variables have quantity codes but are derivative quantities
# of other fluid variables
import numpy as np
from common import array_of_strings

# basic variable indices
var_indices = {\
    'vr'    :       1, 
    'vt'    :       2,
    'vp'    :       3,              

    'dvrdr' :       10,
    'dvtdr' :       11,
    'dvpdr' :       12,

    'dvrdt' :       19,
    'dvtdt' :       20,
    'dvpdt' :       21,

    'dvrdp' :       28,
    'dvtdp' :       29,
    'dvpdp' :       30,

    'dvrdT' :       37,
    'dvtdT' :       38,
    'dvpdT' :       39,

    'dvrdP' :       46,
    'dvtdP' :       47,
    'dvpdP' :       48,

    'omr'   :       301,
    'omt'   :       302,
    'omp'   :       303,

    's'     :       501,
    'p'     :       502,

    'br'    :       801,
    'bt'    :       802,
    'bp'    :       803,

    'dbrdr' :       810,
    'dbtdr' :       811,
    'dbpdr' :       812,

    'dbrdt' :       819,
    'dbtdt' :       820,
    'dbpdt' :       821,

    'dbrdp' :       828,
    'dbtdp' :       829,
    'dbpdp' :       830,

    'dbrdT' :       837,
    'dbtdT' :       838,
    'dbpdT' :       839,

    'dbrdP' :       846,
    'dbtdP' :       847,
    'dbpdP' :       848,

    'jr'    :       1001,
    'jt'    :       1004,
    'jp'    :       1007}

rootlabels = {'v': r'$v$', 'b': r'$B$', 'om': r'$\omega$', 'j': r'$\mathcal{J}$'}
direclabels = {'r': r'$r$', 't': r'$\theta$', 'p': r'$\phi$', 'l': r'$\lambda$', 'z': r'$z$', 'T': r'$\Theta$', 'P': r'$\Phi$'}

def get_varprops(varname):
    # get boolean values deriv, prime (az average subtracted), and sph
    # (sph average subtracted)
    deriv = False
    primevar = False
    sphvar = False
    if 'dr' in varname or 'dt' in varname or 'dp' in varname or 'dT' in varname or 'dP' in varname:
        deriv = True
    if varname[-5:] == 'prime': # prime appears at end to modify varname
        primevar = True
        varname = varname[:-5]
    elif varname[-3:] == 'sph':
        sphvar = True
        varname = varname[:-3]
    return varname, deriv, primevar, sphvar

def get_label(varname):
    varname, deriv, primevar, sphvar = get_varprops(varname)
    if deriv:
        derivdir = varname[-1]
        varname = varname[1:] # remove prepending d
        varname = varname[:-2] # remove appending d?

    # now get root label
    # start with thermal vars
    if varname == 'p':
        label = r'$P$'
    elif varname == 's':
        label = r'$S$'
    # now field variables, with a vector component
    if 'v' in varname or 'b' in varname or 'om' in varname or 'j' in varname:
        rootname = varname[:-1]
        direction = varname[-1]
        label = rootlabels[rootname] + r'$_$' + direclabels[direction]
    if primevar:
        label += r'$^\prime$'
    elif sphvar:
        label += r'$^{\prime\prime}$'
    if deriv:
        label = r'${\partial}$' + label + r'$/$' + r'${\partial}$' + direclabels[derivdir]
    return label.replace('$$', '')

# groups of quantities
def get_quantity_group(tag, magnetism):
    di_out = dict({})
    di_out['totsig'] = None # the default

    # set default qvals: velocity + vorticity, Pressure/ entropy
    # then possibly B, del x B
    if tag == 'default':
        qvals = [1, 2, 3, 301, 302, 303, 501, 502]
        titles = [r'$v_r$', r'$v_\theta$', r'$v_\phi$',\
                r'$\omega_r$', r'$\omega_\theta$', r'$\omega_\phi$',\
                r'$S$', r'$P$']

        if magnetism:
            qvals += [801, 802, 803, 1001, 1002, 1003]
            titles += [r'$B_r$', r'$B_\theta$', r'$B_\phi$',\
                r'$\mathcal{J}_r$', r'$\mathcal{J}_\theta$', r'$\mathcal{J}_\phi$']
        ncol = 4 + 3*magnetism

    if tag == 'torque':
        qvals = [3, 1801, 1802, 1803, 1804, 1819]
        titles = [r'$v_\phi$', r'$-\tau_{\rm{rs}}$', r'$-\tau_{\rm{mc,v_\phi}}$',r'$\tau_{\rm{mc,\Omega_0}}$', r'$\tau_{\rm{v}}$',\
          r'$\mathcal{L}_z$']

        if magnetism:
            qvals += [1805, 1806]
            titles += [r'$\tau_{\rm{mm}}$', r'$\tau_{\rm{ms}}$']
        ncol = 3 + magnetism
        totsig = np.ones(len(titles))
        totsig[0] = totsig[5] = 0; totsig[1] = totsig[2] = -1
        di_out['totsig'] = totsig

    if tag == 'induct':
        qvals = [801, 802, 803]            
        for j in range(1, 31):
            qvals.append(1600 + j)
        titles = array_of_strings(qvals)
        ncol = 11

    if tag == 'v':
        qvals = [1, 2, 3]            
        titles = [get_label('vr'), get_label('vt'), get_label('vp')]
        ncol = 3

    if tag == 'b':
        qvals = [801, 802, 803]            
        titles = [get_label('br'), get_label('bt'), get_label('bp')]
        ncol = 3

    if tag == 'forcer': # linear forces, radial
        qvals = [1201, 1219, 1237, 1216, 1228]
        titles= [r'$-(\mathbf{f}_{\rm{adv}})_r$', r'$(\mathbf{f}_{\rm{cor}})_r$',  r'$(\mathbf{f}_{\rm{p}})_r$', r'$(\mathbf{f}_{\rm{buoy}})_r$', r'$(\mathbf{f}_{\rm{v}})_r$']
        if magnetism:
            qvals += [1248]
            titles += [r'$(\mathbf{f}_{\rm{mag}})_r$']
        ncol = 3
        totsig = np.ones(len(titles))
        totsig[0] = -1
        di_out['totsig'] = totsig

    if tag == 'forcet': # linear forces, theta
        qvals = [1202, 1220, 1238, 1229]
        titles= [r'$-(\mathbf{f}_{\rm{adv}})_\theta$', r'$(\mathbf{f}_{\rm{cor}})_\theta$',  r'$(\mathbf{f}_{\rm{p}})_\theta$',  r'$(\mathbf{f}_{\rm{v}})_\theta$']
        if magnetism:
            qvals += [1249]
            titles += [r'$(\mathbf{f}_{\rm{mag}})_\theta$']
        ncol = 3
        totsig = np.ones(len(titles))
        totsig[0] = -1
        di_out['totsig'] = totsig

    if tag == 'forcep': # linear forces, phi
        qvals = [1203, 1221, 1239, 1230]
        titles= [r'$-(\mathbf{f}_{\rm{adv}})_\phi$', r'$(\mathbf{f}_{\rm{cor}})_\phi$',  r'$(\mathbf{f}_{\rm{p}})_\phi$',  r'$(\mathbf{f}_{\rm{v}})_\phi$']
        if magnetism:
            qvals += [1250]
            titles += [r'$(\mathbf{f}_{\rm{mag}})_\phi$']
        ncol = 3
        totsig = np.ones(len(titles))
        totsig[0] = -1
        di_out['totsig'] = totsig

    if tag == 'efr': # energy fluxes, r
        qvals = [1455, 1458, 1470, 1935, 1923]
        titles = [r'$(\mathbf{\mathcal{F}}_{\rm{enth}})_r$',\
          r'$(\mathbf{\mathcal{F}}_{\rm{enth,pp}})_r$',\
          r'$(\mathbf{\mathcal{F}}_{\rm{cond}})_r$',\
          r'$-(\mathbf{\mathcal{F}}_{\rm{visc}})_r$',\
          r'$(\mathbf{\mathcal{F}}_{\rm{KE}})_r$']
        if magnetism:
            qvals += [2001]
            titles += [r'$(\mathbf{\mathcal{F}}_{\rm{Poynt}})_r$']
        ncol = 3
        totsig = np.ones(len(titles))
        totsig[3] = -1
        di_out['totsig'] = totsig

    if tag == 'eft': # energy fluxes, theta
        qvals = [1456, 1459, 1471, 1936, 1924]
        titles = [r'$(\mathbf{\mathcal{F}}_{\rm{enth}})_\theta$',\
          r'$(\mathbf{\mathcal{F}}_{\rm{enth,pp}})_\theta$',\
          r'$(\mathbf{\mathcal{F}}_{\rm{cond}})_\theta$',\
          r'$-(\mathbf{\mathcal{F}}_{\rm{visc}})_\theta$',\
          r'$(\mathbf{\mathcal{F}}_{\rm{KE}})_\theta$']
        if magnetism:
            qvals += [2002]
            titles += [r'$(\mathbf{\mathcal{F}}_{\rm{Poynt}})_\theta$']
        ncol = 3
        di_out['totsig'] = 'sumrow'
        totsig = np.ones(len(titles))
        totsig[3] = -1
        di_out['totsig'] = totsig

    if tag == 'efp': # energy fluxes, phi
        qvals = [1457, 1460, 1937, 1925]
        titles = [r'$(\mathbf{\mathcal{F}}_{\rm{enth}})_\phi$',\
          r'$(\mathbf{\mathcal{F}}_{\rm{enth,pp}})_\phi$',\
          r'$-(\mathbf{\mathcal{F}}_{\rm{visc}})_\phi$',\
          r'$(\mathbf{\mathcal{F}}_{\rm{KE}})_\phi$']
        if magnetism:
            qvals += [2003]
            titles += [r'$(\mathbf{\mathcal{F}}_{\rm{Poynt}})_\phi$']
        ncol = 3
        totsig = np.ones(len(titles))
        totsig[3] = -1
        di_out['totsig'] = totsig

    if tag == 'indr': # induction, r
        qvals = [1601, 1602, 1603, 1604, 1605, 801]
        titles = [r'$[\left\langle\mathbf{B}\cdot\nabla\mathbf{v}\right\rangle]_r$',\
            r'$-\left\langleB_r(\nabla\cdot\mathbf{v})\right\rangle$',\
            r'$-[\left\langle\mathbf{v}\cdot\nabla\mathbf{B}\right\rangle]_r$',\
            r'$[\nabla\times(\left\langle\mathbf{v}\times\mathbf{B}\right\rangle)]_r$',\
            r'$-[\nabla\times(\eta\nabla\times\langle\mathbf{B}\rangle)]_r$',\
            r'$\left\langle B_r\right\rangle$']
        ncol = 3
        totsig = np.zeros(len(titles))
        totsig[3] = tosig[4] = 1
        di_out['totsig'] = totsig

    if tag == 'indrmean': # energy fluxes, r, mean
        qvals = [1616, 1617, 1618, 1619, 1620, 801]
        titles = [r'$[\left\langle\mathbf{B}\right\rangle\cdot\nabla\left\langle\mathbf{v}\right\rangle]_r$',\
    r'$-\left\langleB_r\right\rangle(\nabla\cdot\left\langle\mathbf{v}\right\rangle)$',\
    r'$-[\left\langle\mathbf{v}\right\rangle\cdot\nabla\left\langle\mathbf{B}\right\rangle]_r$',\
    r'$[\nabla\times(\left\langle\mathbf{v}\right\rangle\times\left\langle\mathbf{B}\right\rangle)]_r$',\
    r'$-[\nabla\times(\eta\nabla\times\langle\mathbf{B}\rangle)]_r$',\
    r'$\left\langleB_r\right\rangle$']
        ncol = 3
        totsig = np.zeros(len(titles))
        totsig[3] = tosig[4] = 1
        di_out['totsig'] = totsig

    if tag == 'indt': # energy fluxes, theta
        qvals = [1606, 1607, 1608, 1609, 1610, 802]
        titles = [r'$[\left\langle\mathbf{B}\cdot\nabla\mathbf{v}\right\rangle]_\theta$',\
            r'$-\left\langleB_\theta(\nabla\cdot\mathbf{v})\right\rangle$',\
            r'$-[\left\langle\mathbf{v}\cdot\nabla\mathbf{B}\right\rangle]_\theta$',\
            r'$[\nabla\times(\left\langle\mathbf{v}\times\mathbf{B}\right\rangle)]_\theta$',\
            r'$-[\nabla\times(\eta\nabla\times\langle\mathbf{B}\rangle)]_\theta$',\
            r'$\left\langle B_\theta\right\rangle$']
        ncol = 3
        totsig = np.zeros(len(titles))
        totsig[3] = tosig[4] = 1
        di_out['totsig'] = totsig

    if tag == 'indtmean': # energy fluxes, theta, mean
        qvals = [1621, 1622, 1623, 1624, 1625, 802]
        titles = [r'$[\left\langle\mathbf{B}\right\rangle\cdot\nabla\left\langle\mathbf{v}\right\rangle]_\theta$',\
            r'$-\left\langleB_\theta\right\rangle(\nabla\cdot\left\langle\mathbf{v}\right\rangle)$',\
            r'$-[\left\langle\mathbf{v}\right\rangle\cdot\nabla\left\langle\mathbf{B}\right\rangle]_\theta$',\
            r'$[\nabla\times(\left\langle\mathbf{v}\right\rangle\times\left\langle\mathbf{B}\right\rangle)]_\theta$',\
            r'$-[\nabla\times(\eta\nabla\times\langle\mathbf{B}\rangle)]_\theta$',\
            r'$\left\langleB_\theta\right\rangle$']
        ncol = 3
        totsig = np.zeros(len(titles))
        totsig[3] = tosig[4] = 1
        di_out['totsig'] = totsig

    if tag == 'indp': # energy fluxes, phi
        qvals = [1611, 1612, 1613, 1614, 1615, 803]

        titles = [r'$[\left\langle\mathbf{B}\cdot\nabla\mathbf{v}\right\rangle]_\phi$',\
        r'$-\left\langleB_\phi(\nabla\cdot\mathbf{v})\right\rangle$',\
        r'$-[\left\langle\mathbf{v}\cdot\nabla\mathbf{B}\right\rangle]_\phi$',\
        r'$[\nabla\times(\left\langle\mathbf{v}\times\mathbf{B}\right\rangle)]_\phi$',\
        r'$-[\nabla\times(\eta\nabla\times\langle\mathbf{B}\rangle)]_\phi$',\
        r'$\left\langle B_\phi\right\rangle$']
        ncol = 3
        totsig = np.zeros(len(titles))
        totsig[3] = tosig[4] = 1
        di_out['totsig'] = totsig

    if tag == 'indpmean': # energy fluxes, phi, mean
        qvals = [1626, 1627, 1628, 1629, 1630, 803]
        titles = [r'$[\left\langle\mathbf{B}\right\rangle\cdot\nabla\left\langle\mathbf{v}\right\rangle]_\phi$',\
            r'$-\left\langleB_\phi\right\rangle(\nabla\cdot\left\langle\mathbf{v}\right\rangle)$',\
            r'$-[\left\langle\mathbf{v}\right\rangle\cdot\nabla\left\langle\mathbf{B}\right\rangle]_\phi$',\
            r'$[\nabla\times(\left\langle\mathbf{v}\right\rangle\times\left\langle\mathbf{B}\right\rangle)]_\phi$',\
            r'$-[\nabla\times(\eta\nabla\times\langle\mathbf{B}\rangle)]_\phi$',\
            r'$\left\langleB_\phi\right\rangle$']
        ncol = 3
        totsig = np.zeros(len(titles))
        totsig[3] = tosig[4] = 1
        di_out['totsig'] = totsig

    if tag == 'ke':
        qvals = [402, 403, 404, 410, 411, 412]
        titles =\
            [r'$\overline{\rm{KE}}_{\rm{r}}$',\
            r'$\overline{\rm{KE}}_{\rm{\theta}}$',\
            r'$\overline{\rm{KE}}_{\rm{\phi}}$',\
            r'$\rm{KE}^\prime_r$',\
            r'$\rm{KE}^\prime_\theta$',\
            r'$\rm{KE}^\prime_\phi$']
        ncol = 3
        totsig = 'sumrow'
    
    if tag == 'me':
        qvals = [1102, 1103, 1104, 1110, 1111, 1112]
        titles =\
            [r'$\overline{\rm{ME}}_{\rm{r}}$',\
            r'$\overline{\rm{ME}}_{\rm{\theta}}$',\
            r'$\overline{\rm{ME}}_{\rm{\phi}}$',\
            r'$\rm{ME}^\prime_r$',\
            r'$\rm{ME}^\prime_\theta$',\
            r'$\rm{ME}^\prime_\phi']
        ncol = 3
        totsig = 'sumrow'

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
        ncol = 4

    if 'meprodshear' in tag:
        nq = 15 # (r, t, p) x (br (d/dr), bt (d/dt), bp (d/dp), curv1, curv2
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
            titles += ['br (d/dr)' + app, 'bt (d/dT)' + app, 'bp (d/dP)' + app, 'curv1' + app, 'curv2' + app]
        ncol = 5
        di_out['totsig'] = 'sumrow'

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
            for j in range(ncol):
                titles_loc.append(basetitles[j] + ' (' + ext + ', ' +\
                        direc + ')')
            titles += titles_loc
            count += 1
        totsig = np.zeros(ncol)
        totsig[0] = 1
        if ncol == 5: # include diffusion
            totsig[4] = 1
        di_out['totsig'] = totsig

    if tag in ['indralt', 'indraltnum']:
        ind_off = 0
        qvals = np.arange(15) + ind_off
        titles = ['(d/dT)(vr*bt)', '-(d/dT)(vt*br)', '-(d/dP)(vp*br)', '(d/dP)(vr*bp)', 'curv.']

    if tag in ['indtalt', 'indtaltnum']:
        ind_off = 15
        qvals = np.arange(15) + ind_off
        titles = ['(d/dP)(vt*bp)', '-(d/dP)(vp*bt)', '-(d/dr)(vr*bt)', '(d/dr)(vt*br)', 'curv.']

    if tag in ['indpalt', 'indpaltnum']:
        ind_off = 30
        qvals = np.arange(15) + ind_off
        titles = ['(d/dr)(vp*br)', '-(d/dr)(vr*bp)', '-(d/dT)(vt*bp)', '(d/dT)(vp*bt)', 'curv.']

    if 'ind' in tag and 'alt' in tag:
        for ext in ['mm', 'pp']:
            more_titles = []
            for j in range(5):
                more_titles.append(titles[j] + '_' + ext)
            titles += more_titles
        ncol = 5 
        di_out['totsig'] = 'sumrow'

    if tag == 'magtorquemm':
        titles = [r'$\tau_{\rm{mm,r}}$', r'$\frac{r\sin\theta}{4\pi}\left\langle B_r\right\rangle\left\langle\frac{\partial B_\phi}{\partial r}\right\rangle$', r'$\frac{r\sin\theta}{4\pi}\left\langle B_\phi\right\rangle\left\langle\frac{\partial B_r}{\partial r}\right\rangle$', r'$\frac{3\sin\theta}{4\pi}\langle B_r\rangle\langle B_\phi\rangle$',\
        r'$\tau_{\rm{mm,\theta}}$', r'$\frac{\sin\theta}{4\pi}\left\langle B_\theta\right\rangle \left\langle\frac{\partial B_\phi}{\partial \theta}\right\rangle$', r'$\frac{\sin\theta}{4\pi}\left\langle B_\phi\right\rangle\left\langle\frac{\partial B_\theta}{\partial \theta}\right\rangle$', r'$\frac{2\cos\theta}{4\pi}\langle B_\theta\rangle\langle B_\phi\rangle$']
        qvals = [2, 6, 7, 8, 3, 9, 10, 11]
        ncol = 4

    if tag == 'magtorquems':
        titles = [r'$\tau_{\rm{ms,r}}$', r'$\frac{r\sin\theta}{4\pi}\left\langle B_r^\prime\frac{\partial B_\phi^\prime}{\partial r}\right\rangle$', r'$\frac{r\sin\theta}{4\pi}\left\langle B_\phi^\prime\frac{\partial B_r^\prime}{\partial r}\right\rangle$', r'$\frac{3\sin\theta}{4\pi}\langle B_r^\prime B_\phi^\prime\rangle$',\
        r'$\tau_{\rm{ms,\theta}}$', r'$\frac{\sin\theta}{4\pi}\left\langle B_\theta^\prime\frac{\partial B_\phi^\prime}{\partial \theta}\right\rangle$', r'$\frac{\sin\theta}{4\pi}\left\langle B_\phi^\prime\frac{\partial B_\theta^\prime}{\partial \theta}\right\rangle$', r'$\frac{2\cos\theta}{4\pi}\langle B_\theta^\prime B_\phi^\prime\rangle$']
        qvals = [4, 12, 13, 14, 5, 15, 16, 17]
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
        ncol = 5
        di_out['totsig'] = np.array([1, 0, 0, 0, 1])

    di_out['qvals'] = np.array(qvals)
    di_out['titles'] = np.array(titles)
    di_out['ncol'] = ncol

    return di_out
