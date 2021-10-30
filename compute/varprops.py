# Author: Loren Matilsky
# Created: 12/19/2018
# Stores properties of various Rayliegh output variables;
# quantity codes, LaTex variable names, etc.
# also groups of quantities, like torques, induction, etc.
# Add to these lists as need be.
# Not all the variables have quantity codes but are derivative quantities
# of other fluid variables
import numpy as np
import sys
from common import array_of_strings, is_an_int
from lut import *

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

def is_basic(varname_full):
    # first split according to arithmetic operations
    varlist = []
    for ele1 in varname_full.split('*'):
        for ele2 in ele1.split('/'):
            for ele3 in ele2.split('+'):
                for ele4 in ele3.split('-'):
                    # strip prime or sph
                    varlist.append(ele4)
    if len(varlist) == 1: # no splitting occurred---basic var
        basic = True
    else:
        basic = False

    # make sure all these basic variables are "kosher"
    kosher = True
    for varname in varlist:
        if varname[-5:] == 'prime': # prime appears at end to modify varname
            varname = varname[:-5]
        elif varname[-3:] == 'sph':
            varname = varname[:-3]
        if not varname in var_indices and not is_an_int(varname):
            kosher = False
    if kosher:
        return basic
    else:
        print ("ERROR! %s is not a valid variable name" %varname_full)
        print ("exiting")
        sys.exit()

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

fancy_op_labels = dict({'+': r'$+$', '-': r'$-$', '*': r'$\times$', '+': r'$\div$'})
simple_op_labels = dict({'+': 'plus', '-': 'minus', '*': 'times', '+': 'div' })

def get_label(varname):
    # first test if varname is valid and if it's a basic variable
    basic = is_basic(varname)

    if basic:
        varname, deriv, primevar, sphvar = get_varprops(varname)
        if deriv:
            derivdir = varname[-1]
            varname = varname[1:] # remove prepending d
            varname = varname[:-2] # remove appending d?

        # now get root label
        # start with thermal vars
        if is_an_int(varname):
            label = varname
        elif varname == 'p':
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
    else:
        operations = []
        for char in varname:
            if char in ['*', '/', '+', '-']:
                operations.append(char)
        label = ''
        simple_label = ''
        count = 0
        for ele1 in varname.split('*'):
            for ele2 in ele1.split('/'):
                for ele3 in ele2.split('+'):
                    for ele4 in ele3.split('-'):
                        label += get_label(ele4)
                        simple_label += ele4
                        if count < len(operations):
                            label += (r' ' + fancy_op_labels[operations[count]] + r' ')
                            simple_label += ('_' + simple_op_labels[operations[count]] + '_')
                        count += 1
        return label, simple_label

# groups of quantities
def get_quantity_group(groupname, magnetism):
    di_out = dict({'groupname': groupname})
    ncol = None
    totsig = None
    qvals = None

    # set default qvals: velocity + vorticity, Pressure/ entropy
    # then possibly B, del x B
    if groupname == 'default':
        qvals = [1, 2, 3, 301, 302, 303, 501, 502]
        if magnetism:
            qvals += [801, 802, 803, 1001, 1004, 1007]

    if groupname == 'torque':
        qvals = [1819, 1801, 1802, 1803, 1804]

        if magnetism:
            qvals += [1805, 1806]
        totsig = np.ones(len(qvals))
        totsig[0] = 0; totsig[1] = totsig[2] = -1

    if groupname == 'induct':
        qvals = [801, 802, 803]            
        for j in range(1, 31):
            qvals.append(1600 + j)

    if groupname == 'v':
        qvals = [1, 2, 3]            
        ncol = 3

    if groupname == 'b':
        qvals = [801, 802, 803]            
        ncol = 3

    if groupname == 'forcer': # linear forces, radial
        qvals = [lookup('rhov_r'), 1201, 1219, 1237, 1216, 1228]
        if magnetism:
            qvals += [1248]
        ncol = 3
        totsig = np.ones(len(qvals))
        totsig[0] = 0
        totsig[1] = -1

    if groupname == 'forcet': # linear forces, theta
        qvals = [lookup('rhov_theta'), 1202, 1220, 1238, 1229]
        if magnetism:
            qvals += [1249]
        ncol = 3
        totsig = np.ones(len(qvals))
        totsig[0] = 0
        totsig[1] = -1

    if groupname == 'forcep': # linear forces, phi
        qvals = [lookup('rhov_phi'), 1203, 1221, 1239, 1230]
        if magnetism:
            qvals += [1250]
        ncol = 3
        totsig = np.ones(len(qvals))
        totsig[0] = 0
        totsig[1] = -1

    if groupname == 'efr': # energy fluxes, r
        qvals = [1455, 1458, 1470, 1935, 1923]
        if magnetism:
            qvals += [2001]
        ncol = 3
        totsig = np.ones(len(qvals))
        totsig[3] = -1

    if groupname == 'eft': # energy fluxes, theta
        qvals = [1456, 1459, 1471, 1936, 1924]
        if magnetism:
            qvals += [2002]
        ncol = 3
        totsig = np.ones(len(qvals))
        totsig[3] = -1

    if groupname == 'efp': # energy fluxes, phi
        qvals = [1457, 1460, 1937, 1925]
        if magnetism:
            qvals += [2003]
        ncol = 3
        totsig = np.ones(len(qvals))
        totsig[3] = -1

    if groupname == 'indr': # induction, r
        qvals = [801, 1601, 1602, 1603, 1604, 1605]
        totsig = np.zeros(len(qvals))
        totsig[-2] = totsig[-1] = 1

    if groupname == 'indrmean': # energy fluxes, r, mean
        qvals = [1616, 1617, 1618, 1619, 1620]
        totsig = np.zeros(len(qvals))
        totsig[-2] = totsig[-1] = 1

    if groupname == 'indt': # energy fluxes, theta
        qvals = [802, 1606, 1607, 1608, 1609, 1610]
        totsig = np.zeros(len(qvals))
        totsig[-2] = totsig[-1] = 1

    if groupname == 'indtmean': # energy fluxes, theta, mean
        qvals = [1621, 1622, 1623, 1624, 1625]
        totsig = np.zeros(len(qvals))
        totsig[-2] = totsig[-1] = 1

    if groupname == 'indp': # energy fluxes, phi
        qvals = [803, 1611, 1612, 1613, 1614, 1615]
        totsig = np.zeros(len(qvals))
        totsig[-2] = totsig[-1] = 1

    if groupname == 'indpmean': # energy fluxes, phi, mean
        qvals = [1626, 1627, 1628, 1629, 1630]
        totsig = np.zeros(len(qvals))
        totsig[-2] = totsig[-1] = 1

    if groupname == 'ke':
        qvals = [402, 403, 404, 410, 411, 412]
        ncol = 3
        totsig = 'sumrow'
    
    if groupname == 'me':
        qvals = [1102, 1103, 1104, 1110, 1111, 1112]
        ncol = 3
        totsig = 'sumrow'

    if groupname == 'teq': # thermal equation
        ncol = 5
        qvals = [1401, lookup('advref'), 1421, 1434, 1435]
        titles = ['(-) adv (pert)', '(-) adv (ref)', 'cond', 'Q(r)', 'visc']
        if magnetism:
            ncol +=1
            qvals.append(1436)
            titles.append('joule')

        totsig = np.ones(ncol)
        totsig[0] = totsig[1] = -1

    # default just use the Rayleigh quantity abbreviations (if qvals has been defined by now)
    if not qvals is None:
        titles = parse_quantities(qvals)[1]

    if 'meprod' in [groupname[:-3], groupname[:-4]]: # this is the exact stuff
        baselen = 6
        ext = groupname[baselen:baselen + 3]
        basetitles = ['ME', 'induct', 'shear', 'advec', 'comp', 'diff']

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
            qvals += [1102 + count] # add in the magnetic energy
            qvals_loc = iqstart + np.arange(0, ncol*3, 3) + count
            qvals += qvals_loc.tolist()
            titles_loc = []
            for j in range(ncol+1):
                titles_loc.append(basetitles[j] + ' (' + ext + ', ' +\
                        direc + ')')
            titles += titles_loc
            count += 1
        ncol += 1 # make room for magnetic energy
        totsig = np.zeros(ncol)
        totsig[0] = 0
        totsig[1] = 1
        if ncol == 6: # include diffusion
            totsig[5] = 1

    if groupname[:9] == 'meprodnum':
        nq = 12 # (r, t, p) x (ind, shear, adv, comp)
        baselen = 9
        ext = groupname[baselen:baselen + 3]
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

    if groupname[:11] == 'meprodtheta':
        nq = 5 # theta (t): (ind, shear, adv, comp, diff)
        baselen = 11
        ext = groupname[baselen:baselen + 3]
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
        direc = 'th'
        app = ' (' + direc + ')'
        titles += ['induct' + app, 'shear' + app, 'advec' + app,\
                    'comp' + app, 'diff' + app]
        ncol = 5
        totsig = np.zeros(len(qvals))
        totsig[0] = totsig[4] = 1

    if groupname[:11] in ['meprodshear', 'meprodadvec']:
        nq = 15 # (r, t, p) x (br (d/dr), bt (d/dt), bp (d/dp), curv1, curv2
        # (similar for adv)
        baselen = 11
        ext = groupname[baselen:baselen + 3]
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
            if 'shear' in groupname:
                titles += ['br (d/dr)' + app, 'bt (d/dT)' + app, 'bp (d/dP)' + app, 'curv1' + app, 'curv2' + app]
            else:
                titles += ['-vr (d/dr)' + app, '-vt (d/dT)' + app, '-vp (d/dP)' + app, 'curv1' + app, 'curv2' + app]

        totsig = np.array([1, 1, 1, 1, 0]) # remove curv2 term by default
                # must correct for this when plotting azimuthal shear
        ncol = 5

    if groupname[:10] == 'meprodmean':
        nq = 15 # (shear, adv, comp, ind, diff) x (r, th, ph)
        baselen = 10
        ext = groupname[baselen:baselen + 3]
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
        totsig = np.array([0, 0, 0, 1, 1])

    if groupname[:6] == 'meprod':
        if len(groupname) == baselen + 4: # another extension, 
            # indicating one direction only
            ext2 = groupname[-1]
            if ext2 == 'r':
                qvals = qvals[:ncol]
                titles = titles[:ncol]
            if ext2 == 't':
                qvals = qvals[ncol:2*ncol]
                titles = titles[ncol:2*ncol]
            if ext2 == 'p':
                qvals = qvals[2*ncol:3*ncol]
                titles = titles[2*ncol:3*ncol]

    if groupname in ['indralt', 'indraltnum']:
        ind_off = 0
        qvals = np.arange(15) + ind_off
        titles = ['(d/dT)(vr*bt)', '-(d/dT)(vt*br)', '-(d/dP)(vp*br)', '(d/dP)(vr*bp)', 'curv.']

    if groupname in ['indtalt', 'indtaltnum']:
        ind_off = 15
        qvals = np.arange(15) + ind_off
        titles = ['(d/dP)(vt*bp)', '-(d/dP)(vp*bt)', '-(d/dr)(vr*bt)', '(d/dr)(vt*br)', 'curv.']

    if groupname in ['indpalt', 'indpaltnum']:
        ind_off = 30
        qvals = np.arange(15) + ind_off
        titles = ['(d/dr)(vp*br)', '-(d/dr)(vr*bp)', '-(d/dT)(vt*bp)', '(d/dT)(vp*bt)', 'curv.']

    if 'ind' in groupname and 'alt' in groupname:
        for ext in ['mm', 'pp']:
            more_titles = []
            for j in range(5):
                more_titles.append(titles[j] + '_' + ext)
            titles += more_titles
        ncol = 5 

    if groupname == 'magtorquemm':
        titles = [r'$\tau_{\rm{mm,r}}$', r'$\frac{r\sin\theta}{4\pi}\left\langle B_r\right\rangle\left\langle\frac{\partial B_\phi}{\partial r}\right\rangle$', r'$\frac{r\sin\theta}{4\pi}\left\langle B_\phi\right\rangle\left\langle\frac{\partial B_r}{\partial r}\right\rangle$', r'$\frac{3\sin\theta}{4\pi}\langle B_r\rangle\langle B_\phi\rangle$',\
        r'$\tau_{\rm{mm,\theta}}$', r'$\frac{\sin\theta}{4\pi}\left\langle B_\theta\right\rangle \left\langle\frac{\partial B_\phi}{\partial \theta}\right\rangle$', r'$\frac{\sin\theta}{4\pi}\left\langle B_\phi\right\rangle\left\langle\frac{\partial B_\theta}{\partial \theta}\right\rangle$', r'$\frac{2\cos\theta}{4\pi}\langle B_\theta\rangle\langle B_\phi\rangle$']
        qvals = [2, 6, 7, 8, 3, 9, 10, 11]
        ncol = 4

    if groupname == 'magtorquems':
        titles = [r'$\tau_{\rm{ms,r}}$', r'$\frac{r\sin\theta}{4\pi}\left\langle B_r^\prime\frac{\partial B_\phi^\prime}{\partial r}\right\rangle$', r'$\frac{r\sin\theta}{4\pi}\left\langle B_\phi^\prime\frac{\partial B_r^\prime}{\partial r}\right\rangle$', r'$\frac{3\sin\theta}{4\pi}\langle B_r^\prime B_\phi^\prime\rangle$',\
        r'$\tau_{\rm{ms,\theta}}$', r'$\frac{\sin\theta}{4\pi}\left\langle B_\theta^\prime\frac{\partial B_\phi^\prime}{\partial \theta}\right\rangle$', r'$\frac{\sin\theta}{4\pi}\left\langle B_\phi^\prime\frac{\partial B_\theta^\prime}{\partial \theta}\right\rangle$', r'$\frac{2\cos\theta}{4\pi}\langle B_\theta^\prime B_\phi^\prime\rangle$']
        qvals = [4, 12, 13, 14, 5, 15, 16, 17]
        ncol = 4

    if groupname == 'ferraro':
        ncol = 5
        qvals = np.arange(15)
        # order: flux_r, flux_t, torque_r, torque_t, torque
        # x tot, mm, pp
        titles = []
        bases = ['flux_r', 'flux_t', 'torque_r', 'torque_t', 'torque']
        suffixes = [' (full)', ' (mm)', ' (pp)']
        for j in range(3):
            for k in range(5):
                titles.append(bases[k] + suffixes[j])

    di_out['qvals'] = np.array(qvals)
    di_out['titles'] = np.array(titles)
    di_out['ncol'] = ncol
    di_out['totsig'] = totsig

    return di_out
