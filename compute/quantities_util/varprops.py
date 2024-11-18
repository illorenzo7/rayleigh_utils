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
from common import array_of_strings, is_an_int, dotdict
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

    'absomr'   :       301,
    'absomt'   :       302,
    'absomp'   :       303,

    'pvr'   :       301,
    'pvt'   :       302,
    'pvp'   :       303,

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
    'jp'    :       1007
}

var_cyl = ['vl', 'vz', 'oml', 'omz', 'absoml', 'absomz', 'pvl', 'pvz', 'bl', 'bz', 'jl', 'jz']

rootlabels = {'v': r'$v$', 'b': r'$B$', 'om': r'$\omega$', 'absom': r'$(\mathbf{\omega} + 2\mathbf{\Omega_0})$', 'pv': r'$Q$', 'j': r'$\mathcal{J}$'}
direclabels = {'r': r'$r$', 't': r'$\theta$', 'p': r'$\phi$', 'l': r'$\lambda$', 'z': r'$z$', 'T': r'$\Theta$', 'P': r'$\Phi$'}

def is_basic(varname_full):
    # first split according to arithmetic operations
    varlist = []
    for ele1 in varname_full.split('*'):
        for ele2 in ele1.split('/'):
            for ele3 in ele2.split('+'):
                for ele4 in ele3.split('='):
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
        if not varname in var_indices and not is_an_int(varname) and not varname in var_cyl:
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

fancy_op_labels = dict({'+': r'$+$', '=': r'$-$', '*': r'$\times$', '/': r'$\div$'})
simple_op_labels = dict({'+': 'plus', '=': 'minus', '*': 'times', '/': 'div' })

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
            if char in ['*', '/', '+', '=']:
                operations.append(char)
        label = ''
        simple_label = ''
        count = 0
        for ele1 in varname.split('*'):
            for ele2 in ele1.split('/'):
                for ele3 in ele2.split('+'):
                    for ele4 in ele3.split('='):
                        label += get_label(ele4)
                        simple_label += ele4
                        if count < len(operations):
                            label += (r' ' + fancy_op_labels[operations[count]] + r' ')
                            simple_label += ('_' + simple_op_labels[operations[count]] + '_')
                        count += 1
        return label, simple_label

# groups of quantities
def get_quantity_group(groupname, magnetism):
    di_out = dotdict({'groupname': groupname})
    ncol = None
    totsig = None
    qvals = None
    titles = None

    # set default qvals: velocity + vorticity, Pressure/ entropy
    # then possibly B, del x B
    if groupname == 'default':
        qvals = [1, 2, 3, 301, 302, 303, 501, 502]
        if magnetism:
            qvals += [801, 802, 803, 1001, 1004, 1007]

    if groupname == 'thermal':
        qvals = [501,502, 507,508, 513,514]

    if groupname == 'torque':
        qvals = [1819, 1801, 1802, 1803, 1804]
        titles = ['L_z', 'tau_rs', 'tau_mcdr', 'tau_mc0', 'tau_v']

        if magnetism:
            qvals += [1805, 1806]
            titles += ['tau_mm', 'tau_ms']
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
        qvals = [1201, 1219, 1237, 1216, 1228]
        if magnetism:
            qvals += [1248]
        totsig = np.ones(len(qvals))
        totsig[0] = -1

    if groupname == 'forcet': # linear forces, theta
        qvals = [1202, 1220, 1238, 1229]
        if magnetism:
            qvals += [1249]
        totsig = np.ones(len(qvals))
        totsig[0] = -1

    if groupname == 'forcep': # linear forces, phi
        qvals = [1203, 1221, 1239, 1230]
        if magnetism:
            qvals += [1250]
        totsig = np.ones(len(qvals))
        totsig[0] = -1

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

    if groupname == 'indraltdiff': # induction, r, old diffusion
        qvals = [801, 1601, 1602, 1603, 1604, 2911]
        totsig = np.zeros(len(qvals))
        totsig[-2] = totsig[-1] = 1

    if groupname == 'indrmean': # induction, r, mean
        qvals = [1616, 1617, 1618, 1619, 1620]
        totsig = np.zeros(len(qvals))
        totsig[-2] = totsig[-1] = 1

    if groupname == 'indt': # induction, theta
        qvals = [802, 1606, 1607, 1608, 1609, 1610]
        totsig = np.zeros(len(qvals))
        totsig[-2] = totsig[-1] = 1

    if groupname == 'indtaltdiff': # induction, theta, old diffusion
        qvals = [802, 1606, 1607, 1608, 1609, 2912]
        totsig = np.zeros(len(qvals))
        totsig[-2] = totsig[-1] = 1

    if groupname == 'indtmean': # induction, theta, mean
        qvals = [1621, 1622, 1623, 1624, 1625]
        totsig = np.zeros(len(qvals))
        totsig[-2] = totsig[-1] = 1

    if groupname == 'indp': # induction, phi
        qvals = [803, 1611, 1612, 1613, 1614, 1615]
        totsig = np.zeros(len(qvals))
        totsig[-2] = totsig[-1] = 1

    if groupname == 'indpaltdiff': # induction, phi, old diffusion
        qvals = [803, 1611, 1612, 1613, 1614, 2913]
        totsig = np.zeros(len(qvals))
        totsig[-2] = totsig[-1] = 1
   
    if groupname == 'indpmean': # energy fluxes, phi, mean
        qvals = [1626, 1627, 1628, 1629, 1630]
        totsig = np.zeros(len(qvals))
        totsig[-2] = totsig[-1] = 1

    if groupname == 'indr29': # induction, r 29
        qvals = [801, 1601, 1602, 1603, 1604, 2911]
        totsig = np.zeros(len(qvals))
        totsig[-2] = totsig[-1] = 1

    if groupname == 'indt29': # induction, theta
        qvals = [802, 1606, 1607, 1608, 1609, 2912]
        totsig = np.zeros(len(qvals))
        totsig[-2] = totsig[-1] = 1

    if groupname == 'indp29': # induction, phi
        qvals = [803, 1611, 1612, 1613, 1614, 2913]
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
        #qvals = [701, 1401, 1402, 1421, 1434, 1435]
        qvals = [701, 1401, 1402, 1421, 1435]
        titles = ['rho*T*S', 'adv (tot)', '-adv (fluc)', 'cond', 'visc']
        #titles = ['rho*T*S', 'adv (tot)', '-adv (fluc)', 'cond', 'Q(r)', 'visc']
        if magnetism:
            ncol +=1
            qvals.append(1436)
            titles.append('joule')

        totsig = np.ones(ncol)
        totsig[0] = totsig[2] = 0.
        totsig[1] = -1

    # default just use the Rayleigh quantity abbreviations (if qvals has been defined by now)
    if titles is None and not qvals is None:
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

    if 'indmeprod' in [groupname, groupname[:-1]]: # sum over the different components of the production by induction
        baselen = 9
        basetitles = ['tot', 'pmp', 'ppm', 'mmm', 'mpp', 'ppp']
        ncol = 6

        # do r, theta, then phi, production terms, in that order
        qvals = []
        titles = []
        count = 0
        for direc in ['r', 'th', 'ph']:
            qvals += [2201 + count] # tot. ind.
            qvals += [2216 + count] # pmp
            qvals += [2228 + count] # ppm
            qvals += [2240 + count] # mmm
            qvals += [2255 + count] # mpp
            qvals += [2264 + count] # ppp

            titles_loc = []
            for j in range(ncol):
                titles_loc.append(basetitles[j] + ' (' + direc + ')')
            titles += titles_loc
            count += 1
        totsig = np.ones(ncol)
        totsig[0] = 0

    if 'diffmeprod' in [groupname, groupname[:-1]]: # sum over the different components of the production by induction
        baselen = 10
        basetitles = ['tot', 'mm', 'pp']
        ncol = 3

        # do r, theta, then phi, production terms, in that order
        qvals = []
        titles = []
        count = 0
        for direc in ['r', 'th', 'ph']:
            qvals += [2213 + count] # tot. diff.
            qvals += [2252 + count] # mm
            qvals += [2276 + count] # pp

            titles_loc = []
            for j in range(ncol):
                titles_loc.append(basetitles[j] + ' (' + direc + ')')
            titles += titles_loc
            count += 1
        totsig = np.ones(ncol)
        totsig[0] = 0

    if 'meprodalt' in [groupname[:-3], groupname[:-4]]:
        # exact forms of "alternate" induction terms
        baselen = 9
        ext = groupname[baselen:baselen + 3]

        basetitles = ['ME', 'induct', 'shear', 'advec', 'comp', 'diff']

        ncol = 5

        # do r, theta, then phi, production terms, in that order
        qvals = []
        titles = []
        count = 0
        iqstart = 2301 # these are the "alternate" terms
        for direc in ['r', 'th', 'ph']:
            qvals += [1102 + count] # add in the magnetic energy
            qvals_loc = iqstart + np.arange(0, 12, 3) + count
            qvals_loc = qvals_loc.tolist()
            qvals_loc += [2213 + count] # diffusive terms haven't changed
            qvals += qvals_loc
            titles_loc = []
            for j in range(ncol+1):
                titles_loc.append(basetitles[j] + ' (tot, ' +\
                        direc + ')')
            titles += titles_loc
            count += 1
        ncol += 1 # make room for magnetic energy
        totsig = np.zeros(ncol)
        totsig[0] = 0
        totsig[1] = totsig[5] = 1

    if 'shearmeprodalt' in [groupname, groupname[:-1]]: # sum over the different components of the production by shear production
        ncol = 3

        # do r, theta, then phi, production terms, in that order
        qvals = []
        titles = []
        count = 0
        for direc in ['r', 'th', 'ph']:
            if direc == 'r':
                perp1 = 'th'
                perp2 = 'ph'
            if direc == 'th':
                perp1 = 'r'
                perp2 = 'ph'
            if direc == 'ph':
                perp1 = 'r'
                perp2 = 'th'

            qvals += [2304 + count] # tot. shear
            qvals += [2401 + count] # shear1
            qvals += [2404 + count] # shear2

            titles += ['tot unbroken',\
                    'B_' + perp1 + ' dotgrad v_' + direc,\
                    'B_' + perp2 + ' dotgrad v_' + direc]
            count += 1
        totsig = np.ones(ncol)
        totsig[0] = 0

    if 'advmeprodalt' in [groupname, groupname[:-1]]: # sum over the different components of the production by shear production
        ncol = 4
        # do r, theta, then phi, production terms, in that order
        qvals = []
        titles = []
        count = 0
        for direc in ['r', 'th', 'ph']:
            app = 'B_' + direc
            titles += ['tot unbroken', '-vr (d/dr)' + app, '-vt (d/dT)' + app, '-vp (d/dP)' + app]

            qvals += [2307 + count] # tot. adv
            qvals += [2407 + count] 
            qvals += [2410 + count]
            qvals += [2413 + count] 
            count += 1
        totsig = np.ones(ncol)
        totsig[0] = 0

    if 'compmeprodalt' in [groupname, groupname[:-1]]: # sum over the different components of the production by shear production
        ncol = 3

        # do r, theta, then phi, production terms, in that order
        qvals = []
        titles = []
        count = 0
        for direc in ['r', 'th', 'ph']:
            if direc == 'r':
                perp1 = 'th'
                perp2 = 'ph'
            if direc == 'th':
                perp1 = 'r'
                perp2 = 'ph'
            if direc == 'ph':
                perp1 = 'r'
                perp2 = 'th'

            qvals += [2310 + count] # tot. comp.
            qvals += [2416 + count] # comp1
            qvals += [2419 + count] # comp2

            titles += ['tot unbroken',\
                    '-B_' + direc + ' div v_' + perp1,\
                    '-B_' + direc + ' div v_' + perp2]

            count += 1
        totsig = np.ones(ncol)
        totsig[0] = 0

    di_out['qvals'] = np.array(qvals)
    di_out['titles'] = np.array(titles)
    di_out['ncol'] = ncol
    di_out['totsig'] = totsig

    return di_out
