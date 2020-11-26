# Author: Loren Matilsky
# Created: 02/02/2020
#
# Purpose: read in equation_coefficients (or reference/transport pair)
# into single class with all the equation coefficients, having 
# human-readable attributes:
# radius
# density (rho)
# dlnrho
# d2lnrho
# temperature (T)
# dlnT
# pressure (P)
# gravity (g)
# dsdr
# heating (Q)
# nu
# dlnu
# kappa
# dlnkappa
# eta (if magnetism = True)
# dlneta (if magnetism = True)

import numpy as np
import sys, os

sys.path.append(os.environ['rapp'])
sys.path.append(os.environ['raco'])

from reference_tools import equation_coefficients
from rayleigh_diagnostics import ReferenceState, TransportCoeffs
from common import c_P, thermo_R
from get_parameter import get_parameter

class eq_human_readable:
    """Rayleigh Universal Equation Coefficients Structure
    ----------------------------------
    self.nr          : number of radial points
    self.radius      : radial coordinates
    self.density     : density
    self.rho         : density
    self.dlnrho      : logarithmic derivative of density
    self.d2lnrho     : d_by_dr of dlnrho
    self.temperature : temperature
    self.T           : temperature
    self.dlnT        : logarithmic derivative of temperature
    self.pressure    : pressure (rho*R*T)
    self.P           : pressure (rho*R*T)
    self.dsdr        : radial entropy gradient
    self.gravity     : gravity 
    self.heating     : volumetric heating (Q) 
    self.nu          : momentum diffusivity (kinematic viscosity)
    self.dlnu        : logarithmic derivative of the viscosity
    self.kappa       : temperature diffusivity (thermometric conductivity)
    self.dlkappa     : logarithmic derivative of the temp. diffusivity
    self.eta :       : magnetic diffusivity 
    self.dlneta      : logarithmic derivative of magnetic diffusivity
    self.lum         : (scalar) stellar luminosity
    """
    def __init__(self, nr):
        self.nr = nr
        self.density = np.zeros(nr)
        self.rho = np.zeros(nr) # same as density
        self.dlnrho = np.zeros(nr)
        self.d2lnrho = np.zeros(nr)
        self.temperature = np.zeros(nr)
        self.T = np.zeros(nr) # same as temperature
        self.dlnT = np.zeros(nr)
        self.pressure = np.zeros(nr)
        self.P = np.zeros(nr)
        self.gravity = np.zeros(nr) 
        self.g = np.zeros(nr) # same as gravity
        self.dsdr = np.zeros(nr)
        self.heating = np.zeros(nr)
        self.Q = np.zeros(nr) # same as heating
        self.nu = np.zeros(nr)
        self.dlnu = np.zeros(nr)
        self.kappa = np.zeros(nr)
        self.dlnkappa = np.zeros(nr)
        self.eta = np.zeros(nr) # these should stay zero 
        self.dlneta = np.zeros(nr) # if magnetism = False
        self.lum = 0.0 # luminosity

def get_eq(dirname, fname='equation_coefficients'): # return an eq_human_readable class associated with
    # [dirname], either using equation_coefficients or 
    # transport/reference files
    if os.path.exists(dirname + '/' + fname):
        # by default, get info from equation_coefficients (if file exists)
        eq = equation_coefficients()
        eq.read(dirname + '/' + fname)
        eq_hr = eq_human_readable(eq.nr)

        eq_hr.radius = eq.radius
        eq_hr.density = eq.functions[0]
        eq_hr.rho = eq_hr.density
        eq_hr.dlnrho = eq.functions[7]
        eq_hr.d2lnrho = eq.functions[8]
        eq_hr.temperature = eq.functions[3]
        eq_hr.T = eq_hr.temperature
        eq_hr.pressure = thermo_R*eq_hr.rho*eq_hr.T
        eq_hr.P = eq_hr.pressure
        eq_hr.dlnT = eq.functions[9]
        eq_hr.gravity = eq.functions[1]/eq_hr.rho*c_P
        eq_hr.g = eq_hr.gravity
        eq_hr.dsdr = eq.functions[13]
        eq_hr.heating = eq.constants[9]*eq.functions[5]
        eq_hr.Q = eq_hr.heating
        eq_hr.nu = eq.constants[4]*eq.functions[2]
        eq_hr.dlnu = eq.functions[10]
        eq_hr.kappa = eq.constants[5]*eq.functions[4]
        eq_hr.dlnkappa = eq.functions[11]
        eq_hr.eta = eq.constants[6]*eq.functions[6] # these are built-in to
        eq_hr.dlneta = eq.functions[12] # equation_coefficients as "zero"
        eq_hr.lum = eq.constants[9]
        # if magnetism = False
    else:
        ref = ReferenceState(dirname + '/reference')
        eq_hr = eq_human_readable(ref.nr)

        eq_hr.radius = ref.radius
        eq_hr.density = ref.density
        eq_hr.rho = eq_hr.density
        eq_hr.dlnrho = ref.dlnrho
        eq_hr.d2lnrho = ref.d2lnrho
        eq_hr.temperature = ref.temperature
        eq_hr.T = eq_hr.temperature
        eq_hr.dlnT = ref.dlnt
        eq_hr.pressure = thermo_R*eq_hr.rho*eq_hr.T
        eq_hr.P = eq_hr.pressure
        eq_hr.gravity = ref.gravity
        eq_hr.g = eq_hr.gravity
        eq_hr.dsdr = ref.dsdr
        eq_hr.heating = eq_hr.rho*eq_hr.T*ref.heating
        eq_hr.Q = eq_hr.heating
        # 'transport' didn't always used to exist, so only read it if possible
        if os.path.exists(dirname + '/transport'):
            trans = TransportCoeffs(dirname + '/transport')
            eq_hr.nu = trans.nu
            eq_hr.dlnu = trans.dlnu
            eq_hr.kappa = trans.kappa
            eq_hr.dlnkappa = trans.dlnkappa
            try:
                eq_hr.eta = trans.eta
                eq_hr.dlneta = dlneta # this will fail for hydro cases
                # "trans" will not have attributes eta, dlneta
            except: # if it failed, just keep the arrays zero             
                pass # (magnetism = False)
        else:
            print ("get_eq(): neither 'equation_coefficients' nor 'transport' found")
            print ("nu, dlnu, etc. will be zero")
        eq_hr.lum = get_parameter(dirname, 'luminosity')
    return eq_hr
