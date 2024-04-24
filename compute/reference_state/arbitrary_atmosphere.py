import numpy as np
from scipy.integrate import simps

def arbitrary_atmosphere(r, s, dsdr, d2sdr2, g, dgdr, r0, T0, p0, cp, gamma):
    cv = cp/gamma
    R = cp*(gamma - 1.0)/gamma
    integrand = g*np.exp(-s/cp)
    ir0 = np.argmin(np.abs(r - r0))
    s0 = s[ir0]
    
    integral = np.zeros_like(r)
    for ir in range(len(r)):
        if ir <= ir0:
            integral[ir] = -simps(integrand[ir:ir0 + 1], r[ir:ir0 + 1])
        else:
            integral[ir] = simps(integrand[ir0:ir + 1], r[ir0:ir + 1])
            
    T = -(1.0/cp)*np.exp(s/cp)*integral + T0*np.exp((s - s0)/cp)   
    p = p0*np.exp(-(s - s0)/R)*(T/T0)**(gamma/(gamma - 1.0))
    rho0 = p0/(R*T0)
    rho = rho0*np.exp(-(s - s0)/R)*(T/T0)**(1.0/(gamma - 1.0))
    
    dlnT = dsdr/cp - g/(cp*T)
    dlnrho = (1.0/(gamma - 1.0))*(dlnT - dsdr/cv)
    dlnp = dlnT + dlnrho
    
    d2lnT = d2sdr2/cp - dgdr/(cp*T) + g/(cp*T) * dlnT
    d2lnrho = (1.0/(gamma - 1.0))*(d2lnT - d2sdr2/cv)
    
    return T, rho, p, dlnT, dlnrho, dlnp, d2lnrho

def arbitrary_atmosphere_nd(rr, dsdr, g, rr0, rho0, tmp0, Di, gamma):
    # compute entropy(r) and higher derivatives
    entr = indefinite_integral(dsdr, rr, rr0)
    d2sdr2 = np.gradient(dsdr, rr) 
    dgdr = np.gradient(g, rr)

    # compute rho(r) and T(r)
    tmp = np.exp(entr)*(tmp0 - Di*indefinite_integral(g*np.exp(-entr), rr, rr0))
    rho = rho0*np.exp(-gamma*entr/(gamma-1.))*(tmp/tmp0)**(1./(gamma-1.))

    # compute higher derivatives of rho and T
    dtdr = dsdr*tmp - Di*g
    dlnt = dtdr/tmp
    d2tdr2 = d2sdr2*tmp + dsdr*dtdr - Di*dgdr
    d2lnt = d2tdr2/tmp - dlnt**2

    dlnrho = -gamma*dsdr/(gamma-1.) + dlnt/(gamma-1.)
    d2lnrho = -gamma*d2sdr2/(gamma-1.) + d2lnt/(gamma-1.)
    
    return rho, tmp, dlnrho, d2lnrho, dlnt
