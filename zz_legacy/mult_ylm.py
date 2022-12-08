import sys,os
sys.path.append(os.environ['rapp'])
from sympy.physics.quantum.cg import CG
import numpy as np
from spectral_utils import SHT

def my_CG(l1,m1, l2,m2, l,m, sht):
    ''' kind of 'clebsch gordon coefficient'
    to figure how two spherical harmonics contribute to Ylm'''
    tmp = np.zeros((sht.n_l, sht.nm), dtype='complex')
    tmp2 = np.zeros((sht.n_l, sht.nm), dtype='complex')
    tmp[l1,m1] = 1.+0.j
    tmp2[l2,m2] = 1.+0.j
    return (sht.to_spectral(sht.to_physical(tmp)*sht.to_physical(tmp2))[l,m])

def my_CG_sympy(l1,m1, l2,m2, l,m):
    return CG(l1,m1, l2,m2, l,m).doit().evalf()
  
def prodlm(A, B, ntheta):
    sht = SHT(ntheta)
    return sht.to_spectral(sht.to_physical(A)*sht.to_physical(B))
