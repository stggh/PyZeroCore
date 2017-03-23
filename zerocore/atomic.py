# -*- coding: utf-8 -*-
import numpy as np
import scipy.constants as C

############################################################################
#
# atomic.py - atomic photodetachment calculation using the zero-core
# contribution method of Stehman and Woo.
# 
# Stephen.Gibson@anu.edu.au
# February 2017
#
############################################################################


class ZeroCore:
    """Photodetachment cross section and anisotropy parameter calculation
       based on Stehman and Woo Phys. Rev. A. 20, 281-290 (1979).

    """

    def __init__(self, en, EA, r0, cosdelta):
        """
        Parameters
        ----------
        en: numpy 1d-array of floats
            energies (in eV) at which to evaluate cross section and 
            anisotropy parameter
    
        EA: float
            electron affinity or threshold energy in cm-1 or eV
        
        r0: float
            zero core radius in metres
    
        cosdelta: float
            cosine of phase shift between outgoing partial waves
  

        """
        self.xsbeta(en, EA, r0, cosdelta)
       

    def _N2(self, gam, r0): 
        return 2*gam*np.exp(2*gam*r0)/(1+2/(gam*r0))

    def _K(self, E): 
        return np.sqrt(2*C.m_e*E)/C.hbar

    def _Rsp(self, k, gam, r0, N): 
       gk = gam*k
       gk2 = gam**2 + k**2
       kr0 = k*r0
       CT = (3*gam**2*k + k**3)/gk2 + gk*r0
       ST =  2*gam**3/gk2 + gam**2*r0
       return N*np.exp(-gam*r0)*(CT*np.cos(kr0)+ST*np.sin(kr0))/(gk*gk2)

    def _Rdp(self, k, gam, r0, N): 
       gk = gam*k
       gk2 = gam**2 + k**2
       kr0 = k*r0
       CT = (k**4 + 6*gk**2 + 3*gam**4)/(gk*k*gk2**2) + r0/gk2
       ST = 3/(gk*k**2*r0) + (3*k**2 + gam**2)/(k*gk2**2) - gam*r0/(k*gk2)
       return N*np.exp(-gam*r0)*(CT*np.cos(kr0) - ST*np.sin(kr0)) 

    def xsbeta(self, en, EA=11784.676, r0=1.5e-10, cosdelta=-0.95):
       EA *= C.e/8065.541 if EA > 10 else 1
       gam = self._K(EA)
       N = self._N2(gam, r0)

       # does not like en = 0. Set 0, and any negative energies to 
       # small value, close enough to threshold
       neg_ens = en <= 0
       en[neg_ens] = 0.01

       k  = self._K(en*C.e)
       omega = en*C.e + EA
       R_sp = self._Rsp(k, gam, r0, N)
       R_dp = self._Rdp(k, gam, r0, N)
       xs = 2*R_dp**2 + R_sp**2
       beta = 2*(R_dp**2 - 2*cosdelta*R_dp*R_sp)/xs
       xs *= (8*C.pi*C.e**2*C.m_e*k*omega/9/C.hbar**2/C.c)

       self.xs = xs
       self.beta = beta
       self.Rsp = R_sp
       self.Rdp = R_dp
       # in general include l(l-1) etc factors, see Cooper+Zare


class O(ZeroCore):
    """O- photodetachment cross section and anisotropy parameter calculation
       based on Stehman and Woo Phys. Rev. A. 20, 281-290 (1979).

    """

    def __init__(self, en, EA=11784.676, r0=1.5e-10, cosdelta=-0.95):

        super().__init__(en, EA, r0, cosdelta)
