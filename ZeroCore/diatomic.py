# -*- coding: utf-8 -*-
import numpy as np
import scipy.constants as const

from scipy.integrate.quadpack import dblquad 
from scipy.integrate.quadrature import *
from scipy.special import j0, j1, jn

############################################################################
#
# diatomic.py - diatomicatomic photodetachment calculation using the zero-core
# contribution method of Stehman and Woo.
# 
# Stephen.Gibson@anu.edu.au
# February 2017
#
############################################################################

class ZeroCore:
    """ Evaluation of beta Eq.(36) and xs Eq.(34)
        from Stehman and Woo PRA 23, 2866 (1981).

    """

    def __init__(self, energy, vib, EA, Evp, Re, R0, Roo, FCF, verbose=True):

        self.EA = EA
        self.Evp = Evp
        self.Re = Re
        self.R0 = R0
        self.Roo = Roo
        self.FCF = FCF

        if verbose:
            print(" v    eKE       beta      xs")

        BETA = []
        XS = []
        for v in vib:
            BETA.append([])
            XS.append([])
            for eKE in energy:
                self.xsbeta(eKE, v)
                if verbose:
                    print("{:2d}  {:5.2f}    {:8.3f}  {:8.3f}".
                          format(v, eKE, self.beta, self.xs))
                BETA[-1].append(self.beta) 
                XS[-1].append(self.xs) 
            if verbose:
                print()

        self.beta = BETA
        self.xs = XS

    def Gamma(self, Evv):  
        return np.sqrt(2*const.m_e*Evv*const.e)/const.hbar

    def psi_d(self, r, theta, gamma):
        gr = gamma*r;
        sc = np.sin(theta)*np.cos(theta)
        polyn = (1+(1+1/gr)*3/gr)
        expn = np.exp(-gr)/r
        psi = polyn*expn*sc
        return psi 

    def psi_d2(self, r, theta, gamma):
        return  (self.psi_d(r, theta, gamma)*r)**2 * np.sin(theta)

    def I0I(self, r, theta, k, beta, gamma):
        st = np.sin(theta)
        ct = np.cos(theta)
        sb = np.sin(beta)
        cb = np.cos(beta)
        J0 = j0(k*r*sb*st)
        return st*st*np.sin(k*r*cb*ct)*J0*self.psi_d(r, theta, gamma)*r*r*r

    def I1I(self, r, theta, k, beta, gamma):
        st = np.sin(theta)
        ct = np.cos(theta)
        sb = np.sin(beta)
        cb = np.cos(beta)
        J1 = j1(k*r*sb*st)
        return st*ct*np.cos(k*r*cb*ct)*J1*self.psi_d(r,theta,gamma)*r*r*r

    def I2I(self, r, theta, k, beta, gamma):
        st = np.sin(theta)
        ct = np.cos(theta)
        sb = np.sin(beta)
        cb = np.cos(beta)
        J2 = jn(2, k*r*sb*st)
        return st*st*np.sin(k*r*cb*ct)*J2*self.psi_d(r,theta,gamma)*r*r*r

    def Rinit(self, theta):
        return np.sqrt((self.Re/2+self.R0*abs(np.cos(theta)))**2 +\
               (self.R0*np.sin(theta))**2)*const.angstrom

    def IOf(self, beta, k, v, func): 
        # double integral of I_n function
        Evv = (self.EA + self.Evp[v]-self.Evp[0])/8065.541
        Gamv = self.Gamma(Evv)
        # integration  theta=0-pi r=0.01-Roo Angtsroms
        Iv, Iv_err = dblquad(func, 0, np.pi, self.Rinit,
                             lambda x: self.Roo*const.angstrom,
                             args=(k, beta, Gamv))
        # wavefunction normalization outsize core region
        N, N_err = dblquad(self.psi_d2, 0, np.pi, self.Rinit,
                   lambda x: self.Roo*const.angstrom, args=(Gamv,))
        return Iv/N

    def xsbeta(self, eKE, v):
        CONST = (1/(4*const.pi*const.epsilon_0))*(2*const.pi/3)*\
                (const.e*const.e/const.hbar/const.c)*\
                (const.m_e/const.hbar)*0.25*1.0e20
        fudgefactor = 30.0e-9

        # Euler range for integration  0-pi
        A = np.arange(0, np.pi, 0.1)
        CI = np.zeros_like(A) 
        SI = np.zeros_like(A)

        if eKE < 1.0e-5:
            eKE = 0.01 # special case, no energy=0
        k = self.Gamma(eKE)
        lam = 1.0e7/(eKE*8065.541 + self.EA+self.Evp[v]-self.Evp[0]) # nm
        omega = 2*const.pi*const.c/(lam*1.0e-9)

        for i, angle in enumerate(A):
           I0 = self.IOf(angle, k, v, self.I0I)
           I1 = self.IOf(angle, k, v, self.I1I)
           I2 = self.IOf(angle, k, v, self.I2I)
           ci = (((I0-I2)*np.sin(angle)/2 + I1*np.cos(angle))**2)*\
                np.sin(angle)
           si = ((I0+I2)**2/4 +\
                ((I0-I2)*np.cos(angle)/2-I1*np.sin(angle))**2)*\
                np.sin(angle)
           CI[i] = ci
           SI[i] = si

        C = simps(CI, A);  # Eq. (28)
        S = simps(SI, A)   # Eq. (29)

        self.beta = (2*C - S)/(C + S)
        self.xs = CONST*self.FCF[v]*k*omega*(C + S)*fudgefactor


class O2(ZeroCore):

    def __init__(self, energy, vib=0, EA=3613, Re=1.341, R0=0.96, Roo=50.0,
                 verbose=True):

        # O2 vibrational energy levels in cm-1, relative to X-state minimum
        Evp = [787.39802, 2343.75905, 3876.57119, 5386.16439,
               6872.50036, 8335.76750]
        # Franck-Condon factors v'=0-6, calculated with CSE program
        FCF = [6.5223e-02, 2.0532e-01, 2.8625e-01, 2.3695e-01,
               1.3223e-01, 5.3448e-02, 1.6210e-02]
      
        super().__init__(energy, vib, EA, Evp, Re, R0, Roo, FCF, verbose)


class S2(ZeroCore):

    def __init__(self, energy, vib=0, EA=16752.966, Re=1.889, R0=1.51,
                 Roo=50.0, verbose=True):

        # S2 vibrational energy levels in cm-1, relative to X-state minimum
        Evp = [362.24185, 1082.19506, 1796.75653, 2505.85869,
               3209.46163, 3907.07971]
        # Franck-Condon factors v'=0-6, calculated with CSE program
        # tmp these are O2's
        FCF = [6.5223e-02, 2.0532e-01, 2.8625e-01, 2.3695e-01,
               1.3223e-01, 5.3448e-02, 1.6210e-02]

        super().__init__(energy, vib, EA, Evp, Re, R0, Roo, FCF, verbose)
