import numpy as np
import matplotlib.pyplot as plt
import ZeroCore.wigner3 as w3j
from cmath import *
from scipy.interpolate import interp1d
from scipy.integrate import tplquad
from scipy.special import sph_harm, spherical_jn 
import scipy.constants as C

class ZeroCore:
      """
      polyatomic photodetachment
      based on Stehman and Woo Phys. Rev. A. 20, 281-290 (1979) and
                               Phys. Rev. A. 28, 760-765 (1983)
      
      """

    def __init__(self, eKE, EA, En, G, R, Rn, Ro, alpha, vib=0):

        const = (1/(4*C.pi*C.epsilon_0))*(2*C.pi/3)*(C.e*C.e/C.hbar/C.c)*\
                (C.m_e/C.hbar)*0.25*1.0e20

        self.EA = EA
        self.En = En
        self.G = G
        self.R = R
        self.Rn = Rn
        self.Ro = Ro
        self.Roo = 50
        self.alpha = alpha
        
        self.xsbeta(eKE, vib)


    def Rinit(self, theta):
        return np.sqrt((self.Re/2+self.R0*abs(np.cos(theta)))**2 +\
               (self.R0*np.sin(theta))**2)*C.angstrom

    def Gamma(self, E):
        return np.sqrt(2*C.m_e*E)/C.hbar

    def N(self, gam, R):
        return np.sqrt(3*gam*np.exp(2*gam*R)*(1+2/gam/R)/2/C.pi)

    def X(self, R, alpha): 
        return R*np.cos(alpha/2)

    def Y(self, R, alpha): 
        return R*np.sin(alpha/2)

    def phi(self, gam, R, r, x):
        return  self.N(self.gam, self.R)*(1 + 1/self.gam/r)*\
                np.exp(-self.gam*r)*x/r**2

    def phi_0x(self, gam, r2, r3, x, X):
        return self.N(gam,  self.Ro)*(self.phi(gam,  self.Ro, r2, x+X) +\
               self.phi(gam,  self.Ro, r3, x+X))

    def phi_0y(self, gam, r2, r3, y, Y):
        return self.N(gam, self.Ro)*(self.phi(gam,  self.Ro, r2, y+Y) -\
               self.phi(gam,  self.Ro, r3, y-Y))

    def phi_0(self, alpha, delta, gam, R, r2, r3, x, y):
        phi0x = self.phi_0x(gam, r2, r3, x, self.X(R, alpha))*np.cos(delta)
        phi0y = self.phi_0y(gam, r2, r3, y, self.Y(R, alpha))*np.sin(delta)
        return phi0x + phi0y

    def phi_i(self, cn, co, alpha, gam, Rn, R0, x, y): 
        return cn*phi(gam, Rn, x, r) + co*self.phi(gam, Ro, x, r)
  
    def Rinit(self, theta):
        return np.sqrt((self.Re/2+self.R0*abs(np.cos(theta)))**2 +\
               (self.R0*np.sin(theta))**2)*const.angstrom

    def Mki(self, k, l, mk, m): 
        def f(r, theta, phi): 
            return spherical_jn(l, k*r)*sph_harm(l, m, theta, phi)\
                   *np.sin(theta)
        overlap, err = tplquad(f, self.Rinit,
                               lambda x: self.Roo*const.angstrom,
                               args=(k, theta))
        return overlap

    def Mii(self, k, m): 
        def f(r, theta, phi): 
            return self.phi(self.gam, *r*sph_harm(1, m, theta, phi)\
                   *phi_i(gam, R, r, x)
        overlap, err = tplquad(f, self.Rinit(), lambda x: self.Roo*C.angstrom,
                               lambda r: 0, lambda r: pi,
                               lambda r: 0, lambda r: 2*pi)
        return overlap

    def Mfi(self, k, l, mk, m):
        return self.Mki(k, l, mk, m) - self.Oki(k, l, mk)*self.Mii(k, m)

    def CplusS(self, k, r):
        sum = 0
        for m in (-1, 0, 1):
            for l in range(0, 4):
                for mk in range(-l, l+1):
                    sum += self.Mfi(k, l, mk, m)**2 
        return sum*2/3/C.pi

    def xsbeta(self, eKE, vib):
        k = self.Gamma(eKE)
        lam = 1.0e7/(eKE*8065.541 + self.EA+self.G(vib, 0)) # nm
        omega = 2*C.pi*C.c/(lam*1.0e-9)
        r = np.arange(self.Ro, 50.0, 0.1)
        self.xs = (16*C.pi**2/3)*(C.e**2/C.hbar/C.c)*(C.m_e*k*omega/C.hbar)*\
                  self.CplusS(k, r)
"""
    def CminusS:
        sum = 0
        for m in (-1, 0, 1):
            for mp in (-1, 0, 1):
                for l in range(0, 4):
                    for mk in range(-l, l+1):
                    
        return
"""

class NO2(ZeroCore):

    def __init__(self, eKE, vib=0):
        EA = 18332.974693

        R = 1.15e-10     # N-O bond length (m)
        Rn = 1.115e-10    # N atom zero core radius (m)
        Ro = 0.96e-10     # O atom zero core radius (m)
        alpha = 119.5*C.pi/180 # O-N-O bond angle (radians)

        G     = lambda v1,v2: 750*(v2+1/2) - 750/2 - 213*v1
        En    = lambda J,K,A,B:  B*J*(J+1)+(A-B)*K*K

        super().__init__(eKE, EA, En, G, R, Rn, Ro, alpha, vib)
