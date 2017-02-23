# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
import scipy.constants as C
from scipy.interpolate import interp1d

from zerocore import atomic

def pwr(x):
    return int(np.floor(np.log10(x))) if x > 0 else 1

def fmt(num, err):
    ''' format xxx+-yyy  number as xxx(yyy) '''
    ints = pwr(num)
    decpl = pwr(err)
    adecpl = abs(decpl)
    f = '{{}}'+"{{}}.{{}}".format(ints+adecpl, adecpl) + "f({{}})"
    return f.format(num, round(err/(10**decpl)))

#---- main ----
# experimental data
eb, xsb = np.loadtxt("data/Branscombxs.dat", unpack=True)
EbP, betaP, ebP = np.loadtxt("data/O-cBeta.dat", unpack=True)
EbP /= 8065.541
ECZ, CZ = np.loadtxt("data/O-CooperZare.dat", unpack=True)
Ebr, betabr = np.loadtxt("data/O-Breyer.dat", unpack=True)

EA = 11784.676   # O-atom electron affinity cm-1
En = np.arange (0.01,5,0.1)

r0P = 1.0e-10 # 0.96e-10

O3P = atomic.O(En, EA, r0P, -0.93) # -0.95)

O3P.xs *= 0.9e65 #5.6e65

P = interp1d(En+EA/8065.541, O3P.xs, bounds_error=False, fill_value=0)

valp = P(3.35)   # Fix me! - this plays with relative 3P/1D scaling
scalefactorp = float(6.2/valp)

O3P.xs *= scalefactorp

Emax = 6.0 
E = np.arange(EA/8065.541, Emax, 0.1)
Pxs = P(E)*scalefactorp

plt.subplot(121)
plt.plot(eb, xsb, 'o', label="Branscomb")
plt.title("cross section")
plt.subplot(122)
plt.errorbar(EbP, betaP, yerr=ebP, fmt='og', label=r"expt $^3\!P$")
plt.plot(ECZ, CZ, 'k', label="CooperZare")

plt.subplot(121)
plt.plot(E, Pxs, 'c-', label="3P")
subr = E > 3.26
plt.legend(numpoints=1, labelspacing=0.1, loc=4, frameon=False,
           prop=FontProperties(size='small'))
plt.xlabel("Photon Energy (eV)")
plt.ylabel(r"$\sigma$ (10$^{18}$ cm$^2$)")
plt.axis(xmax=4, ymax=12)

plt.subplot(122)
plt.axis(xmin=-0.1, xmax=1.5)
plt.plot(En[En<2.1], O3P.beta[En<2.1], '-g', label=r"$^3\!P$ $r_0=1.8$")
plt.plot(Ebr, betabr, 'ro', label=r"Breyer")
plt.xlabel("eKE (eV)")
plt.ylabel(r"$\beta$")
plt.title("anisotropy")
plt.legend(numpoints=1, labelspacing=0.1, loc=0, frameon=False,
           prop=FontProperties(size='small'))
plt.suptitle("O$^{-}$ photodetachment", fontsize=15)

plt.savefig("output/example_O3P.png")
plt.show()
