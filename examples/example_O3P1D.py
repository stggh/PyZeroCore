# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

import scipy.constants as C
from scipy.interpolate import interp1d
from ZeroCore import atomic

def pwr(x):
    return int(np.floor(np.log10(x))) if x > 0 else 1

def fmt(num, err):
    ''' format xxx+-yyy  number as xxx(yyy) '''
    ints = pwr(num)
    decpl = pwr(err)
    adecpl = abs(decpl)
    f = '{{}}'+"{{}}.{{}}".format(ints+adecpl, adecpl) + "f({{}})"
    return f.format(num, round(err/(10**decpl)))

def weighted_ave(x, ex):
    ssq = ex**2
    ssig = (1/ssq).sum()
    wave = (x/ssq).sum()/ssig
    return wave, 1/np.sqrt(ssig)

#---- main ----
EA = 11784.676   # O-atom electron affinity cm-1
EAD = EA + 15867.862   # O(1D2) threshold energy cm-1

# experimental data
ex, xsb = np.loadtxt("data/Branscombxs.dat", unpack=True)
exH, xsH, exsH = np.loadtxt("data/Hlavenka_3P.dat", unpack=True)
exL, xsL, exsL = np.loadtxt("data/Lee_3P.dat", unpack=True)
exW, xsW = np.loadtxt("data/Wuxs.dat", unpack=True)
exZ, xsZ = np.loadtxt("data/Zatxs.dat", unpack=True)
exZ += 1.4611

# anisotropy
EbA, betaA, ebA = np.loadtxt("data/ANU-OP2P32.dat", unpack=True)
EbD, betaD, ebD = np.loadtxt("data/beta_1D.dat", unpack=True)
EbP, betaP, ebP = np.loadtxt("data/beta_3P.dat", unpack=True)
lamH, betaH, ebH = np.loadtxt("data/O-Hanstorp-expt.dat", unpack=True)
EbH = (1.0e7/lamH - EA)/8065.541

EbD = EbD.mean()
betaD, ebD = weighted_ave(betaD, ebD)
print("beta_1D(eKE={:6.2f}) = {:5.3f}+-{:4.3f}".format(EbD, betaD, ebD))

EbP = EbP.mean()
betaP, ebP = weighted_ave(betaP, ebP)
print("beta_3P(eKE={:6.2f}) = {:5.3f}+-{:4.3f}".format(EbP, betaP, ebP))

EbD /= 8065.541
EbP /= 8065.541
ECZ, CZ = np.loadtxt("data/O-CooperZare.dat", unpack=True)

En = np.arange (0.01,5,0.1)

# 0.8e-10 for both looks good for xs, totally misses beta
r0P = 1.5e-10 #0.96e-10
r0D = 1.8e-10

O3P = atomic.O(En, EA, r0P, -0.95)
O1D = atomic.O(En, EAD, r0D, -0.95)

O3P.xs *= 0.9e65 #5.6e65
O1D.xs *= 0.35e64

P = interp1d(En+EA/8065.541, O3P.xs, bounds_error=False, fill_value=0)
D = interp1d(En+EAD/8065.541, O1D.xs, bounds_error=False, fill_value=0)

valp = P(3.35)   # Fix me! - this plays with relative 3P/1D scaling
vald = D(3.9)
scalefactorp = float(6.2/valp)
scalefactord = float((10-valp)/vald)

O3P.xs *= scalefactorp
O1D.xs *= scalefactord

Emax = 6.0 #float(int((En+EAD/8065.541)[-1]))
E = np.arange(EA/8065.541, Emax, 0.1)
Pxs = P(E)*scalefactorp
Dxs = D(E)*scalefactorp

#------- plots ---------------
gs = gridspec.GridSpec(1, 2, left=0.1, bottom=0.12, right=0.95, top=0.85,
                             wspace=0.35, hspace=0.2)
ax0 = plt.subplot(gs[0, 0])
ax1 = plt.subplot(gs[0, 1])

lb, = ax0.plot(ex, xsb, 'C7o', label="Branscomb")
lh = ax0.errorbar(exH, xsH, yerr=exsH, color='C8', marker='o', ls='',
                  label="Hlavenka")
ll = ax0.errorbar(exL, xsL, yerr=exsL, color='C9', marker='o', ls='',
                  label="Lee")
ax0.plot((3.647, 3.647), (7, 10.7), 'C3--', lw=0.4)
ax0.annotate("340nm", (3.2, 11), color='C3', fontsize='small')
ax0.plot((4.66, 4.66),(7, 10.7), 'C3--', lw=0.4)
ax0.annotate("266nm", (4.3, 11), color='C3', fontsize='small')

subr = E <= 3.26
lZC3P, = ax0.plot(E[subr], Pxs[subr], 'C2-', label=r"$^{3}P$ zero core")
subr = E > 3.26
lZC1D, = ax0.plot(E[subr], Pxs[subr]+Dxs[subr], 'C3-', label=r"$^{1}D$ zero core")
lW, = ax0.plot(exW, xsW*0.71, 'C6--', label=r"Wu $\times 0.71$")
lZ, = ax0.plot(exZ, xsZ*0.65, 'C4--', label=r"Zatsarinny $\times 0.65$")

# lh
ax0.legend(handles=[lb, lh, ll, lW, lZ, lZC3P, lZC1D], labelspacing=0.1,
                    loc=4, frameon=False, fontsize='small')
ax0.set_title("cross section")
ax0.set_xlabel("Photon Energy (eV)")
ax0.set_ylabel(r"$\sigma$ (10$^{18}$ cm$^2$)")
ax0.axis(xmax=5, ymax=12)

# anisotropy plot
ax1.errorbar(EbA, betaA, yerr=ebA, fmt='C2o', mfc='w',
             label=r"expt $^3\!P_{2}$")
ax1.errorbar(EbH, betaH, yerr=ebH, fmt='C2s', mfc='w',
             label=r"Hanstorp $^3\!P$")
ax1.errorbar([EbP], [betaP], yerr=[ebP], fmt='C2o', label=r"expt $^3\!P$")
ax1.errorbar([EbD], [betaD], yerr=[ebD], fmt='C3o', label=r"expt $^1\!D$")
ax1.errorbar([3.2], [0.0], yerr=[0.1], fmt='C2+', label=r"$^3\!P$ Domesle")
ax1.errorbar([1.23], [-0.9], yerr=[0.1], fmt='C3+', label=r"$^1\!D$ Domesle")
ax1.plot(ECZ, CZ, 'k', label="CooperZare")

ax1.axis(xmin=-0.1, xmax=3.5)
ax1.plot(En[En<2.1], O3P.beta[En<2.1], 'C2-', label=r"$^{3}P$ $r_{0}=0.96\, \AA$")
ax1.plot(En[En<2.3], O1D.beta[En<2.3], 'C3-', label=r"$^{1}D$ $r_{0}=1.8\, \AA$")
ax1.set_xlabel("eKE (eV)")
ax1.set_ylabel(r"anisotropy parameter $\beta$")
ax1.set_title("anisotropy")
ax1.legend(numpoints=1, labelspacing=0.1, loc=0, frameon=False, fontsize='small')

plt.suptitle("O$^{-}$ photodetachment", fontsize=15)

plt.savefig("output/example_O3P1D.png", dpi=75)
plt.show()
