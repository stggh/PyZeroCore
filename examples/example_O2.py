import numpy as np
from zerocore import diatomic
from scipy.interpolate import splrep, splev

import matplotlib.pyplot as plt
import time

# energy range for the calculation
Em = 0
Ex = 2
dE = 0.1
E = np.arange(Em, Ex+dE/2, dE) 

# vibrational levels
Vm=0
Vx=5
dV=1
Vx += 1
if Vx > 6: 
    Vx = 6
V = np.arange(Vm, Vx, dV)  # vibrational QN

print("Evaluation of beta Eq.(36) and xs Eq.(34)")
print("from Stehman and Woo PRA 23, 2866 (1981).\n")

t0 = time.time()
O2zc = diatomic.O2(E, V)
t1 = time.time()
print("execution time {:.2f} seconds".format(t1-t0))

beta = O2zc.beta
xs = O2zc.xs

col = ['r', 'b', 'g', 'y', 'c', 'm', 'brown'] 
plt.subplot(121)
plt.suptitle(r"Stehman and Woo 1981 R0={:g} $\AA$".format(O2zc.R0))
plt.xlabel(r"eKE(eV)")
plt.ylabel(r"$\beta$")
for v in V:
   plt.plot(E, beta[v], '-', color=col[v], label=str(v))

for v in np.arange(7):
    Xbeta = np.loadtxt("data/O2-/X{:1d}.dat".format(v), unpack=True)
    plt.errorbar(Xbeta[0], Xbeta[1], yerr=Xbeta[2], color=col[v], linestyle='', marker='o')

plt.legend(loc=0, frameon=False, labelspacing=0.1, numpoints=1)

plt.subplot(122)
plt.xlabel(r"$h\nu$(eV)")
plt.ylabel(r"$\sigma$")
Eall = np.arange(Em, 3.5, dE)
xstotal = np.zeros_like(Eall)
for v in V:
    hnu = E + (O2zc.EA +(O2zc.Evp[v]-O2zc.Evp[0]))/8065.541
    spl = splrep(hnu, xs[v], s=0)
    xstotal += splev(Eall, spl, der=0, ext= 3) 
    plt.plot(hnu, xs[v], color=col[v], label=str(v))

plt.plot(Eall, xstotal, 'k-', label="total")

Cosby = np.loadtxt("data/O2-/xs/Cosby.dat", unpack=True)
plt.plot(*Cosby, 'bo', label="Cosby")

plt.legend(loc=0, frameon=False, labelspacing=0.1, numpoints=1)
plt.savefig("output/example_O2.png", dpi=100)
plt.show()
