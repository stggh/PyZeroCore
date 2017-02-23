# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as C
from zerocore import polyatomic

#---- main ----
Em=0.0; Ex=8.0; dE=0.05
try:
    data = input("eKE range min, max, step (eV) [{:g}, {:g}, {:g}]? ".format(Em,Ex,dE))
    Em, Ex, dE = data
except:
    pass

E = np.arange(Em, Ex+dE/2, dE)      # eKE range  (start,end,step)

Vm=0; Vx=5; dV=1
try:
    data = input ("v' range min, max, step (eV) [%d,%d,%d]? "%(Vm,Vx,dV))
    Vm, Vx, dV = data
except:
    pass
Vx += 1
if Vx > 6: Vx = 6
V=range(Vm,Vx,dV)  # vibrational QN


no2 = polyatomic.NO2(E)

plt.plot(no2.xsbeta)
plt.show()

