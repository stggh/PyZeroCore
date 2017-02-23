#######################################################################
#
# O-Scharf.py
#
# Fine-structure intensities for the O-atom
#
# Stephen.Gibson@anu.edu.au
#
########################################################################

import numpy as np
from wigner3 import Wigner3j, Wigner9j

about = """
________________________________________________________________________________

 Scharf and Godefroid  arXiv:0808.3529v1

 Eq (55) 

   \sigma (Ja,Ji) = [Ja,Si,Li] X
                    \Sum_lc [lc] (lc 1 li)^2 A [ja,Si,Li]  X
                                 ( 0 0  0)  
                    N X \Sum_jt [jt] {Sa La Ja}^2
                                     { s li ji}
                                     {Si Li Ji}


 1D2:         La=2, Sa=0, Ja=2
 3PJ:         La=1, Sa=1, Ja=2,1,0
 2P1/2,3/2:   Li=1, Si=1/2, Ji = 3/2, 1/2
________________________________________________________________________________
"""

print(about)

def sigma(Ja, La, Sa, Ji, Li, Si):
    s = 1/2
    li = 1
    wgt = (2*Ja + 1)*(2*Si + 1)*(2*Li + 1)

    sum = 0
    for lc in (0, 1, 2):
        for jt in np.arange(abs(li - 1/2), li + 0.9, 1.0):
           sum += (2*jt + 1)*Wigner9j(Sa, La, Ja, s, li, jt, Si, Li, Ji)**2
        sum *= (2*lc + 1)*Wigner3j(lc, 1, li, 0, 0, 0)**2

    return sum*wgt   #*(2*Ji+1) 

P232 = sigma(Ja=2, La=1, Sa=1, Ji=3/2, Li=1, Si=1/2)
P132 = sigma(Ja=1, La=1, Sa=1, Ji=3/2, Li=1, Si=1/2)
P032 = sigma(Ja=0, La=1, Sa=1, Ji=3/2, Li=1, Si=1/2)

P212 = sigma(Ja=2, La=1, Sa=1, Ji=1/2, Li=1, Si=1/2)
P112 = sigma(Ja=1, La=1, Sa=1, Ji=1/2, Li=1, Si=1/2)
P012 = sigma(Ja=0, La=1, Sa=1, Ji=1/2, Li=1, Si=1/2)

D232 = sigma(Ja=2, La=2, Sa=0, Ji=3/2, Li=1, Si=1/2)
D212 = sigma(Ja=2, La=2, Sa=0, Ji=1/2, Li=1, Si=1/2)

print (" Transition    Intens      ratio/P232")
print ("3P2 <- 2P3/2   {:6.4f}      {:6.3f}".format(P232, P232/P232))
print ("3P1 <- 2P3/2   {:6.4f}      {:6.3f}".format(P132, P132/P232))
print ("3P0 <- 2P3/2   {:6.4f}      {:6.3f}".format(P032, P032/P232))
print ("   sum         {:6.4f}\n".format(P232 + P132 + P032))

print ("3P2 <- 2P1/2   {:6.4f}      {:6.3f}".format(P212, P212/P232))
print ("3P1 <- 2P1/2   {:6.4f}      {:6.3f}".format(P112, P112/P232))
print ("3P0 <- 2P1/2   {:6.4f}      {:6.3f}".format(P012, P012/P232))
print ("   sum         {:6.4f}\n".format(P212 + P112 + P012))

print ("1D2 <- 2P3/2   {:6.4f}      {:6.3f}".format(D232, D232/P232))
print ("1D2 <- 2P1/2   {:6.4f}      {:6.3f}".format(D212, D212/P232))
print ("   sum         {:6.4f}\n".format(D232 + D212))

print ("ratio of 3P2s (3/2)/(1/2) {:5.2f}\n".format(P232/P212))

print ("ratio of 1D2s (3/2)/(1/2) {:5.2f}".format(D232/D212))

# note does not include the Bolztmann factor
# I_1/2 / I_3/2 = 1/2 exp(-177.1 h c 100/kT)
