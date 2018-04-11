PyZeroCoreContribution README
=============================

.. image:: https://travis-ci.org/stggh/PyZeroCore.svg?branch=master
    :target: https://travis-ci.org/stggh/PyZeroCore


Introduction
------------
``PyZeroCoreContribution`` is a python implementation of the Zero-Core-contribution method for the calculation photodetachment cross sections and anisotropy parameters of negative-ions.  Details are in the following publications [1, 2, 3, 4] and others.


Method
------
The photodetachment cross section of a negative ion may be calculated from:

.. math::

   \frac{d\sigma}{d\Omega} = \frac{e^2}{\hbar c}\frac{\omega k^2}{2\pi}|M_{k0}|^2,

where :math:`M_{k0} = < k | \sum_i \hat{e} \cdot \vec{r}_i | 0 >`



Implementation
--------------
This is a Python implementation of the algorithms given in the above publications. *(more detail to follow)*


Reference
---------

1. `R.M. Stehman and S.B. Woo "Zero-core-contribution model and its application to photodetachment of atomic negative-ions" Physical Review A 20, 281-290 (1979) <http://dx.doi.org/10.1103/PhysRevA.20.281>`_

2. `R.M. Stehman and S.B. Woo "Zero-core-contribution calculation of photodetachment cross-sections of O2- and S2-" Physical Review A 28, 2866-2876 (1981) <http://dx.doi.org/10.1103/PhysRevA.23.2866>`_

3. `W.B. Clodius, R.M. Stehman, and S.B. Woo "Zero-core-contribution calculation of a polyatomic photodetachment cross-section - NO2-" Physical Review A 2*, 760-765 (1983) <http://dx.doi.org/10.1103/PhysRevA.28.760>`_

4. `W.B. Clodius, R.M. Stehman, and S.B. Woo, "Zero-core-contribution calculation of photodetachment characteristics of heteronuclear diatomic anions" Physical Review A 28, 1160-1163 (1983) <http://dx.doi.org/10.1103/PhysRevA.28.1160>`_


Citation
--------
If you find PyZeroCore useful in your work please consider citing this project.

.. image:: https://zenodo.org/badge/82896133.svg
   :target: https://zenodo.org/badge/latestdoi/82896133
