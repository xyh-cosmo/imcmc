===================
imcmc
===================
:imcmc: A C++ implemented MCMC library
:Version: 0.0.3
:Author: Youhua Xu
:Homepage: http://lss.bao.ac.cn/~xyh/

Description
============

This code aims to provide easy to use interfece to MCMC simulations.  Currently only the Affine-Invariant Ensemble Sampler algorithm has been implemented, in futher I will add other commonly used algorithms / samplers into this  
library, for example the classical Metropolis-Hastings and Hamiltonian Monte Carlo. Besides I would also like to add an interface to MultiNest.

Dependence
============
1) OpenMPI
2) GSL (use GSL's random number generator)
