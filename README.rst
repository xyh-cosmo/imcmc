=====
imcmc
=====
:imcmc: A C++ implemented MCMC library
:Version: 0.0.4
:Author: Youhua Xu
:Homepage: http://lss.bao.ac.cn/~xyh/

Description
============

This code aims to provide easy to use interfece to MCMC simulations.  Currently only the Affine-Invariant Ensemble Sampler algorithm has been implemented, in futher I will add other commonly used algorithms / samplers into this  
library, for example the classical Metropolis-Hastings and Hamiltonian Monte Carlo. Besides I would also like to add an interface to MultiNest.

Updates
=========
**** 18-Oct-2015 ****
1) Updated the output format, the same as CosmoMC . 
2) keyword 'chain_size' is replaced by sample_step.

Dependence
============
1) OpenMPI
2) GSL (using only GSL's random number generator)
