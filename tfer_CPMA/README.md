# Transfer Function README

The attached functions and script are intended to reproduce the results of
this paper.

## Functions to evaluate transfer functions

This program includes various functions evaluating the transfer
function for the various cases presented in the associated work. They share a
common inputs:

1. m_star - the setpoint mass,

2. m - the masses at which the transfer function will be evaluated,

3. d - the mobility diameter (either as a scalar or as a vector with the same
  length as the masses at which the transfer function is to be evaluated),

4. z - the integer charge state (either as a scalar or as a vector with the same
  length as the masses at which the transfer function is to be evaluated),

5. prop - a struct that contains the properties of the particle mass analyzer
  (a sample script to generate this quantity is include as prop_CPMA.m), and

6. varargin (optional) - name-value pairs to specify either the equivalent
  resolution, inner electrode angular speed, or voltage.

## Main script

This script is included to demonstrate evaluation of the transfer function over
multiple cases.
