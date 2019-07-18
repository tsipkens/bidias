# UBC-aerosol-inversion

A program to invert CPMA-DMA data to find the two-dimensional
mass-mobility distribution.

### Description

This program is organized into several packages and classes. The `main.m`
scripts in the top directory of the code can be called to demonstrate
the code.

#### Classes

###### @Grid

Grid is a class developed to discretize mass-mobility space. This is
done using a simple rectangular grid that can have linear, logarithmic
or custom spaced elements along the edges.

###### @Phantom

Phantom is a class developed to contain the parameters and other information
for the phantom distributions that are used in testing the different inversion
methods.

#### Packages

###### +invert

Contains various functions used to invert the measured data for the desired
two-dimensional distribution. This includes implementations of least-squares,
Tikhonov regularization, Twomey, Twomey-Markowski (including using the method
of Buckley et al. (2017)), and the multiplicative algebraic reconstruction
technique (MART). Also included are functions that, given the true distribution,
can determine the optimal number of iterations or the optimal regularization
parameter.

###### +kernel

Evaluates the transfer function of the DMA and particle mass analyzer (such
as the CPMA or APM). The primary function within the larger program is to
generate a matrix `A` that acts as the forward model. This package references
an imported submodule, `UBC-tfer-PMA`, which contains the package `+tfer_PMA.m`
that is used in evaluating the transfer function of the particle mass
analyzers and some standard reference functions used in `tfer_DMA.m`

###### +tools

A series of utility functions that serve various purposes, including printing
a text-based progress bar (based on code from Samuel Grauer) and a function
to convert mass-mobility distributions to effective density-mobility
distributions.

----------------------------------------------------------------------

#### License

This software is licensed under an MIT license (see the corresponding file
for details).


#### Contact information and acknowledgements

This program was largely written and compiled by Timothy Sipkens
([tsipkens@mail.ubc.ca](mailto:tsipkens@mail.ubc.ca)) while at the
University of British Columbia.

This distribution includes code snippets from the code provided with
the work of Buckley et al.
([https://doi.org/10.1016/j.jaerosci.2017.09.012](https://doi.org/10.1016/j.jaerosci.2017.09.012)),
who used a Twomey-type approach to derive two-dimensional mass-mobility
distributions.

Also included is a reference to code designed to quickly evaluate
the transfer function of particle mass analyzers (e.g. APM, CPMA) by
Sipkens et al. (Under review).

Information on the provided colormaps can be found in an associated
README in the `cmap` folder.
