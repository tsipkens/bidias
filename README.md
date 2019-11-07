# MATLAB tools for 2D inversion of aerosol characteristics (mat-2d-aerosol-inversion)

[![DOI](https://zenodo.org/badge/190667091.svg)](https://zenodo.org/badge/latestdoi/190667091)
[![MIT license](https://img.shields.io/badge/License-MIT-blue.svg)](https://lbesson.mit-license.org/)

A program to invert CPMA-DMA data to find the two-dimensional
mass-mobility distribution associated with [Sipkens et al. (2019)][1].

This program is organzed into several:
classes (folders starting with the @ symbol),
packages (folders starting with the + symbol), and
scripts that form the base of the program.


## Scripts in upper directory

#### Main scripts (`main*.m`)

The `main*.m` scripts in the top directory of the code can be called to
demonstrate use of the code. Scripts to execute this program should be
structured as follows:

1. Optionally, one can define a phantom used to generate synthetic data and a
ground truth. The `@Phantom` class, described below, is designed to
perform this task. The results is an instance of the `@Grid` class, which is
also described below, and a vector, `x_t`, that contains a vectorized form of
the phantom distribution, defined on the output grid.

2. One must now generate a model matrix, `A`, which relates the distribution,
`x`, to the data, `b`, such that **Ax** = **b**. This requires one to compute
the transfer functions of all of the devices involved in the measurement
as well as the grids on which `x` and `b` are to be defined.
For phantom distributions, the grid for `x` can generated using
the `@Phantom` class. In all other cases, the grid for `x` and `b` can be
generated by creating an instance of the `@Grid` class described below.

3. Next, one must define some set of data in an appropriate format. For
simulated data, this is straightforward: `b = A*x;`. For experimental data, one
must have defined a grid for the data as part of Step 2, and then the data
should be imported and vectorized to match that grid. Also in this step, one
should include some definition of the expected uncertainties in each point
in `b`, encoded in the matrix `Lb`. For those cases involving counting noise,
this can be approximated as `Lb = theta*diag(sqrt(b));`, where `theta` is
related to the total number of particle counts as described in
[Sipkens et al. (2019)][1].

4. With this information, one can proceed to implement various inversion
approaches, such as those available in the `+invert` package described below.
Preset groupings on inversion approaches are available in the
`run_inversions*.m` scripts, also described below.

5. Finally, one can post-process and visualize the results as desired. The
`@Grid` class allows for a simple visualization of the inferred distribution
by calling the `plot2d_marg` method of this class. This plots both the
retrieved distribution as well as the marginalized distribution on each of
the axes, taking the reconstruction (e.g. `x_tk1`) as an input.

Of particular note, the `main_jas19.m` script is designed to replicate the
results in the associated paper [Sipkens et al. (2019)][1].

#### Scripts to run a series of inversion methods (`run_inversions*.m`)
As noted above, these scripts are intend to bundle a series of
inversion methods into a single line of code in the `main*.m` scripts.
 This can include optimization routines, included in the `+invert`
 package, which run through several values of the regularization parameters.

The lettered scripts roughly perform as follows:

`run_inversions_a.m` - Attempts to optimize the regularization parameter in
the Tikhonov, MART, Twomey, and Twomey-Markowski approaches.

`run_inversions_b.m` - Re-runs inversion at the set of optimal parameters
produced by *run_inversions_a.m*. Can be modified to adjust the optimization approach used in the Tikhonov solvers (e.g. specifying the `'non-neg'` option).

`run_inversions_c.m` - A simple set of the Tikhonov and Twomey approaches where
the user must explicitly set the regularization parameter of the Tikhonov
schemes.

`run_inversions_d.m` - Run the inversion methods multiple times and time the
length of time required to produce a reconstruction.

## Classes

#### @Grid

Grid is a class developed to discretize mass-mobility space. This is
done using a simple rectangular grid that can have linear, logarithmic
or custom spaced elements along the edges. Methods are designed
to make it easier to deal with gridded data, allowing users to reshape
vectorized data back to a 2D grid (`reshape`) or vice versa. Other
methods allow for plotting the 2D representation of vector data (`plot2d`) or
calculate the gradient of vector data (`grad`). More information is available
in the class definition.

Both the **b** and **x** vectors are defined with respect to an instance of
this class. The vectors are arranged such that the first entry corresponds
to the smallest mass and mobility diameter. The vector proceeds, first with
increasing mass and then with increasing mobility diameter. Vectorizing the
2D gridded data can be done using the colon operand, i.e. `x(:)`, or using
the `vectorize` method.

#### @Phantom

Phantom is a class developed to contain the parameters and other information
for the phantom distributions that are used in testing the different inversion
methods. Currently, the phantom class is programmed to primarily produce
bivariate lognormal distributions and secondarily distributions
that are lognormal with mobility and conditional normal for mass
following [Buckley et al. (2017)][3]. The four sample phantoms from
[Sipkens et al. (2019)][1] can be called using strings encompassing
the distribution numbers or names from that work (e.g. the demonstration phantom
can be generated using `'1'`).

The Phantom class parameterizes the aerosol distribution in two
possible ways:

1. Most generally, the class parameterized the distribution
using a mean, `mu`, and covariance matrix, `Sigma`. For lognormal-lognormal
distributions, the mean and covariance are given in
[log<sub>10</sub>*m*, log<sub>10</sub>*d*]<sup>T</sup>
space. For phantoms of the form provided by [Buckley et al. (2017)][3]
are lognormal in mobility diameter space and conditionally normally
distributed in mass space.

2. The distribution is parameterized using the typical
mass-mobility parameters, stored in the `p` field. This includes
parameters, such as the geometric mean diameter, `dg`;
mass-mobility exponent, `Dm`; and the effective density of particles
with a mobility diameter of 100 nm `rho_100`.

Methods to convert between these parameterizations are provided
with the Phantom class as the `p2cov` and `cov2p` methods.

For experimental data, the Phantom class can also be used to derive
morphological parameters from the reconstructions. Of particular note,
the `fit` method of the Phantom class, takes a reconstruction, `x` and
the grid on which it is defined and creates a bivariate lognormal
phantom that most resembles the data. This done using least squares
analysis. The `p` properties of the Phantom class then contains many of the
morphological parameters of interest to practitioners measuring
mass-mobility distributions.

## Packages

#### +invert

The invert package contains various functions used to invert the measured data
for the desired two-dimensional distribution. This includes implementations of
least-squares, Tikhonov regularization, Twomey, Twomey-Markowski (including using
the method of [Buckley et al. (2017)][3]), and the multiplicative algebraic
reconstruction technique (MART).

Details on these approaches to inversion are provided in the
associated paper, [Sipkens et al. (2019)][1].

Development is underway on the use of an
exponential covariance function to correlate pixel values and reduce
reconstruction errors [Sipkens et al. (Under preparation)][4]..

#### +optimize

This package mirrors the content of the +inver package but,
given the true distribution, aims to determine the optimal number of
iterations for the Twomey and MART schemes or the optimal regularization
parameter for the Twomey-Markowski and Tikhonov methods.


#### +tfer_PMA

This is imported from a package distributed with [Sipkens et al. (2019)][2].
This package is used in evaluating the transfer function of the particle mass
analyzers (PMAs), such as the aerosol particle mass analyzer (APM) and centrifugal
particle mass analyzer (CPMA). The package also contains some standard reference
functions used in `tfer_DMA.m`.

The original repository can be found at
[https://github.com/tsipkens/mat-tfer-pma](https://github.com/tsipkens/mat-tfer-pma).
The current implementation corresponds to v1.3 of that code.

#### +kernel

This package is used to evaluate the transfer function of the DMA and
particle mass analyzer (such as the CPMA or APM). The primary function
within the larger program is to generate a matrix `A` that acts as the
forward model. This package references the `+tfer_PMA` package, noted
above.

#### +tools

A series of utility functions that serve various purposes, including printing
a text-based progress bar (based on code from
[Samuel Grauer](https://www.researchgate.net/profile/Samuel_Grauer))
and a function to convert mass-mobility distributions to effective
density-mobility distributions.

----------------------------------------------------------------------

#### License

This software is licensed under an MIT license (see the corresponding file
for details).

#### How to cite

This work can be cited in two ways.

1. If the methods are used, but the code is not,
please cite [Sipkens et al. (2019)][1].
Note that if the Twomey-Markowski approach is used,
one should also cite [Buckley et al. (2017)][3].

2. If this code is used directly, cite both this code
(including the DOI, included at the top) and the associated paper.

#### Contact information and acknowledgements

This program was largely written and compiled by Timothy Sipkens
([tsipkens@mail.ubc.ca](mailto:tsipkens@mail.ubc.ca)) while at the
University of British Columbia.

This distribution includes code snippets from the code provided with
the work of [Buckley et al. (2017)][3],
who used a Twomey-type approach to derive two-dimensional mass-mobility
distributions. Much of the code from that work has been significantly
modified in this distribution.

Also included is a reference to code designed to quickly evaluate
the transfer function of particle mass analyzers (e.g. APM, CPMA) by
[Sipkens et al. (2019)][2].

The authors would also like to thank Samuel Grauer
for consulting on small pieces of this code (such as
the MART code).

Information on the provided colormaps can be found in an associated
README in the `cmap` folder.

#### References

1. [Sipkens, T. A., Olfert, J. S., & Rogak, S. N. (2019). Inversion methods to determine two-dimensional aerosol mass-mobility distributions: A critical comparison of established methods. *J. Aerosol Sci.*][1]
2. [Sipkens, T. A., Olfert, J. S., & Rogak, S. N. (2019). New approaches to calculate the transfer function of particle mass analyzers. *Aerosol Sci. Technol.* doi: 10.1080/02786826.2019.1680794][2]
3. [Buckley, D. T., Kimoto, S., Lee, M. H., Fukushima, N., Hogan Jr, C. J. (2017). Technical note: A corrected two dimensional data inversion routine for tandem mobility-mass measurements. *J. Aerosol Sci.* 114, 157-168.][3]
4. [Sipkens, T. A., Olfert, J. S., & Rogak, S. N. (Under preparation). Inversion methods to determine two-dimensional aerosol mass-mobility distributions: Existing and novel Bayesian methods.][4]

[1]: https://doi.org/10.1016/j.jaerosci.2019.105484
[2]: https://doi.org/10.1080/02786826.2019.1680794
[3]: https://doi.org/10.1016/j.jaerosci.2017.09.012
[4]: N/A
