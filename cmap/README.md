# Perceptually improved colormaps for MATLAB

*Last updated: November 19, 2019*

These `.mat` files contain the colormaps from four primary sources:

1. *cmocean* - Kristen M. Thyng, Chad A. Greene, Robert D. Hetland, Heather M. Zimmerle,
and Steven F. DiMarco. True colors of oceanography: Guidelines for effective
and accurate colormap selection. Oceanography, September 2016.  
http://dx.doi.org/10.5670/oceanog.2016.66 (More information is available at
https://matplotlib.org/cmocean/).

2. *matplotlib* - Colormaps designed by Stéfan van der Walt and
Nathaniel Smith. (More information is available at https://bids.github.io/colormap/).

3. *colorbrewer2* - Colormaps by Cynthia Brewer and Mark Harrower. (More information
  available at http://colorbrewer2.org/).

4. *turbo* - A. Mikhailov. Turbo, An Improved Rainbow Colormap for Visualization.
(More information is available at https://ai.googleblog.com/2019/08/turbo-improved-rainbow-colormap-for.html).

When loaded directly, the colormaps will appear as the variable `cm` in the
workspace. Otherwise `load_cmap` can be used to load the colormap specified
by a string, `str`, containing the colormap name. The function `load_cmap`
also takes `n` as a second input, which can be used reduce or increase (by interpolation)
the number of colors in the colormap, while still respecting the color order.

It is also noted that the *deep*, *dense*, *matter*, and *tempo* colormaps
are reversed from their original order, such that the darker color is
always first.

The colormaps, and swages indicating their color progression, are included below.

## Sequantial colormaps

#### From mpl colormaps:

![viridis](imgs/viridis.jpg) *viridis*

![inferno](imgs/inferno.jpg) *inferno*

![plasma](imgs/plasma.jpg) *plasma*

![magma](imgs/magma.jpg) *magma*

#### From cmocean:

![thermal](imgs/thermal.jpg) *thermal*

![haline](imgs/haline.jpg) *haline*

![ice](imgs/ice.jpg) *ice*

![deep](imgs/deep.jpg) *deep*

![dense](imgs/dense.jpg) *dense*

![matter](imgs/matter.jpg) *matter*

![tempo](imgs/tempo.jpg) *tempo*

![speed](imgs/speed.jpg) *speed*

#### From colorbrewer2:

![YlGnBu](imgs/YlGnBu.jpg) *YlGnBu*

![BuPu](imgs/BuPu.jpg) *BuPu*

![RdPu](imgs/RdPu.jpg) *RdPu*

## Divergent colormaps

#### From cmocean:

![balance](imgs/balance.jpg) *balance*

![delta](imgs/delta.jpg) *delta*

![curl](imgs/curl.jpg) *curl*

#### From colorbrewer2:

![PuOr](imgs/PuOr.jpg) *PuOr*

![RdBu](imgs/RdBu.jpg) *RdBu*

![PrGn](imgs/PrGn.jpg) *PrGn*

## Rainbow colormaps

![turbo](imgs/turbo.jpg) *turbo* (dedicated source)
