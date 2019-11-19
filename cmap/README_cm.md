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
also takes `n` as a second input, which can be used reduce the number of
colors in the colormap, while still respecting the color order.

It is also noted that the *deep*, *dense*, *matter*, and *tempo* colormaps
are reversed from their original order, such that the darker color is
always first.

### Included Colormaps

#### Sequantial

From mpl colormaps:

1. *viridis*
2. *inferno*
3. *plasma*
4. *magma*

From cmocean:

5. *thermal*
6. *haline*
7. *ice*
8. *deep*
9. *dense*
10. *matter*
11. *tempo*
12. *speed* - Yellow, green colormap

From colorbrewer2:

13. *YlGnBu*
14. *BuPu*
15. *RdPu*

#### Divergent colormaps

From cmocean:

1. *balance*
2. *delta*
3. *curl*

From colorbrewer2:

4. *PuOr*
5. *RdBu*
6. *PrGn*

#### Rainbow colormaps

1. *turbo* (dedicated source)

From cmocean:
2. *phase*
