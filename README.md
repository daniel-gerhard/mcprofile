
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)

------------------------------------------------------------------------

[![minimal R version](https://img.shields.io/badge/R%3E%3D-3.1.0-6666ff.svg)](https://cran.r-project.org/) [![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/mcprofile)](https://cran.r-project.org/package=mcprofile) [![Downloads](https://cranlogs.r-pkg.org/badges/mcprofile)](https://cranlogs.r-pkg.org/)

------------------------------------------------------------------------

[![Last-changedate](https://img.shields.io/badge/last%20change-2017--12--17-yellowgreen.svg)](/commits/master)

mcprofile
=========

Overview
--------

mcprofile allows you to calculate signed root deviance profiles for linear combinations of parameters in a glm. Based on these profiles, multiplicity adjusted p-values and simultaneous confidence intervals controlling the family-wise error rate are calculated.

Installation
------------

``` r
# You can install mcprofile from CRAN by
install.packages("mcprofile")

# Or the the development version from GitHub:
# install.packages("devtools")
devtools::install_github("daniel-gerhard/mcprofile")
```
