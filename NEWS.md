# mcprofile 1.0-2

## Changes in results

* `hoa()` now computes the modified likelihood root correctly for contrasts
  with more than one nonzero coefficient, such as averages, trends or
  Tukey-type comparisons in models with correlated coefficients. Previously
  the nuisance information was taken from a sub-block of the covariance
  matrix, which is only correct for contrasts that select a single
  coefficient. Adjusted statistics, intervals and p-values change for such
  contrasts, in some cases by as much as the whole higher-order correction.
  Pairwise contrasts in cell-means models (for example the Dunnett example in
  the vignette) were nearly unaffected.

* `hoa()` also uses the working weights at the final constrained fit, which
  changes its results slightly (up to about 1e-4) for all contrasts.

* `hoa()` no longer discards every profile point for larger data sets.

## Bug fixes

* `plot()` failed for intervals from `confint(..., alternative = "less")`.

* `expit()` is now an exported generic, so `expit(confint(...))` transforms
  intervals to the probability scale as documented.

* `mcprofile()` failed when all grid points lay on one side of the estimate.

* `mcprofile()` for `lm` objects failed when the data, weights or formula
  were defined inside a function, or when the data had the same name as an
  argument of `mcprofile()` (for example `CM`).

* `confint()` now checks `adjust` and `alternative`. A misspelled
  `alternative` previously gave a confusing error, or could silently return
  an unrelated object called `ci` from the workspace as the intervals.
  Unambiguous abbreviations such as `alternative = "g"` are now accepted.

* `mcprofile()` now gives a clear error for models with aliased (`NA`)
  coefficients, and for a contrast matrix with the wrong number of columns.

## Other changes

* Profiling is faster: each constrained fit solves the quadratic program once
  per iteration, and the starting values on each side of the estimate are
  chosen better.

* Requires ggplot2 3.4.0 or later. `plot()` for intervals no longer uses the
  deprecated `geom_errorbarh()`.

* The citation uses `bibentry()`.

* Added a test suite.
