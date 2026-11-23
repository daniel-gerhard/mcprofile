# Changelog

## mcprofile 1.1-0

CRAN release: 2026-11-23

### Changes in results

- [`hoa()`](https://daniel-gerhard.github.io/mcprofile/reference/hoa.md)
  now computes the modified likelihood root correctly for contrasts with
  more than one nonzero coefficient, such as averages, trends or
  Tukey-type comparisons in models with correlated coefficients.
  Previously the nuisance information was taken from a sub-block of the
  covariance matrix, which is only correct for contrasts that select a
  single coefficient. Adjusted statistics, intervals and p-values change
  for such contrasts, in some cases by as much as the whole higher-order
  correction. Pairwise contrasts in cell-means models (for example the
  Dunnett example in the vignette) were nearly unaffected.

- [`hoa()`](https://daniel-gerhard.github.io/mcprofile/reference/hoa.md)
  also uses the working weights at the final constrained fit, which
  changes its results slightly (up to about 1e-4) for all contrasts.

- [`hoa()`](https://daniel-gerhard.github.io/mcprofile/reference/hoa.md)
  no longer discards every profile point for larger data sets.

### Bug fixes

- [`plot()`](https://rdrr.io/r/graphics/plot.default.html) failed for
  intervals from `confint(..., alternative = "less")`.

- [`expit()`](https://daniel-gerhard.github.io/mcprofile/reference/expit.mcpCI.md)
  is now an exported generic, so `expit(confint(...))` transforms
  intervals to the probability scale as documented.

- [`mcprofile()`](https://daniel-gerhard.github.io/mcprofile/reference/mcprofile.md)
  failed when all grid points lay on one side of the estimate.

- [`mcprofile()`](https://daniel-gerhard.github.io/mcprofile/reference/mcprofile.md)
  for `lm` objects failed when the data, weights or formula were defined
  inside a function, or when the data had the same name as an argument
  of
  [`mcprofile()`](https://daniel-gerhard.github.io/mcprofile/reference/mcprofile.md)
  (for example `CM`).

- [`confint()`](https://rdrr.io/r/stats/confint.html) now checks
  `adjust` and `alternative`. A misspelled `alternative` previously gave
  a confusing error, or could silently return an unrelated object called
  `ci` from the workspace as the intervals. Unambiguous abbreviations
  such as `alternative = "g"` are now accepted.

- [`summary()`](https://rdrr.io/r/base/summary.html) validates `adjust`
  and `alternative` the same way, and also accepts unambiguous
  abbreviations.

- [`mcprofile()`](https://daniel-gerhard.github.io/mcprofile/reference/mcprofile.md)
  now gives a clear error for models with aliased (`NA`) coefficients,
  and for a contrast matrix with the wrong number of columns.

- When a profile cannot be interpolated (for example with fewer than
  four profile points),
  [`confint()`](https://rdrr.io/r/stats/confint.html),
  [`summary()`](https://rdrr.io/r/base/summary.html) and
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) give one
  warning naming the affected contrasts, whose results are `NA`.
  Previously an error message was printed to the console without a
  warning.

### Other changes

- ggplot2 is now imported instead of attached, so
  [`library(mcprofile)`](https://daniel-gerhard.github.io/mcprofile/) no
  longer puts ggplot2 functions on the search path.
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) still returns
  ggplot objects; call
  [`library(ggplot2)`](https://ggplot2.tidyverse.org) to modify them,
  for example with `+ ylim(-5, 5)`.

- Profiling is faster: each constrained fit solves the quadratic program
  once per iteration, and the starting values on each side of the
  estimate are chosen better.

- Requires ggplot2 3.4.0 or later.
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) for intervals
  no longer uses the deprecated
  [`geom_errorbarh()`](https://ggplot2.tidyverse.org/reference/geom_linerange.html).

- The citation uses
  [`bibentry()`](https://rdrr.io/r/utils/bibentry.html).

- Added a test suite.
