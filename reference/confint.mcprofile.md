# Simultaneous Confidence Intervals for Multiple Contrast Profiles

Calculates simultaneous confidence intervals based on signed root
deviance profiles from function
[`mcprofile`](https://daniel-gerhard.github.io/mcprofile/reference/mcprofile.md).

## Usage

``` r
# S3 method for class 'mcprofile'
confint(
  object,
  parm,
  level = 0.95,
  adjust = c("single-step", "none", "bonferroni"),
  alternative = c("two.sided", "less", "greater"),
  ...
)
```

## Arguments

- object:

  An object of class mcprofile

- parm:

  Not used; intervals are computed for all contrasts in `object`.

- level:

  Simultaneous confidence level (1-alpha), default at 0.95

- adjust:

  a character string specifying the adjustment for multiplicity.
  "single-step" controlling the FWER utilising a multivariate normal- or
  t-distribution; "none" for comparison-wise error rate; "bonferroni"
  applying a Bonferroni correction.

- alternative:

  a character string specifying if two- or one-sided confidence
  intervals should be computed

- ...:

  Further arguments, currently ignored.

## Value

An object of class mcpCI

## See also

[`confint.glm`](https://rdrr.io/pkg/MASS/man/confint.html),
[`mcprofile`](https://daniel-gerhard.github.io/mcprofile/reference/mcprofile.md),
[`confint.glht`](https://rdrr.io/pkg/multcomp/man/methods.html)
