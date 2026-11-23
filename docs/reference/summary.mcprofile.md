# Multiple Testing of General Hypotheses

Multiple contrast testing based on signed root deviance profiles.

## Usage

``` r
# S3 method for class 'mcprofile'
summary(
  object,
  margin = 0,
  adjust = "single-step",
  alternative = c("two.sided", "less", "greater"),
  ...
)
```

## Arguments

- object:

  an object of class mcprofile

- margin:

  test margin, specifying the right hand side of the hypotheses.

- adjust:

  a character string specifying the adjustment for multiplicity.
  "single-step" controlling the FWER utilizing a multivariate normal- or
  t-distribution; "none" for comparison-wise error rate, or any other
  method provided by
  [`p.adjust`](https://rdrr.io/r/stats/p.adjust.html).

- alternative:

  a character string specifying the alternative hypothesis.

- ...:

  ...

## Value

An object of class mcpSummary

## See also

[`mcprofile`](https://daniel-gerhard.github.io/mcprofile/reference/mcprofile.md),
[`summary.glht`](https://rdrr.io/pkg/multcomp/man/methods.html)
