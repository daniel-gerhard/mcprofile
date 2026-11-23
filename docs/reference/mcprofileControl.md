# mcprofile Control Arguments

Control arguments for the mcprofile function

## Usage

``` r
mcprofileControl(maxsteps = 10, alpha = 0.01, del = function(zmax) zmax/5)
```

## Arguments

- maxsteps:

  Maximum number of points to be used for profiling each parameter.

- alpha:

  Highest significance level allowed for the profile t-statistics
  (Bonferroni adjusted)

- del:

  Suggested change on the scale of the profile t-statistics. Default
  value chosen to allow profiling at about 10 parameter values.

## See also

[`mcprofile`](https://daniel-gerhard.github.io/mcprofile/reference/mcprofile.md)
