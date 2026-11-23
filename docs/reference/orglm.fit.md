# Fitting Order-Restricted Generalized Linear Models

`orglm.fit` is used to fit generalized linear models with restrictions
on the parameters, specified by giving a description of the linear
predictor, a description of the error distribution, and a description of
a matrix with linear constraints. The `quadprog` package is used to
apply linear constraints on the parameter vector.

## Usage

``` r
orglm.fit(
  x,
  y,
  weights = rep(1, nobs),
  start = NULL,
  etastart = NULL,
  mustart = NULL,
  offset = rep(0, nobs),
  family = gaussian(),
  control = list(),
  intercept = TRUE,
  constr,
  rhs,
  nec
)
```

## Arguments

- x:

  is a design matrix of dimension `n * p`

- y:

  is a vector of observations of length `n`

- weights:

  an optional vector of ‘prior weights’ to be used in the fitting
  process. Should be `NULL` or a numeric vector.

- start:

  starting values for the parameters in the linear predictor.

- etastart:

  starting values for the linear predictor.

- mustart:

  starting values for the vector of means.

- offset:

  this can be used to specify an *a priori* known component to be
  included in the linear predictor during fitting. This should be `NULL`
  or a numeric vector of length equal to the number of cases. One or
  more [`offset`](https://rdrr.io/r/stats/offset.html) terms can be
  included in the formula instead or as well, and if more than one is
  specified their sum is used. See
  [`model.offset`](https://rdrr.io/r/stats/model.extract.html).

- family:

  a description of the error distribution and link function to be used
  in the model. This can be a character string naming a family function,
  a family function or the result of a call to a family function. (See
  [`family`](https://rdrr.io/r/stats/family.html) for details of family
  functions.)

- control:

  a list of parameters for controlling the fitting process. For
  `orglm.fit` this is passed to
  [`glm.control`](https://rdrr.io/r/stats/glm.control.html).

- intercept:

  logical. Should an intercept be included in the *null* model?

- constr:

  a matrix with linear constraints. The columns of this matrix should
  correspond to the columns of the design matrix.

- rhs:

  right hand side of the linear constraint formulation. A numeric vector
  with a length corresponding to the rows of `constr`.

- nec:

  Number of equality constrints. The first `nec` constraints defined in
  `constr` are treated as equality constraints; the remaining ones are
  inequality constraints.

## Value

An object of class `"glm"` is a list containing at least the following
components:

- coefficients:

  a named vector of coefficients

- residuals:

  the *working* residuals, that is the residuals in the final iteration
  of the IWLS fit. Since cases with zero weights are omitted, their
  working residuals are `NA`.

- fitted.values:

  the fitted mean values, obtained by transforming the linear predictors
  by the inverse of the link function.

- rank:

  the numeric rank of the fitted linear model.

- family:

  the [`family`](https://rdrr.io/r/stats/family.html) object used.

- linear.predictors:

  the linear fit on link scale.

- deviance:

  up to a constant, minus twice the maximized log-likelihood. Where
  sensible, the constant is chosen so that a saturated model has
  deviance zero.

- null.deviance:

  The deviance for the null model, comparable with `deviance`. The null
  model will include the offset, and an intercept if there is one in the
  model. Note that this will be incorrect if the link function depends
  on the data other than through the fitted mean: specify a zero offset
  to force a correct calculation.

- iter:

  the number of iterations of IWLS used.

- weights:

  the *working* weights, that is the weights in the final iteration of
  the IWLS fit.

- prior.weights:

  the weights initially supplied, a vector of `1`s if none were.

- df.residual:

  the residual degrees of freedom of the unconstrained model.

- df.null:

  the residual degrees of freedom for the null model.

- y:

  if requested (the default) the `y` vector used. (It is a vector even
  for a binomial model.)

- converged:

  logical. Was the IWLS algorithm judged to have converged?

- boundary:

  logical. Is the fitted value on the boundary of the attainable values?

## Details

Non-`NULL` `weights` can be used to indicate that different observations
have different dispersions (with the values in `weights` being inversely
proportional to the dispersions); or equivalently, when the elements of
`weights` are positive integers \\w_i\\, that each response \\y_i\\ is
the mean of \\w_i\\ unit-weight observations. For a binomial GLM prior
weights are used to give the number of trials when the response is the
proportion of successes: they would rarely be used for a Poisson GLM. If
more than one of `etastart`, `start` and `mustart` is specified, the
first in the list will be used. It is often advisable to supply starting
values for a [`quasi`](https://rdrr.io/r/stats/family.html) family, and
also for families with unusual links such as `gaussian("log")`. For the
background to warning messages about ‘fitted probabilities numerically 0
or 1 occurred’ for binomial GLMs, see Venables & Ripley (2002, pp.
197–8).

## References

- Dobson, A. J. (1990) *An Introduction to Generalized Linear Models.*
  London: Chapman and Hall.

- Hastie, T. J. and Pregibon, D. (1992) *Generalized linear models.*
  Chapter 6 of *Statistical Models in S* eds J. M. Chambers and T. J.
  Hastie, Wadsworth & Brooks/Cole.

- McCullagh P. and Nelder, J. A. (1989) *Generalized Linear Models.*
  London: Chapman and Hall.

- Venables, W. N. and Ripley, B. D. (2002) *Modern Applied Statistics
  with S.* New York: Springer.

## See also

[`glm`](https://rdrr.io/r/stats/glm.html),
[`solve.QP`](https://rdrr.io/pkg/quadprog/man/solve.QP.html)

## Author

Modification of the original glm.fit by Daniel Gerhard. The original R
implementation of `glm` was written by Simon Davies working for Ross
Ihaka at the University of Auckland, but has since been extensively
re-written by members of the R Core team. The design was inspired by the
S function of the same name described in Hastie & Pregibon (1992).
