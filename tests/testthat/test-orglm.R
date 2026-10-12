test_that("an equality-constrained fit equals the refit with the constraint as offset", {
  set.seed(5)
  d <- data.frame(x = rnorm(40), z = rnorm(40))
  d$y <- rpois(40, exp(0.5 + 0.4 * d$x - 0.3 * d$z))
  X <- model.matrix(~ x + z, d)
  ## fix the x coefficient at 0.2
  f <- orglm.fit(X, d$y, family = poisson(), constr = matrix(c(0, 1, 0), 1), rhs = 0.2, nec = 1)
  ref <- glm(y ~ z, offset = 0.2 * x, data = d, family = poisson())
  expect_equal(unname(f$coefficients[c(1, 3)]), unname(coef(ref)), tolerance = 1e-6)
  expect_equal(f$deviance, deviance(ref), tolerance = 1e-8)
})

test_that("equality and inequality constraints are met", {
  d <- cta_data()
  X <- model.matrix(~ concf - 1, d)
  A <- rbind(c(1, -1, 0, 0, 0, 0, 0, 0),    # mu1 == mu2
             c(0, 0, -1, 1, 0, 0, 0, 0),    # mu4 - mu3 >= 0.5
             c(0, 0, 0, 0, 0, 0, 1, -1))    # mu7 - mu8 >= 0.3
  f <- orglm.fit(X, d$foci, family = poisson(), constr = A, rhs = c(0, 0.5, 0.3), nec = 1)
  Ab <- drop(A %*% f$coefficients)
  expect_equal(Ab[1], 0, tolerance = 1e-8)
  expect_true(Ab[2] >= 0.5 - 1e-8)
  expect_true(Ab[3] >= 0.3 - 1e-8)
})

test_that("an unconstrained-optimum constraint reproduces the glm fit", {
  fm <- cta_fit()
  cvec <- matrix(dunnett(8)[2, ], 1)
  f <- orglm.fit(model.matrix(fm), fm$y, family = poisson(), constr = cvec,
                 rhs = sum(cvec * coef(fm)), nec = 1)
  expect_equal(unname(f$coefficients), unname(coef(fm)), tolerance = 1e-6)
})
