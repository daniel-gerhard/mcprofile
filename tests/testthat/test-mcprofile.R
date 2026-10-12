test_that("mcprofile returns one profile per contrast", {
  fm <- cta_fit()
  mc <- mcprofile(fm, dunnett(8))
  expect_s3_class(mc, "mcprofile")
  expect_length(mc$srdp, 7)
  expect_true(all(sapply(mc$srdp, ncol) == 2))
  expect_false(anyNA(unlist(mc$srdp)))
})

test_that("profiles are signed: negative below the estimate, positive above", {
  fm <- cta_fit()
  cvec <- dunnett(8)[3, ]
  b <- grid_around(fm, cvec)
  z <- mcprofile(fm, matrix(cvec, 1), grid = list(b))$srdp[[1]][, 2]
  expect_true(all(z[b < sum(cvec * coef(fm))] < 0))
  expect_true(all(z[b > sum(cvec * coef(fm))] > 0))
})

test_that("Gaussian profiles equal the Wald statistic exactly", {
  set.seed(1)
  d <- data.frame(y = rnorm(30, rep(0:2, each = 10)), g = gl(3, 10))
  fm <- glm(y ~ g, data = d)
  mc <- mcprofile(fm, diag(3))
  w <- wald(mc)
  for (i in 1:3) expect_equal(mc$srdp[[i]][, 2], w$srdp[[i]][, 2], tolerance = 1e-6)
})

test_that("degrees of freedom are stored only when the dispersion is estimated", {
  expect_null(mcprofile(cta_fit(), dunnett(8))$df)
  expect_equal(mcprofile(cta_fit(quasipoisson()), dunnett(8))$df, 72)
})

test_that("lm and the equivalent Gaussian glm give the same profiles", {
  set.seed(2)
  d <- data.frame(y = rnorm(20), g = gl(2, 10))
  expect_equal(mcprofile(lm(y ~ g, data = d), diag(2))$srdp,
               mcprofile(glm(y ~ g, data = d), diag(2))$srdp)
})

test_that("mcprofile.lm finds data local to a function", {
  f <- function() {
    d <- data.frame(y = rnorm(20), g = gl(2, 10), w = runif(20, 0.5, 2))
    fo <- y ~ g
    m <- lm(fo, data = d, weights = w, subset = 1:18)
    mcprofile(m, diag(2))
  }
  expect_s3_class(f(), "mcprofile")
})

test_that("mcprofile.lm is not confused by data named like its own arguments", {
  set.seed(3)
  CM <- data.frame(y = rnorm(20, rep(c(0, 3), each = 10)), g = gl(2, 10))
  ref <- mcprofile(glm(y ~ g, data = CM), diag(2))
  expect_equal(mcprofile(lm(y ~ g, data = CM), diag(2))$srdp, ref$srdp)
})

test_that("grids on one side of the estimate work and match independent fits", {
  fm <- cta_fit()
  cvec <- dunnett(8)[6, , drop = FALSE]
  est <- sum(cvec * coef(fm))
  X <- model.matrix(fm)
  cold <- function(b) sapply(b, function(bb) {
    f <- orglm.fit(X, fm$y, etastart = fm$linear.predictors, family = poisson(),
                   control = fm$control, constr = cvec, rhs = bb, nec = 1)
    sign(bb - est) * sqrt(abs(f$deviance - fm$deviance))
  })
  for (b in list(seq(est + 0.1, est + 2, length = 6),
                 seq(est - 2, est - 0.1, length = 6),
                 c(est - 0.5, seq(est + 0.1, est + 1, length = 4)))) {
    z <- mcprofile(fm, cvec, grid = list(b))$srdp[[1]][, 2]
    expect_equal(z, cold(b), tolerance = 1e-8)
  }
})

test_that("invalid input gives informative errors", {
  fm <- cta_fit()
  expect_error(mcprofile(fm, diag(3)), "Equal number")
  expect_error(mcprofile(fm, dunnett(8)[1:2, ], grid = list(1:3)), "differ")
  expect_error(mcprofile(fm, dunnett(8)[1, , drop = FALSE], grid = list(c(2, 1, 3))), "increasing")
  expect_error(mcprofile(fm, dunnett(8)[1, , drop = FALSE], grid = list(1)), "At least 2")
})

test_that("models with aliased coefficients are rejected", {
  set.seed(4)
  d <- data.frame(g = gl(3, 10), y = rpois(30, 5))
  d$g2 <- d$g
  fm <- glm(y ~ g + g2, data = d, family = poisson())
  expect_error(mcprofile(fm, diag(5)), "aliased \\(NA\\) coefficients: g22, g23")
})
