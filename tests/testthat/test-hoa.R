## r* for a contrast c'beta, computed directly from its definition for a
## canonical GLM: r* = r + log(q / r) / r, q = Wald statistic * rho,
## rho = sqrt(|j_ll(MLE)| / |j_ll(constrained)|), |j_ll| ~ det(J) * c' J^-1 c
rstar_direct <- function(mc, cvec, weightfun) {
  fm <- mc$object
  X <- model.matrix(fm)
  jll <- function(beta) {
    J <- crossprod(X * sqrt(weightfun(drop(X %*% beta))))
    det(J) * drop(t(cvec) %*% solve(J) %*% cvec)
  }
  b <- mc$srdp[[1]][, 1]
  r <- mc$srdp[[1]][, 2]
  est <- sum(cvec * coef(fm))
  q <- (b - est) / sqrt(drop(t(cvec) %*% vcov(fm) %*% cvec))
  rho <- sqrt(jll(coef(fm)) / apply(mc$optpar[[1]], 1, jll))
  r + log(rho * q / r) / r
}

test_that("hoa matches the direct r* formula for a weighted contrast", {
  fm <- cta_fit()
  cvec <- c(1, -0.5, -0.5, 0, 0, 0, 0, 0)
  mc <- mcprofile(fm, matrix(cvec, 1), grid = list(grid_around(fm, cvec)))
  expect_equal(hoa(mc)$srdp[[1]][, 2], rstar_direct(mc, cvec, exp), tolerance = 1e-6)
})

test_that("hoa matches the direct r* formula for a single coefficient", {
  fm <- cta_fit()
  cvec <- c(1, 0, 0, 0, 0, 0, 0, 0)
  mc <- mcprofile(fm, matrix(cvec, 1), grid = list(grid_around(fm, cvec)))
  expect_equal(hoa(mc)$srdp[[1]][, 2], rstar_direct(mc, cvec, exp), tolerance = 1e-6)
})

test_that("hoa does not depend on how the model is parametrised", {
  set.seed(3)
  d <- data.frame(g = gl(4, 10), x = rnorm(40))
  d$y <- rbinom(40, 12, plogis(-0.3 + 0.3 * as.numeric(d$g) + 0.6 * d$x))
  fm <- glm(cbind(y, 12 - y) ~ g + x, data = d, family = binomial())
  X <- model.matrix(fm)
  for (cvec in list(c(0, 1, -1, 0, 0), c(1, 0, 0, 0.5, 2), c(0, 1, 1, 1, -1))) {
    ## same model, reparametrised so that the contrast is the first coefficient
    A <- rbind(cvec, diag(5)[-which(cvec != 0)[1], ])
    Xr <- X %*% solve(A)
    fr <- glm(cbind(d$y, 12 - d$y) ~ Xr - 1, family = binomial())
    grid <- list(grid_around(fm, cvec))
    h1 <- hoa(mcprofile(fm, matrix(cvec, 1), grid = grid))$srdp[[1]][, 2]
    h2 <- hoa(mcprofile(fr, matrix(c(1, 0, 0, 0, 0), 1), grid = grid))$srdp[[1]][, 2]
    expect_equal(h1, h2, tolerance = 1e-8)
  }
})

test_that("hoa keeps the profile for large samples", {
  set.seed(1)
  d <- data.frame(g = gl(4, 1000))
  d$y <- rpois(4000, 20)
  fm <- glm(y ~ g - 1, data = d, family = poisson())
  h <- hoa(mcprofile(fm, diag(4)[1:2, ]))
  expect_true(all(sapply(h$srdp, nrow) >= 20))
  expect_equal(unlist(h$adjestimates), unname(coef(fm)[1:2]), tolerance = 1e-3)
})
