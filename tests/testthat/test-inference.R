## For a Gaussian model the profile is linear in the parameter, so unadjusted
## profile intervals and tests must reproduce the classical t-based results.

gauss_setup <- function() {
  set.seed(10)
  d <- data.frame(y = rnorm(24, rep(c(1, 2, 4), each = 8)), g = gl(3, 8))
  fm <- lm(y ~ g, data = d)
  list(fm = fm, mc = mcprofile(fm, diag(3)))
}

test_that("unadjusted Gaussian intervals equal confint.lm", {
  s <- gauss_setup()
  ci <- confint(s$mc, adjust = "none")$confint
  expect_equal(unname(as.matrix(ci)), unname(confint(s$fm)), tolerance = 1e-6)
})

test_that("unadjusted Gaussian p-values equal the t-tests of summary.lm", {
  s <- gauss_setup()
  p <- summary(s$mc, adjust = "none")$p.values
  expect_equal(p, unname(coef(summary(s$fm))[, 4]), tolerance = 1e-6)
})

test_that("Bonferroni intervals use the Bonferroni-adjusted level", {
  s <- gauss_setup()
  ci <- confint(s$mc, adjust = "bonferroni")$confint
  expect_equal(unname(as.matrix(ci)), unname(confint(s$fm, level = 1 - 0.05 / 3)), tolerance = 1e-6)
})

test_that("one-sided intervals have a bound on the correct side only", {
  mc <- mcprofile(cta_fit(), dunnett(8))
  set.seed(1); two <- confint(mc)$confint
  set.seed(1); less <- confint(mc, alternative = "less")$confint
  set.seed(1); greater <- confint(mc, alternative = "greater")$confint
  expect_named(less, "upper")
  expect_named(greater, "lower")
  ## one-sided limits are tighter than the two-sided ones
  expect_true(all(less$upper < two$upper))
  expect_true(all(greater$lower > two$lower))
})

test_that("one-sided unadjusted p-values add up to one", {
  mc <- mcprofile(cta_fit(), dunnett(8))
  pl <- summary(mc, adjust = "none", alternative = "less")$p.values
  pg <- summary(mc, adjust = "none", alternative = "greater")$p.values
  expect_equal(pl + pg, rep(1, 7))
})

test_that("single-step intervals agree with multcomp for a Gaussian model", {
  skip_if_not_installed("multcomp")
  s <- gauss_setup()
  K <- rbind(c(0, 1, 0), c(0, 0, 1), c(0, -1, 1))
  set.seed(1)
  ci <- confint(mcprofile(s$fm, K))$confint
  set.seed(1)
  ref <- confint(multcomp::glht(s$fm, linfct = K))$confint
  expect_equal(unname(as.matrix(ci)), unname(ref[, 2:3]), tolerance = 1e-3)
})

test_that("confint rejects invalid alternative and adjust, without using a stray 'ci'", {
  mc <- mcprofile(cta_fit(), dunnett(8))
  ci <- data.frame(lower = 1:7, upper = 8:14)
  expect_error(confint(mc, alternative = "two-sided"), "should be one of")
  expect_error(confint(mc, adjust = "holm"), "should be one of")
})

test_that("confint accepts unambiguous abbreviations", {
  mc <- mcprofile(cta_fit(), dunnett(8))
  set.seed(1); a <- confint(mc, alternative = "g", adjust = "bonf")
  set.seed(1); b <- confint(mc, alternative = "greater", adjust = "bonferroni")
  expect_identical(a, b)
})

test_that("summary rejects an invalid alternative", {
  mc <- mcprofile(cta_fit(), dunnett(8))
  expect_error(summary(mc, alternative = "two-sided"), "alternative has to be one of")
})

test_that("an interval bound outside the grid is NA", {
  fm <- cta_fit()
  cvec <- dunnett(8)[7, , drop = FALSE]
  est <- sum(cvec * coef(fm))
  mc <- mcprofile(fm, cvec, grid = list(seq(est + 0.1, est + 2, length = 8)))
  ci <- confint(mc)$confint
  expect_true(is.na(ci$lower))
  expect_false(is.na(ci$upper))
})
