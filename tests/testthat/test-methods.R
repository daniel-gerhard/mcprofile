test_that("plot methods return ggplot objects for all alternatives", {
  mc <- mcprofile(cta_fit(), dunnett(8))
  expect_s3_class(plot(mc), "ggplot")
  for (a in c("two.sided", "less", "greater")) {
    p <- plot(confint(mc, alternative = a))
    expect_s3_class(p, "ggplot")
    expect_silent(ggplot2::ggplot_build(p))
  }
})

test_that("plots of one-sided intervals are open on the correct side", {
  mc <- mcprofile(cta_fit(), dunnett(8))
  lo <- ggplot2::layer_data(plot(confint(mc, alternative = "less")))
  gr <- ggplot2::layer_data(plot(confint(mc, alternative = "greater")))
  expect_true(all(lo$xmin == -Inf) && all(is.finite(lo$xmax)))
  expect_true(all(gr$xmax == Inf) && all(is.finite(gr$xmin)))
})

test_that("expit transforms intervals to the probability scale", {
  d <- aphidlight
  d$lightf <- factor(d$light)
  fm <- glm(cbind(green, black) ~ lightf - 1, data = d, family = binomial())
  set.seed(1)
  ci <- confint(mcprofile(fm, diag(nlevels(d$lightf))))
  pci <- expit(ci)
  expect_s3_class(pci, "mcpCI")
  expect_equal(as.matrix(pci$estimate), plogis(as.matrix(ci$estimate)), ignore_attr = TRUE)
  expect_equal(as.matrix(pci$confint), plogis(as.matrix(ci$confint)), ignore_attr = TRUE)
  expect_true(all(pci$confint >= 0 & pci$confint <= 1))
})

test_that("expit default and exp methods work", {
  expect_equal(expit(c(-Inf, 0, Inf)), c(0, 0.5, 1))
  set.seed(1)
  ci <- confint(mcprofile(cta_fit(), dunnett(8)))
  eci <- exp(ci)
  expect_equal(as.matrix(eci$confint), exp(as.matrix(ci$confint)), ignore_attr = TRUE)
})

test_that("print methods produce output", {
  mc <- mcprofile(cta_fit(), dunnett(8))
  expect_output(print(mc), "Multiple Contrast Profiles")
  expect_output(print(confint(mc)), "Confidence Intervals")
  expect_output(print(summary(mc)), "Multiple Testing")
})
