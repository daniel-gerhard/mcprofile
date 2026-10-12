## shared data and contrast matrices for the tests

cta_data <- function() {
  d <- cta
  d$concf <- factor(d$conc, levels = unique(d$conc))
  d
}

cta_fit <- function(family = poisson()) {
  glm(foci ~ concf - 1, data = cta_data(), family = family)
}

## Dunnett-type contrasts (each group vs the first) for k groups
dunnett <- function(k) {
  cm <- cbind(-1, diag(k - 1))
  rownames(cm) <- paste0("C", 2:k, "-C1")
  cm
}

## grid of n points around the estimate of contrast cvec, in standard errors
grid_around <- function(fm, cvec, at = c(-3, -2, -1, -0.5, 0.5, 1, 2, 3)) {
  est <- sum(cvec * coef(fm))
  se <- sqrt(drop(t(cvec) %*% vcov(fm) %*% cvec))
  est + at * se
}
