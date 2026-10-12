## log determinant of the nuisance information |j_ll| at a constrained fit,
## log det(J) + log(c' J^-1 c) with J = X'WX and c the constraint (up to a constant)
logdetjll <- function(object){
  ## working weights at the final coefficients (object$weights are from the start of the last IRLS step)
  fam <- object$family
  w <- object$prior.weights * fam$mu.eta(object$linear.predictors)^2 / fam$variance(object$fitted.values)
  xm <- object$X * sqrt(w)
  J <- crossprod(xm)
  A <- object$constr
  as.numeric(determinant(J, logarithm = TRUE)$modulus) + log(drop(A %*% solve(J, t(A))))
}
