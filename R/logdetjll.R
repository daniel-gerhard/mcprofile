## log determinant of the nuisance information |j_ll| at a constrained fit,
## log det(J) + log(c' J^-1 c) with J = X'WX and c the constraint (up to a constant)
logdetjll <- function(object){
  xm <- object$X * sqrt(object$weights)
  J <- crossprod(xm)
  A <- object$constr
  as.numeric(determinant(J, logarithm = TRUE)$modulus) + log(drop(A %*% solve(J, t(A))))
}
