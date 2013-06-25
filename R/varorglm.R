varorglm <- function(object){
  require(quadprog)
  eta <- object$X %*% object$coefficients
  y <- object$y
  z <- as.vector(eta  + (y - object$family$linkinv(eta))/object$family$mu.eta(eta))
  w <- sqrt(object$weights)
  xm <- object$X
  lmf <- lm.fit(xm * w, z * w)
  p <- sum(apply(xm, 2, function(x) any(x != 0)))
  chol2inv(lmf$qr$qr[1:p, 1:p, drop = FALSE])  
}
