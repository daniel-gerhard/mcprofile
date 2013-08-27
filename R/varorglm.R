varorglm <- function(object){
  require(quadprog)
  eta <- object$X %*% object$coefficients
  y <- object$y
  z <- as.vector(eta  + (y - object$family$linkinv(eta))/object$family$mu.eta(eta))
  w <- sqrt(object$weights)
  Kp <- as.vector(object$constr) != 0
  K0 <- as.vector(object$constr) == 0
  if (sum(Kp) == 1){
    xm <- object$X
    xm[,Kp] <- 0
    p <- ncol(xm) - 1
  } else {
    xm <- cbind(apply(object$X[,Kp,drop=FALSE], 1, sum), object$X[,K0,drop=FALSE])
    p <- sum(K0) + 1
  }  
  lmf <- lm.fit(xm * w, z * w)  
  chol2inv(lmf$qr$qr[1:p, 1:p, drop = FALSE])  
}
