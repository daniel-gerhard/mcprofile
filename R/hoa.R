
hoa <- function(object, removecenter=0.6){
  model <- object$object
  CM <- object$CM
  sumobj <- summary(model)
  Kuv <- sumobj$cov.unscaled
  r <- lapply(object$srdp, function(sr) sr[, 2])
  q <- lapply(wald(object)$srdp, function(sr) sr[, 2])
  vc <- diag(CM %*% vcov(model) %*% t(CM))
  vcc <- 1 / vc
  j <- (1/det(Kuv))
  j.1 <- j/vcc
  
  j.0 <- lapply(object$cobject, function(vm) sapply(vm, function(x) {
    xt <- try(1/det(mcprofile:::varorglm(x)), silent = TRUE)
    if (class(xt)[1] == "try-error") 
      xt <- NA
    return(xt)
  }))
  rho <- lapply(1:length(j.0), function(i) sqrt(j.1[i]/j.0[[i]]))
  adjr <- lapply(1:length(rho), function(i) r[[i]] + log(rho[[i]] * q[[i]]/r[[i]])/r[[i]])
  adjr <- lapply(adjr, function(x){
    x[is.infinite(x)] <- NA
    return(x)
  })
  
  if (!is.null(removecenter)){
    est <- as.vector(CM %*% model$coefficients)
    lwr <- est - removecenter*sqrt(vc)
    upr <- est + removecenter*sqrt(vc)
  }
  
  rsrdp <- lapply(1:length(adjr), function(i) {
    srdpi <- object$srdp[[i]]
    b <- srdpi[, 1]
    rstar <- adjr[[i]]    
    if (!is.null(removecenter)){
      rstar[b > lwr[i] & b < upr[i]] <- NA
    }
    srdpi[, 2] <- rstar
    return(na.omit(srdpi))
  })  
  
  rest <- sapply(rsrdp, function(x) {
    x <- na.omit(x)
    ispl <- try(interpSpline(x[, 1], x[, 2]), silent = TRUE)
    pfun <- function(xc, obj) predict(obj, xc)$y
    try(uniroot(pfun, range(predict(ispl)$x), obj = ispl)$root, 
        silent = TRUE)
  })
  object$srdp <- rsrdp
  object$adjestimates <- rest
  return(object)        
}


             
             