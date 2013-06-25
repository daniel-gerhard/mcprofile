glm_profiling <-
function(X, Y, W, etastart, O, fam, glmcontrol, est, OriginalDeviance, DispersionParameter, K, b){    
  z <- numeric(length(b))
  pmi <- matrix(nrow=length(b), ncol=length(est))
  obj <- list()
  Kest <- as.vector(K %*% est)
  signs <- rep(1, length(b))
  signs[Kest > b] <- -1

  for (i in 1:length(b)){
    fm <- try(orglm.fit(x=X, y=Y, weights=W, etastart=etastart, offset=O, family=fam, control=glmcontrol, constr=K, rhs=b[i], nec=1), silent=TRUE)
    if (class(fm)[1] != "try-error"){
      zz <- signs[i]*sqrt(abs(fm$deviance - OriginalDeviance)/DispersionParameter)
      z[i] <- zz
      pmi[i,] <- fm$coefficients
      obj[[i]] <- fm
    } else {
      z[i] <- NA
      obj[[i]] <- NA
    }
  }
  
  list(stats=cbind(b, z), param=pmi, object=obj)
}

