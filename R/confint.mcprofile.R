#' Simultaneous Confidence Intervals for Multiple Contrast Profiles
#' 
#' Calculates simultaneous confidence intervals based on signed root deviance profiles from function \code{\link{mcprofile}}.
#' 
#' @param object An object of class mcprofile
#' @param parm Just ignore this...
#' @param level Simultaneous confidence level (1-alpha), default at 0.95
#' @param adjust a character string specifying the adjustment for multiplicity. "single-step" controlling the FWER utilising a multivariate normal- or t-distribution; "none" for comparison-wise error rate; "bonferroni" applying a Bonferroni correction.
#' @param alternative a character string specifying if two- or one-sided confidence intervals should be computed
#' @param ... ...
#' @return An object of class mcpCI
#' @seealso \code{\link[MASS]{confint.glm}}, \code{\link{mcprofile}}, \code{\link[multcomp]{confint.glht}}
#' @keywords htest

confint.mcprofile <-
function(object, parm, level=0.95, adjust=c("single-step","none","bonferroni"), alternative=c("two.sided","less","greater"), ...){
  adjust <- match.arg(adjust)
  alternative <- match.arg(alternative)
  CM <- object$CM
  df <- object$df
  cr <- NULL
  sdlist <- object$srdp  
  spl <- lapply(sdlist, function(x){
    x <- na.omit(x)
    try(interpSpline(x[,1], x[,2]))
  })
  
  if (adjust == "none" | nrow(CM) == 1){
    if (alternative == "two.sided") alpha <- (1-level)/2 else alpha <- 1-level
    if (is.null(df)) quant <- qnorm(1-alpha) else quant <- qt(1-alpha, df=df)
  }
  if (adjust == "bonferroni"){
    if (alternative == "two.sided") alpha <- (1-level)/2 else alpha <- 1-level
    if (is.null(df)) quant <- qnorm(1-alpha/nrow(CM)) else quant <- qt(1-alpha/nrow(CM), df=df)
  }
  if (adjust == "single-step" & nrow(CM) > 1){
    vc <- vcov(object$object)
    VC <- CM %*% vc %*% t(CM)
    d <- 1/sqrt(diag(VC))
    dd <- diag(d)
    cr <- dd %*% VC %*% dd
    if (alternative == "two.sided"){
      if (is.null(df)) quant <- qmvnorm(level, corr=cr, tail="both.tails")$quantile else quant <- qmvt(level, df=df, corr=cr, tail="both.tails")$quantile
    } else {
      if (is.null(df)) quant <- qmvnorm(level, corr=cr, tail="lower.tail")$quantile else quant <- qmvt(level, df=df, corr=cr, tail="lower.tail")$quantile
    }
  }
   
  if (alternative == "two.sided"){
    ci <- data.frame(t(sapply(spl, function(x, quant){
      pfun <- function(xc, obj, quant) predict(obj, xc)$y-quant
      upper <- try(uniroot(pfun, range(predict(x)$x), obj=x, quant=quant)$root, silent=TRUE)
      if (inherits(upper, "try-error")) upper <- NA
      lower <- try(uniroot(pfun, range(predict(x)$x), obj=x, quant=-quant)$root, silent=TRUE)
      if (inherits(lower, "try-error")) lower <- NA
      c(lower, upper)
    }, quant=quant)))
    names(ci) <- c("lower", "upper")
  }
  if (alternative == "less"){
    ci <- data.frame(sapply(spl, function(x, quant){
      pfun <- function(xc, obj, quant) predict(obj, xc)$y-quant
      upper <- try(uniroot(pfun, range(predict(x)$x), obj=x, quant=quant)$root, silent=TRUE)
      if (inherits(upper, "try-error")) upper <- NA
      cbind(c(upper))
    }, quant=quant))
    names(ci) <- "upper"
  }
  if (alternative == "greater"){
    ci <- data.frame(sapply(spl, function(x, quant){
      pfun <- function(xc, obj, quant) predict(obj, xc)$y-quant
      lower <- try(uniroot(pfun, range(predict(x)$x), obj=x, quant=-quant)$root, silent=TRUE)
      if (inherits(lower, "try-error")) lower <- NA
      cbind(c(lower))
    }, quant=quant))
    names(ci) <- "lower"
  }  

  out <- list()
  out$estimate <- data.frame(Estimate=CM %*% coefficients(object$object))
  out$confint <- ci
  out$cr <- cr
  out$CM <- CM
  out$quant <- quant
  out$alternative <- alternative
  out$level <- level
  out$adjust <- adjust
  class(out) <- "mcpCI"
  return(out)
}

