#' Multiple Testing of General Hypotheses
#' 
#' Multiple contrast testing based on signed root deviance profiles.
#' 
#' @param object an object of class mcprofile
#' @param margin test margin, specifying the right hand side of the hypotheses.
#' @param adjust a character string specifying the adjustment for multiplicity. "single-step" controlling the FWER utilizing a multivariate normal- or t-distribution; "none" for comparison-wise error rate, or any other method provided by \code{\link{p.adjust}}.
#' @param alternative a character string specifying the alternative hypothesis.
#' @param ... ...
#' 
#' @return An object of class mcpSummary
#' 
#' @seealso \code{\link{mcprofile}}, \code{\link[multcomp]{summary.glht}}
#' 
#' @keywords htest

summary.mcprofile <-
function(object, margin=0, adjust=c("single-step", "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"), alternative=c("two.sided","less","greater"), ...){
  adjust <- match.arg(adjust)
  alternative <- match.arg(alternative)
  CM <- object$CM
  est <- coefficients(object$object)
  df <- object$df
  ##### Find stats
  sdlist <- object$srdp
  spl <- lapply(sdlist, function(x){
    x <- na.omit(x)
    try(interpSpline(x[,1], x[,2]))
  })
  ptest <- function(spl, delta){
    pst <- try(predict(spl, delta)$y, silent=TRUE)
    if (inherits(pst, "try-error")) NA else pst   
  }
  stat <- sapply(spl, function(x) -1*ptest(x, delta=margin))
  naid <- is.na(stat)
  stat[is.na(stat)] <- 0
  if (is.null(df)){
    switch(alternative, less = {
      praw <- pnorm(stat, lower.tail=TRUE)
    }, greater = {
      praw <- pnorm(stat, lower.tail=FALSE)
    }, two.sided = {
      praw <- pmin(1, pnorm(abs(stat), lower.tail=FALSE)*2)
    })
  } else {
    switch(alternative, less = {
      praw <- pt(stat, df=df, lower.tail=TRUE)
    }, greater = {
      praw <- pt(stat, df=df, lower.tail=FALSE)
    }, two.sided = {
      praw <- pmin(1, pt(abs(stat), df=df, lower.tail=FALSE)*2)
    })
  }
  if (length(praw) > 1){
    if (adjust != "single-step") padj <- p.adjust(praw, method=adjust)
    pfct <- function(q) {
      switch(alternative, two.sided = {
       low <- rep(-abs(q), dim)
       upp <- rep(abs(q), dim)
      }, less = {
       low <- rep(q, dim)
       upp <- rep(Inf, dim)
      }, greater = {
       low <- rep(-Inf, dim)
       upp <- rep(q, dim)
      })
      if (is.null(df)){
       pmvnorm(lower = low, upper = upp, corr = cr)
      } else {
       pmvt(lower = low, upper = upp, df=df, corr = cr)
      }
    }
    if (adjust == "single-step"){
      vc <- vcov(object$object)
      VC <- CM %*% vc %*% t(CM)
      d <- 1/sqrt(diag(VC))
      dd <- diag(d)
      cr <- dd %*% VC %*% dd
      dim <- ncol(cr)
      padj <- numeric(length(stat))
      error <- 0
      for (i in 1:length(stat)) {
       tmp <- pfct(stat[i])
       if (error < attr(tmp, "error")) error <- attr(tmp, "error")
       padj[i] <- tmp
      }
      padj <- 1 - padj
      attr(padj, "error") <- error
    }
  } else { padj <- praw }
  stat[naid] <- NA
  padj[naid] <- NA

  out <- list()
  out$CM <- CM
  out$statistic <- stat
  out$p.values <- padj
  out$margin <- margin
  out$estimate = est
  out$alternative <- alternative
  out$adjust <- adjust
  out$df <- df
  class(out) <- "mcpSummary"
  return(out)
}

