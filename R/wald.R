#' Calculate Wald-Profiles
#' 
#' Transforms a signed root deviance profile of a mcprofile object into a profile of Wald-type statistics
#' 
#' @param object An object of class mcprofile
#' 
#' @return An object of class mcprofile with a wald profile in the srdp slot.
#' 
#' @seealso \code{\link{mcprofile}}
#' 
#' @keywords misc
#' 
#' @examples 
#' #######################################
#' ## cell transformation assay example ##
#' #######################################
#' 
#' str(cta)
#' ## change class of cta$conc into factor
#' cta$concf <- factor(cta$conc, levels=unique(cta$conc))
#' 
#' ggplot(cta, aes(y=foci, x=concf)) + 
#'   geom_boxplot() +
#'   geom_dotplot(binaxis = "y", stackdir = "center", binwidth = 0.2) + 
#'   xlab("concentration")
#'   
#'   
#' # glm fit assuming a Poisson distribution for foci counts
#' # parameter estimation on the log link
#' # removing the intercept
#' fm <- glm(foci ~ concf-1, data=cta, family=poisson(link="log"))
#' 
#' ### Comparing each dose to the control by Dunnett-type comparisons
#' # Constructing contrast matrix
#' library(multcomp)
#' CM <- contrMat(table(cta$concf), type="Dunnett")
#' 
#' # calculating signed root deviance profiles
#' (dmcp <- mcprofile(fm, CM))
#' # computing profiles for the modified likelihood root
#' wp <- wald(dmcp)
#' 
#' plot(wp)
#' 
#' # comparing confidence intervals
#' confint(wp)
#' confint(dmcp)

wald <-
  function(object){
    srdp <- object$srdp
    est <- object$CM %*% coefficients(object$object)
    sde <- sqrt(diag(object$CM %*% vcov(object$object) %*% t(object$CM)))
    wsrdp <- lapply(1:length(srdp), function(i){
      srdpi <- srdp[[i]]
      b <- srdpi[,1]
      srdpi[,2] <- (b-est[i])/sde[i]
      srdpi 
    })
    object$srdp <- wsrdp
    return(object)
  }