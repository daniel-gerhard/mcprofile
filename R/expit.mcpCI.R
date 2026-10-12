#' Inverse logit transformation of Confidence Intervals
#'
#' Inverse logit transformation of confidence interval estimates in mcpCI objects.
#'
#' @param x An object of class mcpCI, or a numeric vector, matrix or data.frame for the default method.
#'
#' @return An object of class mcpCI with transformed estimates. The default method returns \code{1/(1 + exp(-x))}.
#'
#' @family confidence interval transformations
#' @seealso \code{\link{exp}}, \code{\link{confint.mcprofile}}
#'
#' @keywords misc
#' @examples
#' ## proportion of green aphids at each light level
#' aphidlight$lightf <- factor(aphidlight$light)
#' fm <- glm(cbind(green, black) ~ lightf - 1, data = aphidlight,
#'           family = binomial(link = "logit"))
#'
#' ## one contrast per light level: the logit of each group's proportion
#' CM <- diag(nlevels(aphidlight$lightf))
#' rownames(CM) <- paste("light", levels(aphidlight$lightf))
#'
#' mc <- mcprofile(fm, CM)
#' (ci <- confint(mc))
#'
#' ## simultaneous confidence intervals on the probability scale
#' (pci <- expit(ci))
#' plot(pci)

expit <- function(x) UseMethod("expit")

#' @rdname expit.mcpCI
expit.default <- function(x) 1/(1 + exp(-x))

#' @rdname expit.mcpCI
expit.mcpCI <-
  function(x){
    x$estimate <- expit(x$estimate)
    x$confint <- expit(x$confint)
    return(x)
  }
