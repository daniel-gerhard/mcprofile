\name{summary.mcprofile}
\alias{summary.mcprofile}
\title{Multiple Testing of General Hypotheses}
\description{Multiple contrast testing based on signed root deviance profiles.}
\usage{
\method{summary}{mcprofile}(object, margin = 0, adjust = "single-step", 
	alternative = c("two.sided", "less", "greater"), ...)
}
\arguments{
  \item{object}{an object of class mcprofile}
  \item{margin}{test margin, specifying the right hand side of the hypotheses.}
  \item{adjust}{a character string specifying the adjustment for
    multiplicity. "single-step" controlling the FWER utilizing a
    multivariate normal- or t-distribution; "none" for comparison-wise
    error rate, or any other method provided by \code{\link{p.adjust}}}
  \item{alternative}{a character string specifying the alternative hypothesis.}
  \item{\dots}{...}
}
\value{An object of class mcpSummary}
\author{Daniel Gerhard}

\seealso{\code{\link{mcprofile}}, \code{\link[multcomp]{summary.glht}}}

\keyword{ htest }
