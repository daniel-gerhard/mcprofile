\name{hoa}
\alias{hoa}
\title{Higher order approximation (Test-Version)}
\description{Transforms a signed root deviance profile to a modified
  deviance profile. This is only a test version. Does not use the
  correct Fisher information matrix...}
\usage{
hoa(object, removecenter=0.6)
}
\arguments{
  \item{object}{An object of class mcprofile}
  \item{removecenter}{omits estimated statistics in the neighborhood around the estimates. default is est +- 0.6 * std.err}
}
\value{An object of class mcprofile with a hoa profile in the srdp slot.}

\author{Daniel Gerhard}

\seealso{\code{\link{mcprofile}}}

\keyword{ misc }
