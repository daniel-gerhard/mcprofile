\name{mcprofileControl}
\alias{mcprofileControl}
\title{mcprofile Control Arguments}
\description{Control arguments for the mcprofile function}
\usage{
mcprofileControl(maxsteps=10, alpha=0.01, del=function(zmax) zmax/5)
}
\arguments{
  \item{maxsteps}{Maximum number of points to be used for profiling each
    parameter.}  
  \item{alpha}{Highest significance level allowed for the profile
    t-statistics (Bonferroni adjusted)}
  \item{del}{Suggested change on the scale of the profile t-statistics. Default value chosen to allow profiling at about 10 parameter values.}  
}
\author{Daniel Gerhard}

\seealso{\code{\link{mcprofile}}}

\keyword{ misc }
