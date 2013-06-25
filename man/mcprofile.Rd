\name{mcprofile}
\alias{mcprofile}
\alias{mcprofile.lm}
\alias{mcprofile.glm}

\title{Construction of Multiple Contrast Profiles}

\description{Calculates signed root deviance profiles given a \code{\link{glm}} or \code{\link{lm}} object. The profiled parameters of interest are defined by providing a contrast matrix.}

\usage{
mcprofile(object, CM, control = mcprofileControl(), grid=NULL)
\method{mcprofile}{lm}(object, CM, control=mcprofileControl(), grid=NULL)
\method{mcprofile}{glm}(object, CM, control=mcprofileControl(), grid=NULL)
}

\arguments{
  \item{object}{An object of class \code{\link{glm}} or \code{\link{lm}}}
  \item{CM}{A contrast matrix for the definition of parameter linear
  combinations (\code{CM \%*\% coefficients(object)}). The number of columns should be equal to the number of estimated parameters. Providing row names is recommendable.}
  \item{control}{A list with control arguments. See \code{\link{mcprofileControl}}.}
  \item{grid}{A matrix or list with profile support coordinates. Each
  column of the matrix or slot in a list corresponds to a row in the
  contrast matrix, each row of the grid matrix or element of a numeric
  vector in each list slot corresponds to a candidate of the contrast
  parameter. If NULL (default), a grid is found automatically similar to
  function \code{\link[MASS]{profile.glm}}.}
}

\value{An object of class mcprofile. The slot \code{srdp} contains the
  profiled signed root deviance statistics. The \code{optpar} slot
  contains a matrix with profiled parameter estimates.}

\details{The profiles are calculates separately for each row of the
  contrast matrix. The profiles are calculated by constrained IRWLS
  optimization, implemented in function \code{orglm}, using the quadratic programming algorithm of package \code{quadprog}.
}

\author{Daniel Gerhard}

\seealso{\code{\link[MASS]{profile.glm}}, \code{\link[multcomp]{glht}},
  \code{\link[multcomp]{contrMat}}, \code{\link{confint.mcprofile}},
  \code{\link{summary.mcprofile}}, \code{\link{orglm.fit}}, \code{\link[quadprog]{solve.QP}}}

\examples{
#######################################
## cell transformation assay example ##
#######################################

str(cta)
boxplot(foci ~ conc, cta, xlab="concentration", col="lightgreen")
## change class of cta$concentration into factor
cta$concf <- factor(cta$conc, levels=unique(cta$conc))

# glm fit assuming a Poisson distribution for foci counts
# parameter estimation on the log link
# estimating marginal means by removing the intercept
fm <- glm(foci ~ concf-1, data=cta, poisson(link="log"))

### Comparing each dose to the control by Dunnett-type comparisons
# Constructing contrast matrix
CM <- contrMat(table(cta$concf), type="Dunnett")

# calculating signed root deviance profiles
(dmcp <- mcprofile(fm, CM))
# plot profiles
plot(dmcp)
# confidence intervals
(ci <- confint(dmcp))
plot(ci)
}

\keyword{ misc }
