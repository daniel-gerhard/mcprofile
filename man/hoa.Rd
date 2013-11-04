\name{hoa}
\alias{hoa}
\title{Higher order asymptotics using the modified likelihood root}
\description{Transforms a signed root deviance profile to a modified
  likelihood root profile.}
\usage{
hoa(object, maxstat=10)
}
\arguments{
  \item{object}{An object of class mcprofile}
  \item{maxstat}{Limits the statistic to a maximum absolute value (default=10)}
}
\value{An object of class mcprofile with a hoa profile in the srdp slot.}

\author{Daniel Gerhard}

\seealso{\code{\link{mcprofile}}}

\examples{
#######################################
## cell transformation assay example ##
#######################################

str(cta)
## change class of cta$conc into factor
cta$concf <- factor(cta$conc, levels=unique(cta$conc))

ggplot(cta, aes(y=foci, x=concf)) + 
   geom_boxplot() +
   geom_dotplot(binaxis = "y", stackdir = "center", binwidth = 0.2) +    
   xlab("concentration")
   

# glm fit assuming a Poisson distribution for foci counts
# parameter estimation on the log link
# removing the intercept
fm <- glm(foci ~ concf-1, data=cta, family=poisson(link="log"))

### Comparing each dose to the control by Dunnett-type comparisons
# Constructing contrast matrix
library(multcomp)
CM <- contrMat(table(cta$concf), type="Dunnett")

# calculating signed root deviance profiles
(dmcp <- mcprofile(fm, CM))
# computing profiles for the modified likelihood root
hp <- hoa(dmcp)

plot(hp)

# comparing confidence intervals
confint(hp)
confint(dmcp)
}

\keyword{ misc }
