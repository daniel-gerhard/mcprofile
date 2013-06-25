\name{aphidlight}
\alias{aphidlight}
\docType{data}
\title{Aphid attraction at different light intensities}
\description{The light intensity (mumol/m^2s) of green LED light
  should be found, which attracts Aphis fabae best. At each of 4
  replicates 20 aphids were put in a lightproof box with only one green
  LED at one end. All aphids that fly to the green light are caught and
  counted after a period of 5h. This procedure was replicated for 9
  increasing light intensities.}
\usage{data(aphidlight)}
\format{
  A data frame with 36 observations on the following 3 variables.
  \describe{
    \item{\code{light}}{a numeric vector denoting the concentration levels}
    \item{\code{black}}{a numeric vector with the number of aphids
      remaining in the box.}
    \item{\code{green}}{a numeric vector with the number of attracted aphids}
  }
}
\references{
Akyazi, G (2009): Zum Einfluss auf Lichtintensitaet und Lichtqualitaet (Hochleistungs-LEDs) auf das Verhalten von Aphis fabae. IPP MSc 19.
}
\keyword{datasets}
