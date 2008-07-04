\name{sat.act}
\alias{sat.act}
\docType{data}
\title{3 Measures of ability: SATV, SATQ, ACT}
\description{Self reported scores on the SAT Verbal, SAT Quantitative and ACT   were collected as part of the Synthetic Aperture Personality Assessment (SAPA) web based personality assessment project.  Age, gender, and education are also reported. The data from 700 subjects are included here as a demonstration set for correlation and analysis.
}
\usage{data(sat.act)}
\format{
  A data frame with 700 observations on the following 6 variables.
  \describe{
    \item{\code{gender}}{males = 1, females = 2}
    \item{\code{education}}{self reported education 1 = high school ... 5 = graduate work}
    \item{\code{age}}{age}
      \item{\code{ACT}}{ACT composite scores may range from 1 - 36.  National norms have a mean of 20. }
    \item{\code{SATV}}{SAT Verbal scores may range from 200 - 800. }
    \item{\code{SATQ}}{SAT Quantitative scores may range from 200 - 800}
  }
}
\details{hese items were collected as part of the  SAPA project to develop online measures of ability.  The score means are higher than national norms suggesting both self selection for people taking on line personality and ability tests and a self reporting bias in scores.

See also the iq.items data set.  
}
\source{\url{http://personality-project.org} 
}
\references{
}
\examples{
data(sat.act)
describe(sat.act)
pairs.panels(sat.act)
}
\keyword{datasets}