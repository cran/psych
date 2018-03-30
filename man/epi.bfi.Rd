\name{epi.bfi}
\alias{epi.bfi}
\docType{data}
\title{13 personality scales from the Eysenck Personality Inventory and Big 5 inventory}
\description{A small data set of 5 scales from the Eysenck Personality Inventory, 5 from a Big 5 inventory, a Beck Depression Inventory, and State and Trait Anxiety measures.  Used for demonstrations of correlations, regressions, graphic displays.
}
\usage{data(epi.bfi)}
\format{
  A data frame with 231 observations on the following 13 variables.
  \describe{
    \item{\code{epiE}}{EPI Extraversion }
    \item{\code{epiS}}{EPI Sociability (a subset of Extraversion items}
    \item{\code{epiImp}}{EPI Impulsivity (a subset of Extraversion items}
    \item{\code{epilie}}{EPI Lie scale}
    \item{\code{epiNeur}}{EPI neuroticism}
    \item{\code{bfagree}}{Big 5 inventory (from the IPIP) measure of Agreeableness}
    \item{\code{bfcon}}{Big 5 Conscientiousness}
    \item{\code{bfext}}{Big 5 Extraversion}
    \item{\code{bfneur}}{Big 5 Neuroticism}
    \item{\code{bfopen}}{Big 5 Openness}
    \item{\code{bdi}}{Beck Depression scale}
    \item{\code{traitanx}}{Trait Anxiety}
    \item{\code{stateanx}}{State Anxiety}
  }
}
\details{Self report personality scales tend to measure the ``Giant 2" of Extraversion and Neuroticism or the ``Big 5" of Extraversion, Neuroticism, Agreeableness, Conscientiousness, and Openness.  Here is a small data set from Northwestern University undergraduates with scores on the Eysenck Personality Inventory (EPI) and a Big 5 inventory taken from the International Personality Item Pool.  
}
\source{Data were collected at the Personality, Motivation, and Cognition Lab (PMCLab) at Northwestern by William Revelle)
}
\references{\url{https://personality-project.org/pmc.html}
}
\examples{
data(epi.bfi)
pairs.panels(epi.bfi[,1:5])
describe(epi.bfi)
}
\keyword{datasets}
