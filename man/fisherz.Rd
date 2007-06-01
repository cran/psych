\name{fisherz}
\alias{fisherz}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Fisher z transform of r}
\description{convert a correlation to a z score using the Fisher transformation.
}
\usage{
fisherz(rho)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rho}{ a Pearson r }
}

\value{ z value corresponding to r
}

\author{
Maintainer: William Revelle  \email{revelle@northwestern.edu }
}

\examples{
\dontrun{
cors <- seq(-.9,.9,.1)
round(fisherz(cors),2)
}

}
\keyword{ multivariate }% at least one, from doc/KEYWORDS
\keyword{ models }% __ONLY ONE__ keyword per line
