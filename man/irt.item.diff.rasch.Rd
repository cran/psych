\name{irt.item.diff.rasch}
\alias{irt.item.diff.rasch}
\alias{irt.discrim}
\title{Simple function to estimate item difficulties using IRT concepts}
\description{Steps toward a  very crude and preliminary IRT program. These two functions  estimate item difficulty and discrimination parameters.  A better procedure is to use \code{\link{irt.fa}} or the ltm package.
}
\usage{
irt.item.diff.rasch(items)
irt.discrim(item.diff,theta,items)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{items}{ a matrix of items }
  \item{item.diff}{a vector of item difficulties (found by irt.item.diff)}
  \item{theta}{ability estimate from irt.person.theta}
}
\details{Item Response Theory (aka "The new psychometrics") models individual responses to items with a logistic function and an individual (theta) and item difficulty (diff) parameter.

irt.item.diff.rasch finds item difficulties with the assumption of theta=0 for all subjects and that all items are equally discriminating.

irt.discrim takes those  difficulties and theta estimates from \code{\link{irt.person.rasch}} to find item discrimination (beta) parameters.

A far better package with these features is the ltm package.  The IRT functions in the psych-package are for pedagogical rather than production purposes.  They are believed to be accurate, but are not guaranteed. They do seem to be slightly more robust to missing data structures associated with SAPA data sets than the ltm package. 

The \code{\link{irt.fa}} function is also an alternative. This will find \code{\link{tetrachoric}} or \code{\link{polychoric}} correlations and then convert to IRT parameters using factor analysis (\code{\link{fa}}).
}
\value{a vector of item difficulties or item discriminations.

}

\author{William Revelle }
\note{ Under development. Not recommended for public consumption.  See \code{\link{irt.fa}} and \code{\link{score.irt}}  for far better options. 
}
\seealso{ \code{\link{irt.fa}}, \code{\link{irt.person.rasch}} }

\keyword{ multivariate}% at least one, from doc/KEYWORDS
\keyword{ models}% __ONLY ONE__ keyword per line
