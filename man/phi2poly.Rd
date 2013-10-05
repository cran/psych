\name{phi2tetra}
\alias{phi2tetra}
\alias{phi2poly}
\title{ Convert a phi coefficient to a tetrachoric correlation }
\description{Given a phi coefficient (a Pearson r calculated on two dichotomous variables), and the marginal frequencies (in percentages), what is the corresponding estimate of the tetrachoric correlation?

Given a two x two table of counts \cr
\tabular{lll}{
\tab a \tab b \cr
\tab c \tab d  \cr
}
The phi coefficient is (a - (a+b)*(a+c))/sqrt((a+b)(a+c)(b+d)(c+c)).

This function reproduces the cell entries for specified marginals and then calls John Fox's polychor function. 
}

\usage{
phi2tetra(ph,m,n=NULL,correct=TRUE)
phi2poly(ph,cp,cc,n=NULL,correct=TRUE) #deprecated
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ph}{phi }
  \item{m}{a vector of the selection ratio and probability of criterion.  In the case where ph is a matrix, m is a vector of the frequencies of the selected cases}
  \item{correct}{When finding tetrachoric correlations, should we correct for continuity for small marginals.  See \code{\link{tetrachoric}} for a discussion.}
  \item{n}{If the marginals are given as frequencies, what was the total number of cases?}
  \item{cp}{ probability of the predictor -- the so called selection ratio }
  \item{cc}{probability of the criterion -- the so called success rate. }
}
\details{requires the mvtnorm package
}
\value{a tetrachoric correlation
}
\author{ William Revelle}
\seealso{ \code{\link{tetrachoric}}, \code{\link{Yule2phi.matrix}}, \code{\link{phi2poly.matrix}} }
\examples{
 phi2tetra(.3,c(.5,.5))
#phi2poly(.3,.3,.7)
}
\keyword{ models }% at least one, from doc/KEYWORDS
\keyword{ models }% __ONLY ONE__ keyword per line
