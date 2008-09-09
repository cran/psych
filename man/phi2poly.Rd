\name{phi2poly}
\alias{phi2poly}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Convert a phi coefficient to a polychoric correlation }
\description{Given a phi coefficient (a Pearson r calculated on two dichotomous variables), and the marginal frequencies (in percentages), what is the corresponding estimate of the polychoric correlation?

Given a two x two table of counts \cr
\tabular{lll}{
\tab a \tab b \cr
\tab c \tab d  \cr
}
The phi coefficient is (a - (a+b)*(a+c))/sqrt((a+b)(a+c)(b+d)(c+c)).

This function reproduces the cell entries for specified marginals and then calls John Fox's polychor function. 
}

\usage{
phi2poly(ph, cp, cc)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ph}{phi }
  \item{cp}{ probability of the predictor -- the so called selection ratio }
  \item{cc}{probability of the criterion -- the so called success rate. }
}
\details{Uses John Fox's polycor package, which in turn requires the mvtnorm package
}
\value{a polychoric correlation
}
\author{ William Revelle}
\seealso{  \code{\link{polychor.matrix}}, \code{\link{Yule2phi.matrix}}, \code{\link{phi2poly.matrix}} }
\examples{
#phi2poly(.3,.5,.5)
#phi2poly(.3,.3,.7)
}
\keyword{ models }% at least one, from doc/KEYWORDS
\keyword{ models }% __ONLY ONE__ keyword per line
