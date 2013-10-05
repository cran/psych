\name{polychor.matrix}
\alias{polychor.matrix}
\alias{Yule2poly.matrix}
\alias{phi2poly.matrix}
\alias{Yule2phi.matrix}
\title{Phi or Yule coefficient matrix  to polychoric coefficient matrix}
\description{A set of deprecated functions that have replaced by \code{\link{Yule2tetra}} and \code{\link{Yule2phi}}. 

Some older correlation matrices were reported as matrices of Phi or of Yule correlations.  That is, correlations were found from the two by two table of counts:
 \cr
\tabular{lll}{
\tab a \tab b \cr
\tab c \tab d  \cr
}
Yule Q is (ad - bc)/(ad+bc). \cr

With marginal frequencies of a+b, c+d, a+c, b+d.

Given a square matrix of such correlations, and the proportions for each variable that are in the a + b cells, it is possible to reconvert each correlation into a two by two table and then estimate the corresponding polychoric correlation (using John Fox's polychor function. 
}
\usage{
Yule2poly.matrix(x, v)  #deprectated
phi2poly.matrix(x, v)     #deprectated
Yule2phi.matrix(x, v)     #deprectated
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a matrix of phi or Yule  coefficients }
  \item{v}{A vector of marginal frequencies  }
}
\details{These functions call \code{\link{Yule2poly}},  \code{\link{Yule2phi}} or \code{\link{phi2poly}} for each cell of the matrix. See those functions for more details.  See \code{\link{phi.demo}} for an example.

}
\value{A matrix of correlations
}

\author{ William Revelle}

\examples{
#demo <- phi.demo() 
#compare the phi (lower off diagonal and polychoric correlations (upper off diagonal)
#show the result from poly.mat
#round(demo$tetrachoric$rho,2)
#show the result from phi2poly
#tetrachorics above the diagonal, phi below the diagonal 
#round(demo$phis,2) 


}
\keyword{ models }% at least one, from doc/KEYWORDS
\keyword{ multivariate }% __ONLY ONE__ keyword per line
