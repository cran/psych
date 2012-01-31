\name{fa.sort}
\Rdversion{1.1}
\alias{fa.sort}

\title{Sort factor analysis or principal components analysis loadings}

\description{Although the print.psych function will sort factor analysis loadings, sometimes it is useful to do this outside of the print function. fa.sort takes the output from the fa or principal functions and sorts the loadings for each factor.  Items are located in terms of their greatest loading.  The new order is returned as an element in the fa list.
}
\usage{
fa.sort(fa.results,polar=FALSE)
}

\arguments{
  \item{fa.results}{The output from a factor analysis or principal components analysis using \code{\link{fa}} or \code{\link{principal}}. }
  \item{polar}{Sort by polar coordinates of first two factors (FALSE)}
}
\details{
The fa.results$loadings are replaced with sorted loadings.
}
\value{ A sorted factor analysis, principal components analysis, or omega loadings matrix. 

These sorted values are used internally by the various diagram functions. 

The values returned are the same as  \code{\link{fa}}, except in sorted order.  In addition, the order is returned as an additional element in the fa list.

 

}

\author{William Revelle
}

\seealso{
 See Also as \code{\link{fa}},\code{\link{print.psych}}, \code{\link{fa.diagram}},
}
\examples{
test.simple <- fa(item.sim(16),2,rotate="oblimin")
fa.sort(test.simple)
}

\keyword{ multivariate }

