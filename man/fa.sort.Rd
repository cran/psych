\name{fa.sort}
\Rdversion{1.1}
\alias{fa.sort}

\title{Sort factor analysis or principal components analysis loadings}

\description{Although the print.psych function will sort factor analysis loadings, sometimes it is useful to do this outside of the print function. fa.sort takes the output from the fa or principal functions and sorts the loadings for each factor.  Items are located in terms of their greatest loading.
}
\usage{
fa.sort(fa.results)
}

\arguments{
  \item{fa.results}{The output from a factor analysis or principal components analysis using \code{\link{fa}} or \code{\link{principal}}. 
}
}
\details{
The fa.results$loadings are replaced with sorted loadings.
}
\value{ A sorted factor analysis, principal components analysis, or omega loadings matrix. 

These sorted values are used internally by the various diagram functions.  

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

