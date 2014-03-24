\name{fa.sort}
\Rdversion{1.1}
\alias{fa.sort}
\alias{fa.organize}


\title{Sort factor analysis or principal components analysis loadings}

\description{Although the print.psych function will sort factor analysis loadings, sometimes it is useful to do this outside of the print function. fa.sort takes the output from the fa or principal functions and sorts the loadings for each factor.  Items are located in terms of their greatest loading.  The new order is returned as an element in the fa list.
}
\usage{
fa.sort(fa.results,polar=FALSE)
fa.organize(fa.results,o=NULL,i=NULL,cn=NULL) 
}


\arguments{
  \item{fa.results}{The output from a factor analysis or principal components analysis using \code{\link{fa}} or \code{\link{principal}}. }
  \item{polar}{Sort by polar coordinates of first two factors (FALSE)}
\item{o}{The order in which to order the factors}
\item{i}{The order in which to order the items}
\item{cn}{new factor names}
}
\details{
The fa.results$loadings are replaced with sorted loadings.

fa.organize takes a factor analysis or components output and reorganizes the factors in the o order. Items are organized in the i order.  This is useful when comparing alternative factor solutions.
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
test.simple <- fa(sim.item(16),2)
fa.sort(test.simple)
fa.organize(test.simple,c(2,1))  #the factors but not the items have been rearranged
}

\keyword{ multivariate }

