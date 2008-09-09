\name{factor.model}
\alias{factor.model}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Find R =  F F' + U2 is the basic factor model }
\description{The basic factor or principal components model is that a correlation or covariance matrix may be reproduced by the product of a factor loading matrix times its transpose.  Find this reproduced matrix.  Used by \code{\link{factor.fit}}, \code{\link{VSS}}, \code{\link{ICLUST}}, etc.
}
\usage{
factor.model(f)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{f}{ A matrix of loadings. }
}
\details{
}
\value{ A correlation or covariance matrix.
}
\references{Gorsuch, Richard, (1983) Factor Analysis. Lawrence Erlebaum Associates. 
\cr
Revelle, W. In preparation) An Introduction to Psychometric Theory with applications in R (\url{http://personality-project.org/r/book/}) } 
\author{ 
 \email{revelle@northwestern.edu } \cr
 \url{http://personality-project.org/revelle.html} \cr
}
\seealso{ \code{\link{ICLUST.graph}},\code{\link{ICLUST.cluster}}, \code{\link{cluster.fit} }, \code{\link{VSS}}, \code{\link{omega} }}
\examples{

f2 <- matrix(c(.9,.8,.7,rep(0,6),.6,.7,.8),ncol=2)
mod <- factor.model(f2)
round(mod,2)
}
\keyword{ multivariate }% at least one, from doc/KEYWORDS
\keyword{ models }% __ONLY ONE__ keyword per line
