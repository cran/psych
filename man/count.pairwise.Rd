\name{count.pairwise}
\alias{count.pairwise}
\title{ Count number of pairwise cases for a data set with missing (NA) data. }
\description{
  When doing cor(x, use= "pairwise"), it is nice to know the number of cases for each pairwise correlation.  This is particularly useful when doing SAPA type analyses.
}
\usage{
count.pairwise(x, y = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ An input matrix, typically a data matrix ready to be correlated. }
  \item{y}{ An optional second input matrix }
}

\value{result = matrix of counts of pairwise observations
}

\author{

Maintainer: William Revelle \email{revelle@northwestern.edu}}

\examples{
\dontrun{
x <- matrix(rnorm(1000),ncol=6)
y <- matrix(rnorm(500),ncol=3)
x[x < 0] <- NA
y[y> 1] <- NA

count.pairwise(x)
count.pairwise(y)
count.pairwise(x,y)
}

    
}
\keyword{ models }% at least one, from doc/KEYWORDS
\keyword{ multivariate }% __ONLY ONE__ keyword per line
