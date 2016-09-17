\name{count.pairwise}
\alias{count.pairwise}
\alias{pairwiseDescribe}
\title{ Count number of pairwise cases for a data set with missing (NA) data. }
\description{
  When doing cor(x, use= "pairwise"), it is nice to know the number of cases for each pairwise correlation.  This is particularly useful when doing SAPA type analyses.
}
\usage{
count.pairwise(x, y = NULL,diagonal=TRUE)
pairwiseDescribe(x,diagonal=FALSE,...) 

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ An input matrix, typically a data matrix ready to be correlated. }
  \item{y}{ An optional second input matrix }
  \item{diagonal}{if TRUE, then report the diagonal, else fill the diagonals with NA}
  \item{...}{Other parameters to pass to describe}
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
y[y > 1] <- NA

count.pairwise(x)
count.pairwise(y)
count.pairwise(x,y)
count.pairwise(x,diagonal=FALSE)
pairwiseDescribe(x,quant=c(.1,.25,.5,.75,.9))
}

    
}
\keyword{ models }% at least one, from doc/KEYWORDS
\keyword{ multivariate }% __ONLY ONE__ keyword per line
