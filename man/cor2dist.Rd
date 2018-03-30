\name{cor2dist}
\alias{cor2dist}
\title{Convert correlations to distances (necessary to do multidimensional scaling of correlation data)}

\description{A minor helper function to convert correlations (ranging from -1 to 1) to distances (ranging from 0 to 2). \eqn{d = \sqrt{(2(1-r))}}. 
}
\usage{
cor2dist(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{If square, then assumed to be a correlation matrix, otherwise the correlations are found first.
}
}

\value{dist: a square matrix of distances.

}
\references{
Revelle, William. (in prep) An introduction to psychometric theory with applications in R. Springer.  Working draft available at \url{https://personality-project.org/r/book/} 

}
\author{William Revelle}

\note{
For an example of doing multidimensional scaling on data that are normally factored, see Revelle (in prep)
}


\keyword{ multivariate }
\keyword{ models}

