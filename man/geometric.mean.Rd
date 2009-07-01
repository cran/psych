\name{geometric.mean}
\alias{geometric.mean}
\title{ Find the geometric mean of a vector or columns of a data.frame. }
\description{The geometric mean is the nth root of n products or e to the mean log of x.
Useful for describing non-normal, i.e., geometric distributions.
}
\usage{
geometric.mean(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ a vector or data.frame}
}
\details{Useful for teaching how to write functions, also useful for showing the different ways of estimating central tendency. 
}
\value{geometric mean(s) of x or x.df.
}

\author{ William Revelle}
\note{ Not particularly useful if there are elements that are <= 0.
}
\seealso{ \code{\link{harmonic.mean}}, \code{\link{mean}} }
\examples{

x <- seq(1,5)
x2 <- x^2
geometric.mean(x)
geometric.mean(x2)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ multivariate }
