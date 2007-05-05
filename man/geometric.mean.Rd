\name{geometric.mean}
\alias{geometric.mean}
%- Also NEED an '\alias' for EACH other topic documented here.
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
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.
x <- seq(1,5)
x2 <- x^2
geometric.mean(x)
geometric.mean(x2)
## The function is currently defined as
function(x) { exp(mean(log(x))) }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ multivariate }
