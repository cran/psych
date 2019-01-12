\name{mat.sort}
\Rdversion{1.1}
\alias{mat.sort}
\alias{matSort}
\title{Sort the elements of a correlation matrix to reflect factor loadings}
\description{To see the structure of a correlation matrix, it is helpful to organize the items so that the similar items are grouped together. One such grouping technique is factor analysis.  mat.sort will sort the items by a factor model (if specified), or any other order, or by the loadings on the first factor (if unspecified)
}
\usage{
mat.sort(m, f = NULL)
matSort(m, f = NULL)
}

\arguments{
  \item{m}{A correlation matrix
}
  \item{f}{A factor analysis output (i.e., one with a loadings matrix) or a matrix of weights
}
}
\details{The factor analysis output is sorted by size of the largest  factor loading for each variable  and then the matrix items are organized by those loadings.  The default is to sort by the loadings on the first factor.  Alternatives allow for ordering based upon any vector or matrix. 
}
\value{
A sorted correlation matrix, suitable for showing with \code{\link{corPlot}}.
}

\author{William Revelle}




\seealso{ \code{\link{fa}}, \code{\link{corPlot}}
}
\examples{data(Bechtoldt.1)
sorted <- mat.sort(Bechtoldt.1,fa(Bechtoldt.1,5))
corPlot(sorted,xlas=2) #vertical xaxis names
}

\keyword{ multivariate}% at least one, from doc/KEYWORDS
\keyword{models }% __ONLY ONE__ keyword per line
