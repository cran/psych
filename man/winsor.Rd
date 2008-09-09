\name{winsor}
\alias{winsor}
\alias{winsor.means}

\title{Find the Winsorized scores or means for a vector, matrix, or data.frame }
\description{Among the  robust estimates of central tendency are trimmed means and Winsorized means.  This function finds the Winsorized mean.  The top and bottom trim values are given values of the trimmed and 1- trimmed quantiles.  Then means are found.    }
\usage{
winsor(x, trim = 0.2, na.rm = TRUE)
winsor.means(x, trim = 0.2, na.rm = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A data vector, matrix or data frame}
  \item{trim}{Percentage of data to move  from the top and bottom of the distributions}
  \item{na.rm}{Missing data are removed }
}
\details{Among the many robust estimates of central tendency, some recommend the Winsorized mean.  Rather than just dropping the top and bottom trim percent, these extreme values are replaced with values at the trim and 1- trim quantiles.
}
\value{A scalar or vector of winsorized scores or winsorized means (depending upon the call).
}
\references{Wilcox, Rand R. (2005) Introduction to robust estimation and hypothesis testing. Elsevier/Academic Press. Amsterdam ; Boston.
  

 }
\author{William Revelle with modifications suggested by Joe Paxton}
\seealso{ \code{\link{interp.median}} }
\examples{
data(sat.act)
winsor.means(sat.act)
x <- matrix(1:100,ncol=5)
winsor(x)
winsor.means(x)
y <- 1:11
winsor(y,trim=.5)
}
\keyword{univar}

