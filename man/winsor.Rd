\name{winsor}
\alias{winsor}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Find the Winsorized mean for a vector, matrix, or data.frame }
\description{Among the  robust estimates of central tendency are trimmed means and Winsorized means.  This function finds the Winsorized mean.  The top and bottom trim values are given values of the trimmed and 1- trimmed quantiles.  Then means are found.    }
\usage{
winsor(x, trim = 0.2, na.rm = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A data vector, matrix or data frame}
  \item{trim}{Percentage of data to move  from the top and bottom of the distributions}
  \item{na.rm}{Missing data are removed }
}
\details{Among the many robust estimates of central tendency, some recommend the Winsorized mean.  Rather than just dropping the top and bottom trim percent, these extreme values are replaced with values at the trim and 1- trim quantiles.
}
\value{A scalar or vector winsorized means
}
\references{Wilcox, Rand R. (2005) Introduction to robust estimation and hypothesis testing. Elsevier/Academic Press. Amsterdam ; Boston.
  

 }
\author{William Revelle}
\seealso{ \code{\link{interp.median}} }
\examples{
data(sat.act)
winsor(sat.act)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{univar}

