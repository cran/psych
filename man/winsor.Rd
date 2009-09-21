\name{winsor}
\alias{winsor}
\alias{winsor.mean}
\alias{winsor.means}
\alias{winsor.sd}
\alias{winsor.var}

\title{Find the Winsorized scores, means, sds or variances  for a vector, matrix, or data.frame }
\description{Among the  robust estimates of central tendency are trimmed means and Winsorized means.  This function finds the Winsorized scores.   The top and bottom trim values are given values of the trimmed and 1- trimmed quantiles.  Then means, sds, and variances are found.    }
\usage{
winsor(x, trim = 0.2, na.rm = TRUE)
winsor.mean(x, trim = 0.2, na.rm = TRUE)
winsor.means(x, trim = 0.2, na.rm = TRUE)  
winsor.sd(x, trim = 0.2, na.rm = TRUE)  
winsor.var(x, trim = 0.2, na.rm = TRUE)  

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A data vector, matrix or data frame}
  \item{trim}{Percentage of data to move  from the top and bottom of the distributions}
  \item{na.rm}{Missing data are removed }
}
\details{Among the many robust estimates of central tendency, some recommend the Winsorized mean.  Rather than just dropping the top and bottom trim percent, these extreme values are replaced with values at the trim and 1- trim quantiles.
}
\value{A scalar or vector of winsorized scores or winsorized means, sds, or variances (depending upon the call).
}
\references{Wilcox, Rand R. (2005) Introduction to robust estimation and hypothesis testing. Elsevier/Academic Press. Amsterdam ; Boston.
  

 }
\author{William Revelle with modifications suggested by Joe Paxton and a further correction added  (January, 2009) to preserve the original order for the winsor case.}
\seealso{ \code{\link{interp.median}} }
\examples{
data(sat.act)
winsor.means(sat.act) #compare with the means of the winsorized scores
y <- winsor(sat.act)
describe(y)
xy <- data.frame(sat.act,y)
pairs.panels(xy) #to see the effect of winsorizing 
x <- matrix(1:100,ncol=5)
winsor(x)
winsor.means(x)
y <- 1:11
winsor(y,trim=.5)
}
\keyword{univar}

