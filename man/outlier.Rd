\name{outlier}
\alias{outlier}
\title{Find and graph Mahalanobis squared distances to detect outliers}
\description{
The Mahalanobis distance is \eqn{D^2 = (x-\mu)' \Sigma^-1 (x-\mu)} where \eqn{\Sigma} is the covariance of the x matrix.  D2 may be used as a way of detecting outliers in distribution.  Large D2 values, compared to the expected Chi Square values indicate an unusual response pattern.  The mahalanobis function in stats does not handle missing data.
}
\usage{
outlier(x, plot = TRUE, bad = 5,na.rm = TRUE, xlab, ylab, ...)
}

\arguments{
  \item{x}{A data matrix or data.frame}

 \item{plot}{Plot the resulting QQ graph}
  \item{bad}{Label the bad worst values}

 \item{na.rm}{Should missing data be deleted}
  \item{xlab}{Label for x axis}
  \item{ylab}{Label for y axis}
  \item{\dots}{More graphic parameters, e.g., cex=.8}
}
\details{
Adapted from the mahalanobis function and help page from stats.
}
\value{The D2 values for each case}
\references{
 Yuan, Ke-Hai and Zhong, Xiaoling, (2008) Outliers, Leverage Observations, and Influential Cases in Factor Analysis: Using Robust Procedures to Minimize Their Effect, Sociological Methodology, 38, 329-368.
 }
\author{William Revelle}

\seealso{
 \code{\link{mahalanobis}} 
}
\examples{
#first, just find and graph the outliers
d2 <- outlier(sat.act)
#combine with the data frame and plot it with the outliers highlighted in blue
sat.d2 <- data.frame(sat.act,d2)
pairs.panels(sat.d2,bg=c("yellow","blue")[(d2 > 25)+1],pch=21)
}

\keyword{multivariate }
\keyword{ models }
