\name{RMSEA}
\alias{RMSEA}
\title{Root Mean Squared Error of Approximation from chisq, df, and n}
\description{
Find the RMSEA from model chi square, degrees of freedom and number of observations.
Show confidence intervals. 
}
\usage{
RMSEA(chisq, dof, n.obs, alpha = 0.1)
}
\arguments{
  \item{chisq}{The Chi square statistic from an analysis}
  \item{dof}{Degrees of freedom of the model}
  \item{n.obs}{Number of observations}
  \item{alpha}{alpha level for confidence intervals}
}
\details{RMSEA is just a chisquare adjusted by its degrees of freedom and the sample size.
\eqn{\sqrt{chisq/(dof * (nobs)) -1/(n.obs -1)}}.  It is given in most of the appropriate functions (e.g.,\code{\link{fa}}) and given here for completeness.}

\value{
  \item{RMSEA}{the estimated value}
  \item{RMSEA-L}{the lower bound}
  \item{RMSEA-U}{the upper bound}
  }
\references{ Steiger, J. H.,and  Lind, J. C. (1980). Statistically based tests for the number of common factors. Paper presented at the Annual Meeting of the Psychometric Society, Iowa City, IA.
}
\author{
William Revelle
}

\seealso{\code{\link{fa}}, \code{\link{omega}} , \code{\link{esem}}

}


\keyword{ multivariate }
\keyword{ models}