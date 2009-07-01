\name{VSS.scree}
\alias{VSS.scree}
\title{  Plot a scree test}
\description{Cattell's scree test is one of most simple ways of testing the number of components in  a correlation matrix. Here we plot the  eigen values of a correlation matrix.
}
\usage{
VSS.scree(rx, main = "scree plot")
}
\details{Among the many ways to choose the optimal number of factors is the scree test.  A better function to show the scree as well as compare it to randomly parallel solutions is found found in \code{\link{fa.parallel}}
}
\arguments{
  \item{rx}{ a correlation matrix or a data matrix. If data, then correlations are found using pairwise deletions. }
  \item{main}{ Title }
}

\references{ \url{http://personality-project.org/r/vss.html}}
\author{William Revelle 
}
\seealso{  \code{\link{fa.parallel}} \code{\link{VSS.plot}},  \code{\link{ICLUST}}, \code{\link{omega}}}
\examples{
#VSS.scree(attitude)
#VSS.scree(cor(attitude)


}
\keyword{ multivariate }% at least one, from doc/KEYWORDS
\keyword{ hplot }% __ONLY ONE__ keyword per line
