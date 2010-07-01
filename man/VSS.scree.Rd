\name{VSS.scree}
\alias{VSS.scree}
\alias{scree}
\title{Plot the successive eigen values for a scree test}
\description{Cattell's scree test is one of most simple ways of testing the number of components or factors in  a correlation matrix. Here we plot the  eigen values of a correlation matrix as well as the eigen values of a factor analysis.
}
\usage{
scree(rx,factors=TRUE,main="Scree plot",add=FALSE) 
VSS.scree(rx, main = "scree plot")

}
\details{Among the many ways to choose the optimal number of factors is the scree test.  A better function to show the scree as well as compare it to randomly parallel solutions is found found in \code{\link{fa.parallel}}
}
\arguments{
  \item{rx}{ a correlation matrix or a data matrix. If data, then correlations are found using pairwise deletions. }
 \item{factors}{If true, draw the scree for factors and components}
 \item{main}{ Title }
 \item{add}{Should multiple plots be drawn?}
}

\references{ \url{http://personality-project.org/r/vss.html}}
\author{William Revelle 
}
\seealso{  \code{\link{fa.parallel}} \code{\link{VSS.plot}},  \code{\link{ICLUST}}, \code{\link{omega}}}
\examples{
scree(attitude)
#VSS.scree(cor(attitude)


}
\keyword{ multivariate }% at least one, from doc/KEYWORDS
\keyword{ hplot }% __ONLY ONE__ keyword per line
