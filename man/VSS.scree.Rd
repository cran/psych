\name{VSS.scree}
\alias{VSS.scree}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{  Plot a scree test}
\description{Cattell's scree test is one of most simple ways of testing the number of components in  a correlation matrix. Here we plot the  eigen values of a correlation matrix.
}
\usage{
VSS.scree(rx, main = "scree plot")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rx}{ a correlation matrix or a data matrix. If data, then correlations are found using pairwise deletions. }
  \item{main}{ Title }
}

\references{ \url{http://personality-project.org/r/vss.html}}
\author{William Revelle \cr
        Department of Psychology \cr
        Northwestern University \cr
        Evanston, Illiniois \cr


Maintainer: William Revelle <revelle@northwestern.edu>

}
\seealso{  \code{\link{VSS.plot}},  \code{\link{ICLUST}}, \code{\link{omega}}}
\examples{
#VSS.scree(attitude)
#VSS.scree(cor(attitude)


}
\keyword{ multivariate }% at least one, from doc/KEYWORDS
\keyword{ hplot }% __ONLY ONE__ keyword per line
