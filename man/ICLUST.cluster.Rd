\name{ICLUST.cluster}
\alias{ICLUST.cluster}

\title{Function to form hierarchical cluster analysis of items }
\description{
  The guts of the  \code{\link{ICLUST}} algorithm.  Called by \code{\link{ICLUST}} See ICLUST for description.
}
\usage{
ICLUST.cluster(r.mat, ICLUST.options)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{r.mat}{ A correlation matrix}
  \item{ICLUST.options}{ A list of options (see \code{\link{ICLUST}}) }
}
\details{
  See \code{\link{ICLUST}}
}
\value{
  A list of cluster statistics, described more fully in \code{\link{ICLUST}}
  \item{comp1 }{Description of 'comp1'}
  \item{comp2 }{Description of 'comp2'}
  ...
}
\references{Revelle, W. 1979, Hierarchical Cluster Analysis and the Internal Structure of Tests. Multivariate Behavioral Research, 14, 57-74. 
\url{http://personality-project.org/revelle/publications/iclust.pdf} \cr
See also  more extensive documentation at 
\url{http://personality-project.org/r/r.ICLUST.html} }
\author{William Revelle \cr
        Department of Psychology \cr
        Northwestern University \cr
        Evanston, Illinois \cr
        \email{  revelle@northwestern.edu           }\cr
 \url{http://personality-project.org/revelle.html}   }
\note{ Although the main code for ICLUST is here in ICLUST.cluster, the more extensive documentation is for \code{\link{ICLUST}}.
}
\seealso{ \code{\link{ICLUST.graph}},\code{\link{ICLUST}}, \code{\link{cluster.fit} }, \code{\link{VSS}}, \code{\link{omega} }}

\keyword{multivariate }
\keyword{cluster }
