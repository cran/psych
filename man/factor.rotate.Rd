\name{factor.rotate}
\alias{factor.rotate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{``Hand" rotate a factor loading matrix }
\description{Given a factor or components matrix, it is sometimes useful to do arbitrary rotations of particular pairs of variables.  This supplements the much more powerful rotation package GRProtation and is meant for specific requirements to do unusual rotations.
}
\usage{
factor.rotate(f, angle, col1, col2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{f}{original loading matrix }
  \item{angle}{ angle (in degrees!) to rotate }
  \item{col1}{ column in factor matrix defining the first variable}
  \item{col2}{ column in factor matrix defining the second variable }
}
\details{Partly meant as a demonstration of how rotation works, factor.rotate is useful for those cases that require specific rotations that are not available in more advanced packages such as GPArotation.
}
\value{the resulting rotated matrix of loadings.
}
\references{ \url{http://personality-project.org/revelle/syllabi/405.syllabus.html} \cr }
\author{
Maintainer: William Revelle \email{revelle@northwestern.edu  }
}
\note{For a complete rotation package, see GPArotation }
\examples{


}
\keyword{ multivariate }% at least one, from doc/KEYWORDS
\keyword{ models }% __ONLY ONE__ keyword per line
