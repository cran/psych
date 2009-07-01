\name{factor.rotate}
\alias{factor.rotate}
\title{``Hand" rotate a factor loading matrix }
\description{Given a factor or components matrix, it is sometimes useful to do arbitrary rotations of particular pairs of variables.  This supplements the much more powerful rotation package GPArotation and is meant for specific requirements to do unusual rotations.
}
\usage{
factor.rotate(f, angle, col1=1, col2=2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{f}{original loading matrix or a data frame  (can be output from a factor analysis function}
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
#using the Harman 24 mental tests, rotate the 2nd and 3rd factors 45 degrees
pc <- principal(Harman74.cor$cov,4,rotate=TRUE)
pcr45 <- factor.rotate(pc,45,2,3)
pcr90 <- factor.rotate(pcr45,45,2,3)
print(factor.congruence(pc,pcr45),digits=3) #poor congruence with original
print(factor.congruence(pc,pcr90),digits=3) #factor 2 and 3 have been exchanged and 3 flipped
}
\keyword{ multivariate }% at least one, from doc/KEYWORDS
\keyword{ models }% __ONLY ONE__ keyword per line
