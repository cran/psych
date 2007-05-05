\name{factor.rotate}
\alias{factor.rotate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ "Hand" rotate a factor loading matrix }
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
\details{Partly meant as a demonstration of how rotation works, factor.rotate is useful for those cases that require specific rotations that are not available in more advanced packages such as GPAroation.
}
\value{the resulting rotated matrix of loadings.
}
\references{ \url{http://personality-project.org/revelle/syllabi/405.syllabus.html} \cr }
\author{
Maintainer: William Revelle \email{revelle@northwestern.edu  }
}
\note{For a complete rotation package, see GPArotation }
\examples{


## The function is currently defined as
function(f,angle,col1,col2)  {
   #hand rotate two factors from a loading matrix
   #see the GPArotation package for much more elegant procedures
     nvar<- dim(f)[2]
     rot<- matrix(rep(0,nvar*nvar),ncol=nvar)
     rot[cbind(1:nvar, 1:nvar)] <- 1
     theta<- 2*pi*angle/360
     rot[col1,col1]<-cos(theta)
     rot[col2,col2]<-cos(theta)
     rot[col1,col2]<- -sin(theta)
     rot[col2,col1]<- sin(theta)
     result <- f \%*\% rot
     return(result) }
}
\keyword{ multivariate }% at least one, from doc/KEYWORDS
\keyword{ models }% __ONLY ONE__ keyword per line
