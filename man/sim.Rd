\name{sim}
\alias{sim}
\title{Functions to simulate psychological/psychometric data}
\description{A number of functions in the psych package will generate simulated data.  These functions include
\code{\link{sim.circ}}, \code{\link{sim.congeneric}}, \code{\link{sim.dichot}}, \code{\link{sim.hierarchical}}, \code{\link{sim.item}}, \code{\link{sim.structural}},  and \code{\link{sim.VSS}}. These functions are separately documented and are listed here for ease of the help function.  See each function for more detailed help.
}
\usage{
sim()
}

\details{Simulation of data structures is a very useful tool in psychometric research and teaching.  By knowing ``truth" it is possible to see how well various algorithms can capture it.  Various simulation functions in psych are:

\code{\link{sim.structural}}  A function to combine a measurement and structural model into one data matrix.  Useful for understanding structural equation models.


\code{\link{sim.congeneric}}   A function to create congeneric items/tests for demonstrating classical test theory.
 
\code{\link{sim.hierarchical}}  A function to create data with a hierarchical (bifactor) structure.  

\code{\link{sim.item}}      A function to create items that either have a simple structure or a circumplex structure.

\code{\link{sim.circ}}    Create data with a circumplex structure.

\code{\link{sim.dichot}}    Create dichotomous item data with a simple or circumplex structure.



}
\references{Revelle, W. (in preparation) An Introduction to Psychometric Theory with applications in R. Springer. at \url{http://personality-project.org/r/book/}  }

\author{William Revelle}

\seealso{ See above}
\examples{
sim.congeneric()
R <- sim.hierarchical()
R
fx <- matrix(c(.9,.8,.7,rep(0,6),c(.8,.7,.6)),ncol=2)
fy <- c(.6,.5,.4)
Phi <- matrix(c(1,0,.5,0,1,.4,0,0,0),ncol=3)
print(sim.structural(fx,fy,Phi),digits=2)  
}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{multivariate}
\keyword{datagen}

