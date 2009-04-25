\name{sim}
\alias{sim}
\alias{sim.simplex}
\title{Functions to simulate psychological/psychometric data.}
\description{A number of functions in the psych package will generate simulated data.  These functions include
\code{link{sim}} (for a factor simplex), and \code{link{sim}}( for a data simplex) \code{\link{sim.circ}} for a circumplex structure, \code{\link{sim.congeneric}} (for a one factor factor congeneric model), \code{\link{sim.dichot}} (to simulate dichotomous items), \code{\link{sim.hierarchical}} (a hierarchical factor model), \code{\link{sim.item}} (general item simulations), \code{\link{sim.structural}} (general simulation of structural models), \code{\link{sim.anova}} for ANOVA and lm simulatins, and \code{\link{sim.VSS}}. These functions are separately documented and are listed here for ease of the help function.  See each function for more detailed help.
}
\usage{
sim(fx=NULL,Phi=NULL,fy=NULL,n=0,mu=NULL,raw=FALSE)
sim.simplex(nvar =12, r=.8,mu=NULL, n=0)
}
\arguments{
  \item{fx}{The measurement model for x. If NULL, a 4 factor model is generated}
  \item{Phi}{The structure matrix of the latent variables}
  \item{fy}{The measurement model for y}
  \item{mu}{The means structure for the fx factors}
  \item{n}{ Number of cases to simulate.  If n=0 or NULL, the population matrix is returned.}
  \item{raw}{if raw=TRUE, raw data are returned as well.}\
  \item{nvar}{Number of variables for a simplex structure}
  \item{r}{the base correlation for a simplex}
 }

\details{Simulation of data structures is a very useful tool in psychometric research and teaching.  By knowing ``truth" it is possible to see how well various algorithms can capture it.  

The default values for \code{\link{sim.structure}} is to generate a 4 factor, 12 variable data set with a simplex structure between the factors.  

Other simulation functions in psych are:

\code{\link{sim.structure}}  A function to combine a measurement and structural model into one data matrix.  Useful for understanding structural equation models.


\code{\link{sim.congeneric}}   A function to create congeneric items/tests for demonstrating classical test theory. This is just a special case of sim.structure.
 
\code{\link{sim.hierarchical}}  A function to create data with a hierarchical (bifactor) structure.  

\code{\link{sim.item}}      A function to create items that either have a simple structure or a circumplex structure.

\code{\link{sim.circ}}    Create data with a circumplex structure.

\code{\link{sim.dichot}}    Create dichotomous item data with a simple or circumplex structure.

\code{\link{sim.anova}}    Simulate a 3 way balanced ANOVA or linear model, with or without repeated measures. 



}
\references{Revelle, W. (in preparation) An Introduction to Psychometric Theory with applications in R. Springer. at \url{http://personality-project.org/r/book/}  }

\author{William Revelle}

\seealso{ See above}
\examples{
simplex <- sim()
round(simplex$model,2)

congeneric <- sim.congeneric()
round(congeneric,2)
R <- sim.hierarchical()
R
fx <- matrix(c(.9,.8,.7,rep(0,6),c(.8,.7,.6)),ncol=2)
fy <- c(.6,.5,.4)
Phi <- matrix(c(1,0,.5,0,1,.4,0,0,0),ncol=3)
print(sim.structure(fx,Phi,fy,),digits=2)  
cor.plot(R) #show it graphically

simp <- sim.simplex()
#show the simplex structure using cor.plot
cor.plot(simp)
}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{multivariate}
\keyword{datagen}

