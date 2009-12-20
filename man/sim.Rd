\name{sim}
\alias{sim}
\alias{sim.simplex}
\alias{sim.minor}
\title{Functions to simulate psychological/psychometric data.}
\description{A number of functions in the psych package will generate simulated data.  These functions include
\code{link{sim}} (for a factor simplex), and \code{\link{sim.simplex}}( for a data simplex) \code{\link{sim.circ}} for a circumplex structure, \code{\link{sim.congeneric}} (for a one factor factor congeneric model), \code{\link{sim.dichot}} (to simulate dichotomous items), \code{\link{sim.hierarchical}} (a hierarchical factor model), \code{\link{sim.item}} (general item simulations),
\code{\link{sim.minor}} (to simulate major and minor factors)
\code{\link{sim.structural}} (general simulation of structural models), \code{\link{sim.anova}} for ANOVA and lm simulatins, and \code{\link{sim.VSS}}. These functions are separately documented and are listed here for ease of the help function.  See each function for more detailed help.
}
\usage{
sim(fx=NULL,Phi=NULL,fy=NULL,n=0,mu=NULL,raw=TRUE)
sim.simplex(nvar =12, r=.8,mu=NULL, n=0)
sim.minor(nvar=12,nfact=3,n=0,fbig=NULL,fsmall = NULL)
}
\arguments{
  \item{fx}{The measurement model for x. If NULL, a 4 factor model is generated}
  \item{Phi}{The structure matrix of the latent variables}
  \item{fy}{The measurement model for y}
  \item{mu}{The means structure for the fx factors}
  \item{n}{ Number of cases to simulate.  If n=0 or NULL, the population matrix is returned.}
  \item{raw}{if raw=TRUE, raw data are returned as well.}
  \item{nvar}{Number of variables for a simplex structure}
  \item{nfact}{Number of large factors to simulate in sim.minor}
  \item{r}{the base correlation for a simplex}
  \item{fbig}{Factor loadings for the main factors.  Default is a simple structure with loadings sampled from (.8,.6) for nvar/nfact variables and 0 for the remaining.  If fbig is specified, then it will be used.}
  \item{fsmall}{nvar/2 small factors are generated with loadings sampled from (-.2,0,.2)}
 }

\details{Simulation of data structures is a very useful tool in psychometric research and teaching.  By knowing ``truth" it is possible to see how well various algorithms can capture it.  

The default values for \code{\link{sim.structure}} is to generate a 4 factor, 12 variable data set with a simplex structure between the factors.  

Two data structures that are particular challenges to exploratory factor analysis are the simplex structure and the presence of minor factors.  Simplex structures \code{\link{sim.simplex}}will typically occur in developmental or learning contexts and have a correlation structure of r between adjacent variables and r^n for variables n apart.  Although just one latent variable (r) needs to be estimated, the structure will have nvar-1 factors.  

Many simulations of factor structures assume that except for the major factors, all residuals are normally distributed around 0.  An alternative, and perhaps more realistic situation, is that the there are a few major (big) factors and many minor (small) factors.  The challenge is thus to identify the major factors. \code{\link{sim.minor}} generates such structures.


Other simulation functions in psych are:

\code{\link{sim.structure}}  A function to combine a measurement and structural model into one data matrix.  Useful for understanding structural equation models.  Combined with \code{\link{structure.diagram}} to see the proposed structure.


\code{\link{sim.congeneric}}   A function to create congeneric items/tests for demonstrating classical test theory. This is just a special case of sim.structure.
 
\code{\link{sim.hierarchical}}  A function to create data with a hierarchical (bifactor) structure.  

\code{\link{sim.item}}      A function to create items that either have a simple structure or a circumplex structure.

\code{\link{sim.circ}}    Create data with a circumplex structure.

\code{\link{sim.dichot}}    Create dichotomous item data with a simple or circumplex structure.

\code{\link{sim.minor}}   Create a factor structure for nvar variables defined by nfact major factors and nvar/2 ``minor" factors.  

\code{\link{sim.anova}}    Simulate a 3 way balanced ANOVA or linear model, with or without repeated measures. 



}
\references{Revelle, W. (in preparation) An Introduction to Psychometric Theory with applications in R. Springer. at \url{http://personality-project.org/r/book/}  }

\author{William Revelle}

\seealso{ See above}
\examples{
simplex <- sim.simplex() #create the default structure
round(simplex,2)  #the correlation matrix



congeneric <- sim.congeneric()
round(congeneric,2)
R <- sim.hierarchical()
R
fx <- matrix(c(.9,.8,.7,rep(0,6),c(.8,.7,.6)),ncol=2)
fy <- c(.6,.5,.4)
Phi <- matrix(c(1,0,.5,0,1,.4,0,0,0),ncol=3)
R <- sim.structure(fx,Phi,fy) 
cor.plot(R$model) #show it graphically

simp <- sim.simplex()
#show the simplex structure using cor.plot
cor.plot(simp,colors=TRUE)
}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{multivariate}
\keyword{datagen}

