\name{sim}
\alias{sim}
\alias{sim.simplex}
\alias{sim.minor}
\alias{sim.omega}
\alias{sim.general}
\alias{sim.parallel}
\alias{sim.rasch}
\alias{sim.irt}
\alias{sim.npl}
\alias{sim.npn}
\alias{sim.poly}
\alias{sim.poly.npl}
\alias{sim.poly.npn}
\alias{sim.poly.ideal}
\alias{sim.poly.ideal.npl}
\alias{sim.poly.ideal.npn}
\title{Functions to simulate psychological/psychometric data.}
\description{A number of functions in the psych package will generate simulated data with particular structures.  These functions include
\code{\link{sim}} for a factor simplex, and \code{\link{sim.simplex}} for a data simplex, \code{\link{sim.circ}} for a circumplex structure, \code{\link{sim.congeneric}} for a one factor factor congeneric model, \code{\link{sim.dichot}} to simulate dichotomous items, \code{\link{sim.hierarchical}} to create a hierarchical factor model, \code{\link{sim.item}} a more general item simulation,
\code{\link{sim.minor}} to simulate major and minor factors,
\code{\link{sim.omega}} to test various examples of omega,
\code{\link{sim.parallel}} to compare the efficiency of various ways of deterimining the number of factors,
\code{\link{sim.rasch}} to create simulated rasch data, 
\code{\link{sim.irt}} to create general 1 to 4 parameter IRT data by calling 
\code{\link{sim.npl}} 1 to 4 parameter logistic IRT or 
\code{\link{sim.npn}} 1 to 4 paramater normal IRT,
\code{\link{sim.poly}} to create polytomous ideas by calling
\code{\link{sim.poly.npn}} 1-4 parameter polytomous normal theory items or
\code{\link{sim.poly.npl}} 1-4 parameter polytomous logistic items, and 
\code{\link{sim.poly.ideal}} which creates data following an ideal point or unfolding model by calling 
\code{\link{sim.poly.ideal.npn}} 1-4 parameter polytomous normal theory ideal point model or 
\code{\link{sim.poly.ideal.npl}} 1-4 parameter polytomous logistic ideal point model.

\code{\link{sim.structural}} a general simulation of structural models,  and \code{\link{sim.anova}} for ANOVA and lm simulations, and \code{\link{sim.VSS}}. Some of these functions are separately documented and are listed here for ease of the help function.  See each function for more detailed help.
}
\usage{
sim(fx=NULL,Phi=NULL,fy=NULL,alpha=.8,lambda = 0,n=0,mu=NULL,raw=TRUE)
sim.simplex(nvar =12, alpha=.8,lambda=0,beta=1,mu=NULL, n=0)
sim.general(nvar=9,nfact =3, g=.3,r=.3,n=0) 
sim.minor(nvar=12,nfact=3,n=0,fbig=NULL,fsmall = c(-.2,.2),bipolar=TRUE) 
sim.omega(nvar=12,nfact=3,n=0,fbig=NULL,fsmall = c(-.2,.2),bipolar=TRUE,om.fact=3,flip=TRUE,option="equal",ntrials=10)
sim.parallel(ntrials=10,nvar = c(12,24,36,48),nfact = c(1,2,3,4,6),
n = c(200,400)) 
sim.rasch(nvar = 5,n = 500, low=-3,high=3,d=NULL, a=1,mu=0,sd=1)
sim.irt(nvar = 5, n = 500, low=-3, high=3,a=NULL,c=0,z=1,d=NULL,mu=0,sd=1,mod="logistic")
sim.npl(nvar = 5, n = 500, low=-3,high=3,a=NULL,c=0,z=1,d=NULL,mu=0,sd=1)
sim.npn(nvar = 5, n = 500, low=-3,high=3,a=NULL,c=0,z=1,d=NULL,mu=0,sd=1)
sim.poly(nvar = 5 ,n = 500,low=-2,high=2,a=NULL,c=0,z=1,d=NULL, mu=0,sd=1,cat=5,mod="logistic") 
sim.poly.npn(nvar = 5 ,n = 500,low=-2,high=2,a=NULL,c=0,z=1,d=NULL, mu=0,sd=1,cat=5) 
sim.poly.npl(nvar = 5 ,n = 500,low=-2,high=2,a=NULL,c=0,z=1,d=NULL, mu=0,sd=1,cat=5) 
sim.poly.ideal(nvar = 5 ,n = 500,low=-2,high=2,a=NULL,c=0,z=1,d=NULL, mu=0,sd=1,cat=5,mod="logistic") 
sim.poly.ideal.npn(nvar = 5 ,n = 500,low=-2,high=2,a=NULL,c=0,z=1,d=NULL, mu=0,sd=1,cat=5) 
sim.poly.ideal.npl(nvar = 5 ,n = 500,low=-2,high=2,a=NULL,c=0,z=1,d=NULL, mu=0,sd=1,cat=5,theta=NULL) 

}
\arguments{
  \item{fx}{The measurement model for x. If NULL, a 4 factor model is generated}
  \item{Phi}{The structure matrix of the latent variables}
  \item{fy}{The measurement model for y}
  \item{mu}{The means structure for the fx factors}
  \item{n}{ Number of cases to simulate.  If n=0 or NULL, the population matrix is returned.}
  \item{raw}{if raw=TRUE, raw data are returned as well.}
  \item{nvar}{Number of variables for a simplex structure}
  \item{nfact}{Number of large factors to simulate in sim.minor,number of group factors in sim.general,sim.omega}
  \item{g}{General factor correlations in sim.general}
  \item{r}{group factor correlations in sim.general}
  \item{alpha}{the base correlation for an autoregressive simplex}
  \item{lambda}{the trait component of a State Trait Autoregressive Simplex}
  \item{beta}{Test reliability of a STARS simplex}
  \item{fbig}{Factor loadings for the main factors.  Default is a simple structure with loadings sampled from (.8,.6) for nvar/nfact variables and 0 for the remaining.  If fbig is specified, then  each factor has loadings sampled from it.}
  \item{bipolar}{if TRUE, then positive and negative loadings are generated from fbig}
  \item{om.fact}{Number of  factors to extract in omega}
  \item{flip}{In omega, should item signs be flipped if negative}
  \item{option}{In omega, for the case of two factors, how to weight them?}
  \item{fsmall}{nvar/2 small factors are generated with loadings sampled from (-.2,0,.2)}
  \item{ntrials}{Number of replications per level}
  \item{low}{lower difficulty for sim.rasch or sim.irt}
  \item{high}{higher difficulty for sim.rasch or sim.irt}
  \item{a}{if not specified as a vector, the descrimination parameter a = \eqn{\alpha} will be set to 1.0 for all items}
  \item{d}{ if not specified as a vector, item difficulties (d = \eqn{\delta}) will range from low to high}
  \item{c}{the gamma parameter: if not specified as a vector, the guessing asymptote is set to 0}
  \item{z}{the zeta parameter: if not specified as a vector, set to 1} 
  \item{sd}{the standard deviation for the underlying latent variable in the irt simulations}
  \item{mod}{which IRT model to use, mod="logistic" simulates a logistic function, otherwise, a normal function}
  \item{cat}{Number of categories to simulate in sim.poly.  If cat=2, then this is the same as simulating t/f items and sim.poly is functionally equivalent to sim.irt}
  \item{theta}{The underlying latent trait value for each simulated subject}
   }

\details{Simulation of data structures is a very useful tool in psychometric research and teaching.  By knowing ``truth" it is possible to see how well various algorithms can capture it.  For a much longer discussion of the use of simulation in psychometrics, see the accompany vignettes.  

The simulations documented here are a miscellaneous set of functions that will be documented in other help files eventually.

The default values for \code{\link{sim.structure}} is to generate a 4 factor, 12 variable data set with a simplex structure between the factors. This, and the simplex of items (\code{\link{sim.simplex}}) can also be converted in a STARS model with an autoregressive component (alpha) and a stable trait component (lambda). 

Two data structures that are particular challenges to exploratory factor analysis are the simplex structure and the presence of minor factors.  Simplex structures \code{\link{sim.simplex}} will typically occur in developmental or learning contexts and have a correlation structure of r between adjacent variables and r^n for variables n apart.  Although just one latent variable (r) needs to be estimated, the structure will have nvar-1 factors.  

An alternative version of the simplex is the State-Trait-Auto Regressive Structure (STARS) which has both a simplex state structure, with autoregressive path alpha and a trait structure with path lambda. This simulated in  \code{\link{sim.simplex}} by specifying a non-zero lambda value.

Many simulations of factor structures assume that except for the major factors, all residuals are normally distributed around 0.  An alternative, and perhaps more realistic situation, is that the there are a few major (big) factors and many minor (small) factors.  The challenge is thus to identify the major factors. \code{\link{sim.minor}} generates such structures.  The structures generated can be thought of as havinga a major factor structure with some small correlated residuals.  

Another structure worth considering is direct modeling of a general factor with several group factors.  This is done using \code{\link{sim.general}}.

Although coefficient \eqn{\omega}{\omega} is a very useful indicator of the general factor saturation of a unifactorial test (one with perhaps several sub factors), it has problems with the case of multiple, independent factors.  In this situation, one of the factors is labelled as ``general'' and  the omega estimate is too large.  This situation may be explored using the \code{\link{sim.omega}} function. 

The four irt simulations, sim.rasch, sim.irt, sim.npl and sim.npn, simulate dichotomous items following the Item Response model.  sim.irt just calls either sim.npl (for logistic models) or sim.npn (for normal models) depending upon the specification of the model. 

The logistic model is \deqn{P(i,j) = \gamma + \frac{\zeta-\gamma}{1+ e^{\alpha(\delta-\theta)}}}{P(i,j) = \gamma + (\zeta-\gamma)/(1+ exp(\alpha(\delta-\theta)))} where \eqn{\gamma} is the lower asymptote or guesssing parameter, \eqn{\zeta} is the upper asymptote (normally 1), \eqn{\alpha} is item discrimination and \eqn{\delta} is item difficulty.  For the 1 Paramater Logistic (Rasch) model, gamma=0, zeta=1, alpha=1 and item difficulty is the only free parameter to specify.

For the 2PL and 2PN models, a = \eqn{\alpha} and  d = \eqn{\delta} are specified. \cr
For the 3PL or 3PN models, items also differ in their guessing parameter c =\eqn{\gamma}. \cr
For the 4PL and 4PN models, the upper asymptote, z= \eqn{\zeta} is also specified.  \cr
(Graphics of these may be seen in the demonstrations for the \code{\link{logistic}} function.)

The normal model (irt.npn calculates the probability using pnorm instead of the logistic function used in irt.npl, but the meaning of the parameters are otherwise the same.  With the a = \eqn{\alpha} parameter = 1.702 in the logistic model the two models are practically identical.

In parallel to the dichotomous IRT simulations are the poly versions which simulate polytomous item models.  They have the additional parameter of how many categories to simulate.  In addition, the \code{\link{sim.poly.ideal}} functions will simulate an ideal point or unfolding model in which the response probability varies by the distance from each subject's ideal point.  Some have claimed that this is a more appropriate model of the responses to personality questionnaires.  It will lead to simplex like structures which may be fit by a two factor model.  The middle items form one factor, the extreme a bipolar factor.

Other simulation functions in psych are:

\code{\link{sim.structure}}  A function to combine a measurement and structural model into one data matrix.  Useful for understanding structural equation models.  Combined with \code{\link{structure.diagram}} to see the proposed structure.  


\code{\link{sim.congeneric}}   A function to create congeneric items/tests for demonstrating classical test theory. This is just a special case of sim.structure.
 
\code{\link{sim.hierarchical}}  A function to create data with a hierarchical (bifactor) structure.  

\code{\link{sim.item}}      A function to create items that either have a simple structure or a circumplex structure.

\code{\link{sim.circ}}    Create data with a circumplex structure.

\code{\link{sim.dichot}}    Create dichotomous item data with a simple or circumplex structure.


\code{\link{sim.minor}}   Create a factor structure for nvar variables defined by nfact major factors and nvar/2 ``minor" factors for n observations.  

Although the standard factor model assumes that K major factors (K << nvar) will account for the correlations among the variables

\deqn{R = FF' + U^2} 
where R is of rank P and F is a P x K matrix of factor coefficients and U is a diagonal matrix of uniquenesses.  However, in many cases, particularly when working with items, there are many small factors (sometimes referred to as correlated residuals) that need to be considered as well.  This leads to a data structure such that 
\deqn{R = FF' + MM' + U^2} 
where R is a P x P matrix of correlations, F is a  P x K factor loading matrix,  M is a P x P/2 matrix of minor factor loadings, and U is a diagonal matrix (P x P) of uniquenesses.  

Such a correlation matrix will have a poor \eqn{\chi^2} value in terms of goodness of fit if just the K factors are extracted, even though for all intents and purposes, it is well fit.  

sim.minor will generate such data sets with big factors with loadings of .6 to .8 and small factors with loadings of -.2 to .2.  These may both be adjusted.

\code{\link{sim.parallel}} Create a number of simulated data sets using sim.minor to show how parallel analysis works.  The general observation is that with the presence of minor factors, parallel analysis is probably best done with component eigen values rather than factor eigen values, even when using the factor model. 

\code{\link{sim.anova}}    Simulate a 3 way balanced ANOVA or linear model, with or without repeated measures. Useful for teaching research  methods and generating teaching examples. 



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
cor.plot(simp,colors=TRUE,main="A simplex structure")
#Show a STARS model 
simp <- sim.simplex(alpha=.8,lambda=.8)
#show the simplex structure using cor.plot
cor.plot(simp,colors=TRUE,main="State Trait Auto Regressive Simplex" )
}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{multivariate}
\keyword{datagen}

