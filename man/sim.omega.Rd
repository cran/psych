\name{sim.omega}
\alias{sim.omega}
\alias{sim.general}
\alias{sim.parallel}
\title{Further functions to simulate psychological/psychometric data.
}
\description{A number of functions in the psych package will generate simulated data with particular structures.  These functions include
\code{\link{sim}} for a factor simplex, and \code{\link{sim.simplex}} for a data simplex, \code{\link{sim.circ}} for a circumplex structure, \code{\link{sim.congeneric}} for a one factor factor congeneric model, \code{\link{sim.dichot}} to simulate dichotomous items, \code{\link{sim.hierarchical}} to create a hierarchical factor model, \code{\link{sim.item}} a more general item simulation,
\code{\link{sim.minor}} to simulate major and minor factors,
\code{\link{sim.omega}} to test various examples of omega,
\code{\link{sim.parallel}} to compare the efficiency of various ways of deterimining the number of factors.   See the help pages for some more simulation functions here: 
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
sim.omega(nvar = 12, nfact = 3, n = 500, g = NULL, sem = FALSE, fbig = NULL, 
     fsmall = c(-0.2, 0.2),bipolar = TRUE, om.fact = 3, flip = TRUE,
      option = "equal", ntrials = 10)
sim.parallel(ntrials=10,nvar = c(12,24,36,48),nfact = c(1,2,3,4,6),
            n = c(200,400)) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
   \item{nvar}{Number of variables}
  \item{nfact}{Number of group factors in sim.general, sim.omega}
  \item{n}{ Number of cases to simulate.  If n=0 or NULL, the population matrix is returned.}
  \item{g}{General factor correlations in sim.general and general factor loadings in sim.omega and sim.minor}
  \item{sem}{Should the sim.omega function do both an EFA omega as well as a CFA omega using the sem package?}
 \item{fbig}{Factor loadings for the main factors.  Default is a simple structure with loadings sampled from (.8,.6) for nvar/nfact variables and 0 for the remaining.  If fbig is specified, then  each factor has loadings sampled from it.}
   \item{fsmall}{nvar/2 small factors are generated with loadings sampled from (-.2,0,.2)}
 \item{bipolar}{if TRUE, then positive and negative loadings are generated from fbig}
 \item{om.fact}{Number of  factors to extract in omega}
 \item{flip}{In omega, should item signs be flipped if negative}
   \item{option}{In omega, for the case of two factors, how to weight them?}
 \item{ntrials}{Number of replications per level}
}
\details{Simulation of data structures is a very useful tool in psychometric research and teaching.  By knowing ``truth" it is possible to see how well various algorithms can capture it.  For a much longer discussion of the use of simulation in psychometrics, see the accompany vignettes.  

The simulations documented here are a miscellaneous set of functions that will be documented in other help files eventually.

The default values for \code{\link{sim.structure}} is to generate a 4 factor, 12 variable data set with a simplex structure between the factors. This, and the simplex of items (\code{\link{sim.simplex}}) can also be converted in a STARS model with an autoregressive component (alpha) and a stable trait component (lambda). 

Two data structures that are particular challenges to exploratory factor analysis are the simplex structure and the presence of minor factors.  Simplex structures \code{\link{sim.simplex}} will typically occur in developmental or learning contexts and have a correlation structure of r between adjacent variables and r^n for variables n apart.  Although just one latent variable (r) needs to be estimated, the structure will have nvar-1 factors.  

An alternative version of the simplex is the State-Trait-Auto Regressive Structure (STARS) which has both a simplex state structure, with autoregressive path alpha and a trait structure with path lambda. This simulated in  \code{\link{sim.simplex}} by specifying a non-zero lambda value.

Many simulations of factor structures assume that except for the major factors, all residuals are normally distributed around 0.  An alternative, and perhaps more realistic situation, is that the there are a few major (big) factors and many minor (small) factors.  The challenge is thus to identify the major factors. \code{\link{sim.minor}} generates such structures.  The structures generated can be thought of as havinga a major factor structure with some small correlated residuals. To make these simulations complete, the possibility of a general factor is considered.  For simplicity, sim.minor allows one to specify a set of loadings to be sampled from for g, fmajor and fminor.  Alternatively, it is possible to specify the complete factor matrix.

Another structure worth considering is direct modeling of a general factor with several group factors.  This is done using \code{\link{sim.general}}.

Although coefficient \eqn{\omega}{\omega} is a very useful indicator of the general factor saturation of a unifactorial test (one with perhaps several sub factors), it has problems with the case of multiple, independent factors.  In this situation, one of the factors is labelled as ``general'' and  the omega estimate is too large.  This situation may be explored using the \code{\link{sim.omega}} function with general left as NULL.  If there is a general factor, then results from \code{\link{sim.omega}} suggests that omega estimated either from EFA or from SEM does a pretty good job of identifying it but that the EFA approach using Schmid-Leiman transformation is somewhat more robust than the SEM approach. 

The four irt simulations, sim.rasch, sim.irt, sim.npl and sim.npn, simulate dichotomous items following the Item Response model.  sim.irt just calls either sim.npl (for logistic models) or sim.npn (for normal models) depending upon the specification of the model. 

The logistic model is \deqn{P(i,j) = \gamma + \frac{\zeta-\gamma}{1+ e^{\alpha(\delta-\theta)}}}{P(i,j) = \gamma + (\zeta-\gamma)/(1+ exp(\alpha(\delta-\theta)))} where \eqn{\gamma} is the lower asymptote or guesssing parameter, \eqn{\zeta} is the upper asymptote (normally 1), \eqn{\alpha} is item discrimination and \eqn{\delta} is item difficulty.  For the 1 Paramater Logistic (Rasch) model, gamma=0, zeta=1, alpha=1 and item difficulty is the only free parameter to specify.

For the 2PL and 2PN models, a = \eqn{\alpha} and  d = \eqn{\delta} are specified. \cr
For the 3PL or 3PN models, items also differ in their guessing parameter c =\eqn{\gamma}. \cr
For the 4PL and 4PN models, the upper asymptote, z= \eqn{\zeta} is also specified.  \cr
(Graphics of these may be seen in the demonstrations for the \code{\link{logistic}} function.)

The normal model (irt.npn calculates the probability using pnorm instead of the logistic function used in irt.npl, but the meaning of the parameters are otherwise the same.  With the a = \eqn{\alpha} parameter = 1.702 in the logistic model the two models are practically identical.

In parallel to the dichotomous IRT simulations are the poly versions which simulate polytomous item models.  They have the additional parameter of how many categories to simulate.  In addition, the \code{\link{sim.poly.ideal}} functions will simulate an ideal point or unfolding model in which the response probability varies by the distance from each subject's ideal point.  Some have claimed that this is a more appropriate model of the responses to personality questionnaires.  It will lead to simplex like structures which may be fit by a two factor model.  The middle items form one factor, the extreme a bipolar factor.

By default, the theta parameter is created in each function as normally distributed with mean mu=0  and sd=1.  In the case where you want to specify the theta to be equivalent from another simulation or fixed for a particular experimental condition, either take the theta object from the output of a previous simulation, or create it using whatever properties are desired. 

The previous functions all assume one latent trait.  Alternatively, we can simulate dichotomous or polytomous items with a particular structure using the sim.poly.mat function.  This takes as input the population correlation matrix, the population marginals, and the sample size.  It returns categorical items with the specified structure.

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


\code{\link{sim.multilevel}}  To understand some of the basic concepts of multilevel modeling, it is useful to create multilevel structures.  The correlations of aggregated data is sometimes called an 'ecological correlation'.  That group level and individual level correlations are independent makes such inferences problematic.  This simulation allows for demonstrations that correlations within groups do not imply, nor are implied by, correlations between group means. 

}
\value{\code{link{sim.parallel}} returns the results of multiple simulations of fa.parallel for various combinations of the number of variables, numbers of factors, and sample size.
  \code{link{sim.general}} returns either a correlation matrix (if n is not specified) or a data matrix with a general factor.
}
\references{ \url{https://personality-project.org/r/r.omega.html} \cr

Revelle, William. (in prep) An introduction to psychometric theory with applications in R. Springer.  Working draft available at \url{https://personality-project.org/r/book/} 

Revelle, W. (1979).  Hierarchical cluster analysis and the internal structure of tests. Multivariate Behavioral Research, 14, 57-74. (\url{https://personality-project.org/revelle/publications/iclust.pdf})

Revelle, W. and Zinbarg, R. E. (2009) Coefficients alpha, beta, omega and the glb: comments on Sijtsma.  Psychometrika, 74, 1, 145-154. (\url{https://personality-project.org/revelle/publications/rz09.pdf}

Waller, N. G. (2017) Direct {Schmid-Leiman} Transformations and Rank-Deficient Loadings Matrices.  Psychometrika.  DOI: 10.1007/s11336-017-9599-0

 Zinbarg, R.E., Revelle, W., Yovel, I., & Li. W.  (2005). Cronbach's Alpha, Revelle's Beta, McDonald's Omega: Their relations with each and two alternative conceptualizations of reliability. Psychometrika. 70, 123-133.  \url{https://personality-project.org/revelle/publications/zinbarg.revelle.pmet.05.pdf}

Zinbarg, R., Yovel, I. & Revelle, W.  (2007).  Estimating omega  for structures containing two group factors:  Perils and prospects.  Applied Psychological Measurement. 31 (2), 135-157.
 
 Zinbarg, R., Yovel, I., Revelle, W. & McDonald, R. (2006).  Estimating generalizability to a universe of indicators that all have one attribute in common:  A comparison of estimators for omega.  Applied Psychological Measurement, 30, 121-144. DOI: 10.1177/0146621605278814 \url{https://journals.sagepub.com/doi/10.1177/0146621605278814}
}
\author{ 
 \url{https://personality-project.org/revelle.html} \cr
Maintainer: William Revelle  \email{  revelle@northwestern.edu  } 
}
\note{A few  of the various sim functions are included here.
}



\seealso{\code{\link{omega}}
}
\examples{
#test <- sim.omega()
  }


\keyword{ multivariate}
\keyword{datagen}

