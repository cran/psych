\name{omegaStats}
\alias{omegaStats}

\title{Omega Statistics from a Factor Loading Matrix}
\description{Among the many ways of estimating internal consistency, the family of model based statistics known as omega (omega_hierarchical, omega_total, omega_subtest, omega_bifactor) provide a number of solutions. 
\code{\link{omegaStats}} provides common output for several different estimation procedures.    The \code{\link{omega}} function will find three versions of omega, as will the \code{\link{CFA}} function.  Because of their different assumptions, the omega_hierarchical and omega_bifactor will provide different estimates of the general factor saturation of a test. Users of other packages (e.g., lavaan, may find the basic statistics by providing a correlation matrix and factor loadings matrix. See the last example.  }
\usage{
omegaStats(r, loads = NULL, Phi = NULL, gload = NULL, group = NULL, h2 = NULL, 
    covar = FALSE, lavaan=NULL)
}

\arguments{
  \item{r}{A correlation matrix}
  \item{loads}{A factor loadings matrix with the first column representing a general factor.  Or alternatively users may specify the general factor and group factor loadings separately.}
  \item{Phi}{Factor intercorrelations}
  \item{gload}{A matrix of loadings on a general factor (if loads is not specified)}
  \item{group}{A matrix of group factor loadings (if loads is not specified)}
  \item{h2}{Communalities of the variables (r and loads are not specified)}
  \item{covar}{If TRUE, treat the r matrix as a covariance matrix to find u2}
  \item{lavaan}{To find omega from lavaan output, specify the lavaan object (see example)}
}

\details{Model based estimates of internal consistency include the family of estimates known as omega.  Omega_hierarchical (found by the \code{\link{omega}} function) estimates the general factor saturation of a test, while  \code{\link{CFA}} finds the amount of variance accounted by the first factor in a bifactor solution.  Although both procedures agree in terms of their estimates for the total amount of reliable variance, they disagree with respect to the amount that is common to all items.  

lavaan based estimates differ yet again (see the final example which is not run).
}
\value{

\item{omega hierarchical}{The \eqn{\omega_h} coefficient}
  \item{omega.lim}{The limit of \eqn{\omega_h} as the test becomes infinitly large}
  \item{omega total}{The \eqn{\omega_t} coefficient}
\item{stats}{as found by \code{\link{fa.stats}} }

\item{...}{Many other statistics as reported by \code{\link{fa}} }
\item{Call}{echoes the call to the function}
}
\references{
Jensen, Arthur R. and Li-Jen  Weng (1994) What is a good g?  Intelligence, 18, 3, 231-258

Revelle, William. (in prep) An introduction to psychometric theory with applications in R. Springer.  Working draft available at \url{https://personality-project.org/r/book/} 

Revelle, W. (1979).  Hierarchical cluster analysis and the internal structure of tests. Multivariate Behavioral Research, 14, 57-74. (\url{https://personality-project.org/revelle/publications/iclust.pdf})

Revelle, W. and Condon, D.M. (2019) Reliability from alpha to omega: A tutorial.  Psychological Assessment, 31, 12, 1395-1411. https://doi.org/10.1037/pas0000754.  \url{https://osf.io/preprints/psyarxiv/2y3w9} Preprint available from PsyArxiv 

Revelle, W. and Zinbarg, R. E. (2009) Coefficients alpha, beta, omega and the glb: comments on Sijtsma.  Psychometrika, 74, 1, 145-154. (\url{https://personality-project.org/revelle/publications/rz09.pdf}

Waller, N. G. (2017) Direct {Schmid-Leiman} Transformations and Rank-Deficient Loadings Matrices.  Psychometrika.  DOI: 10.1007/s11336-017-9599-0

 Zinbarg, R.E., Revelle, W., Yovel, I., & Li. W.  (2005). Cronbach's Alpha, Revelle's Beta, McDonald's Omega: Their relations with each and two alternative conceptualizations of reliability. Psychometrika. 70, 123-133.  \url{https://personality-project.org/revelle/publications/zinbarg.revelle.pmet.05.pdf}

Zinbarg, R., Yovel, I. & Revelle, W.  (2007).  Estimating omega  for structures containing two group factors:  Perils and prospects.  Applied Psychological Measurement. 31 (2), 135-157.
 
 Zinbarg, R., Yovel, I., Revelle, W. & McDonald, R. (2006).  Estimating generalizability to a universe of indicators that all have one attribute in common:  A comparison of estimators for omega.  Applied Psychological Measurement, 30, 121-144. DOI: 10.1177/0146621605278814


}
\author{William Revelle
}
\note{What is interesting is finding omega from the bifactor solution of non-correlated group factors.  Direct CFA.bifactor will find a general factor (even though not there), although conventional fa procedures will fail to find a general factor.

}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{fa}}, \code{\link{CFA}}, \code{\link{omega}}

}
\examples{
#Create a 3 correlated factor model using default values and 10000 observations
v9 <- sim.hierarchical(n=10000)

model <- 'F1=~ V1 + V2 + V3
          F2=~ V4 + V5 + V6
          F3 =~ V7 +V8 + V9'
test <-  CFA.bifactor(model,v9$observed)
test   #show the results
omegaStats(test$r,test$loadings)


om <- omega(v9$observed)
omegaStats(r=om$R,loads=om$schmid$sl[,1:4])

##Now creat a three factor model with orthogonal factors
model <- 'F1=~ V1 + V2 + V3
          F2=~ V4 + V5 + V6
          F3 =~ V7 +V8 + V9'



model9 <- 'F1 =~ .9*V1 + .8*V2 + .7*V3
           F2 =~ .8*V4 + .7*V5 +.6*V6
           F3 =~ .7*V7 + .6*V8 +.5*V9'

set.seed(42)          
v9s <- sim(model9,n=500)
 test <- CFA(model9,v9s$observed )  #do a cfa using Lavaan syntax
 test.bi <- CFA.bifactor(model9,v9s$observed)
 test.hi <- CFA.bifactor(model9,v9s$observed,g=TRUE)
 om <- omega(v9s$observed)
 omegaStats(test.bi$r,test.bi$loadings)

###do not run the following test of lavaan
if(FALSE) {
lav.model <- 'g=~ V1 + V2 + V3+ V4 + V5 + V6+V7 +V8 + V9
          F1=~ V1 + V2 + V3
          F2=~ V4 + V5 + V6
          F3 =~ V7 +V8 + V9'
library(lavaan)
fit <- cfa(lav.model,data=v9$observed,orthogonal=TRUE,std.lv=TRUE, estimator="ULS")
omegaStats(r=v9$r,lavaan=fit)   #from lavaan
fa.congruence(fit,test) 
}  #end of lavaan test

}
\keyword{ multivariate }
\keyword{ models }
