\name{factor.scores}
\alias{factor.scores}
\title{Various ways to estimate factor scores for the factor analysis model}

\description{A fundamental problem with factor analysis is that although the model is defined at the structural level, it is indeterminate at the data level. This problem of factor indeterminancy leads to alternative ways of estimating factor scores, none of which is ideal.  Following Grice (2001) four different methods are available here.
}
\usage{
factor.scores(x, f, Phi = NULL, method = c("Thurstone", "tenBerge", "Anderson", 
       "Bartlett", "Harman","components"),rho=NULL,missing=FALSE,impute="none")}
\arguments{
  \item{x}{Either a matrix of data if scores are to be found, or a correlation matrix if just the factor weights are to be found.}
  \item{f}{The output from the \code{\link{fa}}, \code{\link{irt.fa}},  \code{\link{factor.rotate}} functions, or a factor loading matrix.}
  \item{Phi}{If a pattern matrix is provided, then what were the factor intercorrelations.  Does not need to be specified if f is the output from the  \code{\link{fa}} or \code{\link{irt.fa}} functions.}
  \item{method}{Which of four factor score estimation procedures should be used. Defaults to "Thurstone" or regression based weights.  See details below for the other four methods.}
  \item{rho}{If x is a set of data and rho is specified, then find scores based upon the fa results and the correlations reported in rho.  Used when scoring fa.poly results.}
  \item{missing}{If missing is TRUE, missing items  are imputed using either the median or mean.  If missing is FALSE, the default, scores are found based upon the mean of the available items for each subject.}
  \item{impute}{By default, only complete cases are scored.  But, missing data can be imputed using "median" or "mean".  The number of missing by subject is reported.  If impute = "none", missing data are not scored.
}
}

\details{Although the factor analysis model is defined at the structural level, it is undefined at the data level.  This is a well known but little discussed problem with factor analysis.  

Factor scores represent estimates of common part of the variables and should not be thought of as identical to the factors themselves. If a factor is thought of as a chop stick stuck into the center of an ice cream cone and factor score estimates are represented by straws anywhere along the edge of the cone the problem of factor indeterminacy becomes clear, for depending on the shape of the cone, two straws can be negatively correlated with each other. (The imagery is taken from Niels Waller, adapted from Stanley Mulaik). In a very clear discussion of the problem of factor score indeterminacy, Grice (2001) reviews several alternative ways of estimating factor scores and considers weighting schemes that will produce uncorrelated factor score estimates as well as the effect of using coarse coded (unit weighted) factor weights.

\code{\link{factor.scores}} uses four different ways of estimate factor scores.  In all cases, the factor score estimates are based upon the data matrix, X, times a weighting matrix, W, which weights the observed variables.

For polytomous or dichotmous data, factor scores can be estimated using Item Response Theory techniques (e.g., using \code{link{irt.fa}} and then \code{link{scoreIrt}}.  Such scores are still just factor score estimates, for the IRT model is a latent variable model equivalent to factor analysis.  

\itemize{
\item method="Thurstone" finds the regression based weights: \eqn{W = R^{-1} F} where R is the correlation matrix and F is the factor loading matrix. 

\item method="tenBerge" finds weights such that the correlation between factors for an oblique solution is preserved. Note that  formula 8 in Grice has a typo in the formula for C and should be:
 \eqn{L = F \Phi^(1/2) }
 \eqn{C = R^(-1/2) L (L' R^(-1) L)^(-1/2) }
 \eqn{W = R ^(-1/2) C \Phi^(1/2) }
 
\item method="Anderson" finds weights such that the factor scores will be uncorrelated: \eqn{W = U^{-2}F (F' U^{-2} R  U^{-2} F)^{-1/2}} where U is the diagonal matrix of uniquenesses. The Anderson method works for orthogonal factors only, while the tenBerge method works for orthogonal or oblique solutions.

\item method = "Bartlett"  finds weights given \eqn{W = U^{-2}F (F' U^{-2}F)^{-1}}

 
\item method="Harman" finds weights based upon socalled "idealized" variables: \eqn{W =  F (t(F) F)^{-1}}.

\item method="components" uses weights that are just component loadings.  
}
}
\value{
\itemize{
\item scores (the factor scores if the raw data is given)
\item weights (the factor weights) 
\item r.scores  (The correlations of the factor score estimates.)
\item missing   A vector of the number of missing observations per subject
}
}
\references{
Grice, James W.,2001, Computing and evaluating factor scores,  Psychological Methods, 6,4, 430-450. (note the typo in equation 8)


 ten Berge, Jos M.F.,  Wim P. Krijnen, Tom Wansbeek and Alexander Shapiro (1999) Some new results on correlation-preserving factor scores prediction methods. Linear Algebra and its Applications, 289, 311-318.

Revelle, William. (in prep) An introduction to psychometric theory with applications in R. Springer.  Working draft available at \url{https://personality-project.org/r/book/} 

}
\author{
William Revelle
}




\seealso{
\code{\link{fa}}, \code{\link{factor.stats}}
}
\examples{
f3 <- fa(Thurstone,3 )
f3$weights  #just the scoring weights
f5 <- fa(bfi,5) #this does the factor analyis
my.scores <- factor.scores(bfi,f5, method="tenBerge")
#compare the tenBerge factor score correlation to the factor correlations
cor(my.scores$scores,use="pairwise") - f5$Phi  #compare to the f5$Phi values
#compare the default (regression) score correlations to the factor correlations
cor(f5$scores,use="pairwise")  - f5$Phi
#compare to the f5 solution

}

\keyword{ multivariate }
\keyword{ models}
