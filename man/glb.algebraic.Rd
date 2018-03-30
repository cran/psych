\name{glb.algebraic}
\alias{glb.algebraic}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Find the greatest lower bound to reliability.
}
\description{
The greatest lower bound solves the ``educational testing problem". That is, what is the reliability of a test? (See \code{\link{guttman}} for a discussion of the problem). Although there are many estimates of a test reliability (Guttman, 1945) most underestimate the true reliability of a test.

For a given covariance matrix of items, C, the function finds the greatest lower bound to reliability of the total score using the csdp function from the Rcsdp package.}

\usage{
glb.algebraic(Cov, LoBounds = NULL, UpBounds = NULL)
}

\arguments{
 \item{Cov}{A p * p covariance matrix. Positive definiteness is not checked.}
 \item{LoBounds}{A vector \eqn{l =(l_1, \dots, l_p)}{L = (l1 ... lp)} of length p with lower bounds to the diagonal elements \eqn{x_i}.	The  default l=(0, . . . , 0) does not imply any constraint, because positive semidefiniteness of the matrix \eqn{\tilde{ C} + Diag(x)}{C0 + Diag(x)} implies \eqn{0 \leq x_i}{0 \le xi.}}
 \item{UpBounds}{A vector u =(u1, . . . , up) of length p with upper bounds to the diagonal elements xi. The default is u = v.}
}
\details{
If C is a p * p-covariance matrix, v = diag(C) its diagonal (i. e. the vector of variances \eqn{v_i = c_{ii}}), \eqn{\tilde { C} = C - Diag(v)}{C0 = C - Diag(v)} is the covariance matrix with 0s substituted in the diagonal and x = the vector \eqn{x_1, \dots ,x_n}{(x1, . . . , xp)} the educational testing problem is (see e. g., Al-Homidan 2008)

\deqn{\sum_{i=1}^p x_i \rightarrow \min}{(Sum i = 1 to p  xi) ->  min }

s.t. \deqn{\tilde{ C} + Diag(x) \geq 0}{C0 + Diag(x) >= 0}(i.e. positive semidefinite) and \eqn{x_i \leq v_i, i=1,\dots,p}{xi \le vi, i = 1 ..., p}. This is the same as minimizing the trace of the symmetric matrix
\deqn{\tilde{ C}+diag(x)=\left(\begin{array}{llll}
      x_1    & c_{12} & \ldots & c_{1p} \\
      c_{12} & x_2    & \ldots & c_{2p} \\
      \vdots & \vdots & \ddots & \vdots \\
      c_{1p} & c_{2p} & \ldots & x_p\\
      \end{array}\right)}{C0 + Diag(x)}


s. t. \eqn{\tilde{ C} + Diag(x)}{C0 + Diag(x)} is positive semidefinite and \eqn{x_i \leq v_i}{xi \le vi}. 

The greatest lower bound to reliability is
\deqn{\frac{\sum_{ij} \bar{c_{ij}} + \sum_i x_i}{\sum_{ij}c_{ij}}}{ (sum cij (i \ne j) + sum xi )/ sum cij}

Additionally, function glb.algebraic allows the user to  change the upper bounds \eqn{x_i \leq v_i}{xi \le vi} to
\eqn{x_i \leq u_i}{xi \le ui} and add lower bounds \eqn{l_i \leq x_i}{li \le xi}.

The greatest lower bound to reliability is applicable for tests with non-homogeneous items. It gives a sharp lower bound to the reliability of the total test score.

Caution: Though glb.algebraic gives exact lower bounds for exact covariance matrices, the estimates from empirical matrices may be strongly biased upwards for small and medium sample sizes.

glb.algebraic is wrapper for a call to function csdp of package Rcsdp (see its documentation).

If Cov is the covariance matrix of subtests/items with known lower bounds, rel, to their reliabilities (e. g. Cronbachs \eqn{\alpha}), LoBounds can be used to improve the lower bound to reliability by setting LoBounds <- rel*diag(Cov).

Changing UpBounds can be used to relax constraints \eqn{x_i \leq v_i}{xi \le vi} or to fix \eqn{x_i}{xi}-values by setting LoBounds[i] < -z; UpBounds[i] <- z.
}
\value{
\item{glb }{The algebraic greatest lower bound}
\item{solution}{The vector x of the solution of the semidefinite program. These are the elements on the diagonal of C.}
\item{status}{Status of the solution. See documentation of csdp in package Rcsdp. If status is 2 or greater or equal than 4, no glb and solution is returned. If status is not 0, a warning message is generated.}
\item{Call}{The calling string}
}
\references{Al-Homidan S (2008). Semidefinite programming for the educational testing problem. Central European Journal of Operations Research, 16:239-249.

Bentler PM (1972) A lower-bound method for the dimension-free measurement of internal consistency. Soc Sci Res 1:343-357.

Fletcher R (1981) A nonlinear programming problem in statistics (educational testing). SIAM J Sci Stat Comput 2:257-267.

Shapiro A, ten Berge JMF (2000). The asymptotic bias of minimum trace factor analysis, with applications to the greatest lower bound to reliability. Psychometrika, 65:413-425.

ten Berge, Socan G (2004). The greatest bound to reliability of a test and the hypothesis of unidimensionality. Psychometrika, 69:613-625.
}
\author{Andreas Moltner \cr 
Center of Excellence for Assessment in Medicine/Baden-Wurttemberg\cr
University of Heidelberg\cr 

William Revelle\cr 
Department of Psychology \cr
Northwestern University Evanston, Illiniois \cr https://personality-project.org/revelle.html

}


\seealso{
 For an alternative estimate of the greatest lower bound, see \code{\link{glb.fa}}.  For multiple estimates of reliablity, see \code{\link{guttman}}
}
\examples{

Cv<-matrix(c(215, 64, 33, 22,
              64, 97, 57, 25,
              33, 57,103, 36,
              22, 25, 36, 77),ncol=4)

Cv                    # covariance matrix of a test with 4 subtests
Cr<-cov2cor(Cv)       # Correlation matrix of tests
if(!require(Rcsdp)) {print("Rcsdp must be installed to find the glb.algebraic")} else {
 glb.algebraic(Cv)     # glb of total score
glb.algebraic(Cr)      # glb of sum of standardized scores

 w<-c(1,2,2,1)         # glb of weighted total score
 glb.algebraic(diag(w) \%*\% Cv \%*\% diag(w))  
alphas <- c(0.8,0,0,0) # Internal consistency of first test is known

glb.algebraic(Cv,LoBounds=alphas*diag(Cv))

                      # Fix all diagonal elements to 1 but the first:

lb <- glb.algebraic(Cr,LoBounds=c(0,1,1,1),UpBounds=c(1,1,1,1))
lb$solution[1]        # should be the same as the squared mult. corr.
smc(Cr)[1] 
}                        

}

\keyword{multivariate}

