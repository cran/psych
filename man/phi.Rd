\name{phi}
\alias{phi}

\title{ Find the phi coefficient of correlation between two dichotomous variables }
\description{Given a 1 x 4 vector or a 2 x 2 matrix of frequencies, find the phi coefficient of correlation.  Typical use is in the case of predicting a dichotomous criterion from a dichotomous predictor.
}
\usage{
phi(t, digits = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{t}{a 1 x 4 vector or a 2 x 2 matrix }
  \item{digits}{ round the result to digits }
}
\details{In many prediction situations, a dichotomous predictor (accept/reject) is validated against a dichotomous criterion (success/failure).  Although a polychoric correlation estimates the underlying Pearson correlation as if the predictor and criteria were continuous and bivariate normal variables, and the tetrachoric correlation if both x and y are assumed to dichotomized normal distributions,  the phi coefficient is the Pearson applied to a matrix of 0's and 1s. 

The phi coefficient was first reported by Yule (1912), but should not be confused with the \code{\link{Yule}} Q coefficient. 

For a very useful discussion of various measures of association given a 2 x 2 table, and why one should probably prefer the \code{\link{Yule}} Q coefficient, see Warren (2008). 

Given a two x two table of counts \cr
\tabular{llll}{
\tab a \tab b \tab a+b (R1)\cr
\tab c \tab d  \tab c+d (R2)\cr
\tab a+c(C1) \tab b+d (C2) \tab a+b+c+d (N)
}
convert all counts to fractions of the total and then \\
Phi = [a- (a+b)*(a+c)]/sqrt((a+b)(c+d)(a+c)(b+d) ) =\\
(a - R1 * C1)/sqrt(R1 * R2 * C1 * C2)

This is in contrast to the Yule coefficient, Q, where \\
Q = (ad - bc)/(ad+bc) which is the same as \\
 [a- (a+b)*(a+c)]/(ad+bc)

Since the phi coefficient is just a Pearson correlation applied to dichotomous data, to find a matrix of phis from a data set involves just finding the correlations using cor  or \code{\link{lowerCor}} or \code{\link{corr.test}}.  
}
\value{phi coefficient of correlation
}

\author{William Revelle with modifications by Leo Gurtler }

\references{Warrens, Matthijs (2008), On Association Coefficients for 2x2 Tables and Properties That Do Not Depend on the Marginal Distributions. Psychometrika, 73, 777-789.

Yule, G.U. (1912). On the methods of measuring the association between two attributes. Journal of the Royal Statistical Society, 75, 579-652.}

\seealso{ \code{\link{phi2tetra}} ,\code{\link{Yule}}, \code{\link{Yule.inv}} \code{\link{Yule2phi}}, \code{\link{tetrachoric}} and \code{\link{polychoric}}}
\examples{
phi(c(30,20,20,30))
phi(c(40,10,10,40))
x <- matrix(c(40,5,20,20),ncol=2)
phi(x)


}
\keyword{multivariate }% at least one, from doc/KEYWORDS
\keyword{models }% __ONLY ONE__ keyword per line
