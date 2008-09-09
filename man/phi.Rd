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
\details{In many prediction situations, a dichotomous predictor (accept/reject) is validated against a dichotomous criterion (success/failure).  Although a polychoric correlation estimates the underlying Pearson correlation as if the predictor and criteria were continuous and bivariate normal variables, the phi coefficient is the Pearson applied to a matrix of 0's and 1s.  

Given a two x two table of counts \cr
\tabular{llll}{
\tab a \tab b \tab a+b\cr
\tab c \tab d  \tab c+d \cr
\tab a+c \tab b+d \tab a+b+c+d
}
convert all counts to fractions of the total and then \\
Phi = a- (a+b)*(a+c)/sqrt((a+b)(c+d)(a+c)(b+d) ) 


}
\value{phi coefficient of correlation
}

\author{William Revelle with modifications by Leo Gurtler }

\seealso{ \code{\link{phi2poly}} ,\code{\link{Yule}}, \code{\link{Yule2phi}}}
\examples{
phi(c(30,20,20,30))
phi(c(40,10,10,40))
x <- matrix(c(40,5,20,20),ncol=2)
phi(x)


}
\keyword{multivariate }% at least one, from doc/KEYWORDS
\keyword{models }% __ONLY ONE__ keyword per line
