\name{partial.r}
\alias{partial.r}

\title{ Find the partial correlations for a set (x) of variables with set (y) removed. }
\description{A straightforward application of matrix algebra to remove the effect of the variables in the y set from the x set. Input may be either a data matrix or a correlation matrix.  Variables in x and y are specified by location.
}
\usage{
partial.r(m, x, y)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{m}{A data or correlation matrix}
  \item{x}{The variable numbers associated with the X set. }
  \item{y}{The variable numbers associated with the Y set}
  
}
\details{It is sometimes convenient to partial the effect of a number of variables (e.g., sex, age, education) out of the correlations of another set of variables.  This could be done laboriously by finding the residuals of various multiple correlations, and then correlating these residuals.  The matrix algebra alternative is to do it directly. 
To find the confidence intervals and "significance" of the correlations, use the \code{\link{corr.p}} function with n = n - s where s is the numer of covariates. 
}
\value{The matrix of partial correlations.}
\references{ Revelle, W. (in prep) An introduction to psychometric theory with applications in R. To be published by Springer.  (working draft available at  \url{http://personality-project.org/r/book/} 
}
\author{ William Revelle }

\seealso{ \code{\link{mat.regress}} for a similar application for regression}
\examples{
jen <- make.hierarchical()    #make up a correlation matrix 
round(jen[1:5,1:5],2)
par.r <- partial.r(jen,c(1,3,5),c(2,4))
cp <- corr.p(par.r,n=98)  #assumes the jen data based upon n =100.
print(cp,short=FALSE)  #show the confidence intervals as well
}
\keyword{multivariate}

