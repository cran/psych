\name{partial.r}
\alias{partial.r}

\title{ Find the partial correlations for a set (x) of variables with set (y) removed. }
\description{A straightforward application of matrix algebra to remove the effect of the variables in the y set from the x set. Input may be either a data matrix or a correlation matrix.  Variables in x and y are specified by location.  If x and y are not specified, then the effect of all variables are partialled from all the other correlations.  
}
\usage{
partial.r(data, x, y,use="pairwise",method="pearson")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{A data or correlation matrix}
  \item{x}{The variable names or locations associated with the X set. }
  \item{y}{The variable  names or locations associated with the Y set}
  \item{use}{How should we treat missing data? The default is pairwise complete.}
  \item{method}{Which method of correlation should we use, the default is pearson.}
  
}
\details{There are two ways to use \code{\link{partial.r}}.  One is to find the complete partial correlation matrix (that is, partial all the other variables out of each variable).  This may be done by simply specify the raw data or correlation matrix.  (In the case of raw data, correlations will be found according to use and method.)  In this case, just specify the data matrix.  .

The second usage is to partial a set of variables(y) out of another set (x). It is sometimes convenient to partial the effect of a number of variables (e.g., sex, age, education) out of the correlations of another set of variables.  This could be done laboriously by finding the residuals of various multiple correlations, and then correlating these residuals.  The matrix algebra alternative is to do it directly. 
To find the confidence intervals and "significance" of the correlations, use the \code{\link{corr.p}} function with n = n - s where s is the numer of covariates. 
}
\value{The matrix of partial correlations.}
\references{ Revelle, W. (in prep) An introduction to psychometric theory with applications in R. To be published by Springer.  (working draft available at  \url{http://personality-project.org/r/book/} 
}
\author{ William Revelle }

\seealso{ \code{\link{setCor}} for a similar application for regression. \code{\link{lowerMat}} to neatly show a correlation matrix, and \code{\link{corr.p}} to find the confidence intervals of a correlation. }
\examples{
jen <- make.hierarchical()    #make up a correlation matrix 
lowerMat(jen[1:5,1:5])
par.r <- partial.r(jen,c(1,3,5),c(2,4))
lowerMat(par.r)
cp <- corr.p(par.r,n=98)  #assumes the jen data based upon n =100.
print(cp,short=FALSE)  #show the confidence intervals as well
 #partial all from all correlations.
lowerMat(partial.r(jen)) 
}
\keyword{multivariate}

