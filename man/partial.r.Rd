\name{partial.r}
\alias{partial.r}

\title{ Find the partial correlations for a set (x) of variables with set (y) removed. }
\description{A straightforward application of matrix algebra to remove the effect of the variables in the y set from the x set. Input may be either a data matrix or a correlation matrix.  Variables in x and y are specified by location.  If x and y are not specified, then the effect of all variables are partialled from all the other correlations.  May also be done using formula input which is more convenient when comparing results to regression models
}
\usage{
partial.r(data, x, y, use="pairwise",method="pearson")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{A data or correlation matrix}
  \item{x}{The variable names or locations associated with the X set (or formula input) }
  \item{y}{The variable  names or locations associated with the Y set to be partialled from the X set}
  \item{use}{How should we treat missing data? The default is pairwise complete.}
  \item{method}{Which method of correlation should we use, the default is pearson.}
  
}
\details{There are two ways to use \code{\link{partial.r}}.  One is to find the complete partial correlation matrix (that is, partial all the other variables out of each variable).  This may be done by simply specifying the raw data or correlation matrix.  (In the case of raw data, correlations will be found according to use and method.)  In this case, just specify the data matrix. 

 Alternatively, if we think of the data as an X matrix and Y matrix, then (D = X + Y) with correlations R.  Then the partial correlations of the X predictors with the Y variables partialled out are just the last column of R^(-1). See the \code{\link{Tal.Or}} example below.
 
 The second usage is to partial a set of variables(y) out of another set (x). It is sometimes convenient to partial the effect of a number of variables (e.g., sex, age, education) out of the correlations of another set of variables.  This could be done laboriously by finding the residuals of various multiple correlations, and then correlating these residuals.  The matrix algebra alternative is to do it directly. 
To find the confidence intervals and "significance" of the correlations, use the \code{\link{corr.p}} function with n = n - s where s is the number of covariates. 

 
A perhaps easier format is to use formula input compatible with that used in \code{\link{setCor}}. If using formula input,we specify X and Y with the partialled variables specified by subtraction.    That is X ~ Y - z,
This is useful in the case of multiple regression using  which uses this notation.



Following a thoughtful request from Fransisco Wilheim, I just find the correlations of the variables specified in the call (previously I  had found the entire correlation matrix, which is a waste of time and breaks if some variables are non-numeric).)

In the case of non-positive definite matrices, find the Pinv (pseudo inverse) of the matrix.
   
}
\value{The matrix of partial correlations.}
\references{ Revelle, W. (in prep) An introduction to psychometric theory with applications in R. To be published by Springer.  (working draft available at  \url{https://personality-project.org/r/book/} 
}
\author{ William Revelle }

\seealso{ \code{\link{setCor}} for a similar application for regression. \code{\link{lowerMat}} to neatly show a correlation matrix, and \code{\link{corr.p}} to find the confidence intervals of a correlation. }

\examples{
jen <- make.hierarchical()    #make up a correlation matrix 
lowerMat(jen[1:5,1:5])
par.r <- partial.r(jen,c(1,3,5),c(2,4))
lowerMat(par.r)
#or
R <- jen[1:5,1:5]
par.r <- partial.r(R, y = cs(V2,V4))
lowerMat(par.r)
cp <- corr.p(par.r,n=98)  #assumes the jen data based upon n =100.
print(cp,short=FALSE)  #show the confidence intervals as well
 #partial all from all correlations.
lowerMat(partial.r(jen)) 


#Consider the Tal.Or data set.
lowerCor(Tal.Or)
#partial gender and age from these relations (they hardly change)
partial.r(Tal.Or,1:4,cs(gender,age))
#find the partial correlations between the first three variables and the DV (reaction)
round(partial.r(Tal.Or[1:4])[4,1:3],2) #The partial correlations with the criterion

#Consider the eminence data set from Del Giudice.
if(require("psychTools")) {
data(eminence)
partial.r(reputation ~ works + citations - birth.year, data=eminence)
}


}
\keyword{multivariate}

