\name{psych.misc}
\alias{psych.misc}
\alias{misc}
\alias{lowerCor}
\alias{lowerMat}
\alias{progressBar}
\alias{reflect}
\alias{shannon}
\alias{test.all}
\title{Miscellaneous helper functions for the psych package}
\description{This is a set of minor, if not trivial, helper functions.
lowerCor finds the correlation of x variables and then prints them using 
lowerMat which is a trivial, but useful, function to round off and print the lower triangle of a matrix.
reflect reflects the output of a factor analysis or principal components analysis so that one or more factors is reflected. (Requested by Alexander Weiss.)
progressBar prints out ...  as a calling routine (e.g., \code{\link{tetrachoric}}) works through a tedious calculation.  shannon finds the Shannon index (H) of diversity or of information. test.all tests all the examples in a package.
}
\usage{
psych.misc()
lowerCor(x,digits=2,use="pairwise",method="pearson")
lowerMat(R, digits = 2)
reflect(f,flip=NULL)
progressBar(value,max,label=NULL) 
shannon(x,correct=FALSE,base=2)
test.all(p) 
}
\arguments{
  \item{R}{A rectangular matrix or data frame (probably a correlation matrix)}
  \item{x}{A data matrix or data frame}
  \item{digits}{round to digits}
  \item{use}{Should pairwise deletion be done, or one of the other options to cor}
  \item{method}{"pearson", "kendall", "spearman"} 
  \item{value}{the current value of some looping variable}
  \item{max}{The maximum value the loop will achieve}
  \item{label}{what function is looping}
  \item{f}{The object returned from either a factor analysis (fa) or a principal components analysis (principal) }
  \item{flip}{The factor or components to be reversed keyed (by factor number)}
  \item{correct}{Correct for the maximum possible information in this item}
  \item{base}{What is the base for the log function (default=2, e implies base = exp(1) ) }
  \item{p}{The name of a package to be activated and then have all the examples tested.}
}
\value{ 

The lower triangle of a matrix, rounded to digits with titles abbreviated to digits + 3 (lowerMat)
or a series of dots (progressBar).

lowerCor prints the lower diagonal correlation matrix but returns (invisibly) the full correlation matrix found with the use and method parameters. The default values are for pairwise deletion of variables, and to print to 2 decimal places. 


\code{\link{shannon}} finds Shannon's H index of information.  Used for estimating the complexity or diversity of the distribution of responses in a vector or matrix. \deqn{H = -\sum{p_i log(p_i) }}

\code{link{test.all}} allows one to test all the examples in specified package.  This allows us to make sure that those examples work when other packages are also loaded.  

}

\details{ corFiml is a very short wrapper function that uses the power of the lavaan package.  corFiml calls three lavaan functions: getMissingPatterns, getMissingPatternStats and estimate.moments.fiml.  corFiml is included here to allow for Full Information Maximum Likelihood (FIML) correlation and covariance matrices without digging into the lavaan package.  lavaan in turn requires a number of other packages, including boot, MASS, mnormt, pbivnorm, and quadprog.  

lowerCor prints out the lower off diagonal matrix rounded to digits with column names abbreviated to digits + 3 characters, but also returns the full and unrounded matrix.  By default, it uses pairwise deletion of variables.  It in turn calls

lowerMat which does the pretty printing.  

}

\seealso{\code{\link{corr.test}} to find correlations, count the pairwise occurrences, and to give significance tests for each correlation.  \code{\link{r.test}} for a number of tests of correlations, including tests of the difference between correlations.  \code{\link{lowerUpper}} will display the differences between two matrices.}
\examples{
lowerMat(Thurstone)
lowerCor(bfi[1:10])
#fiml <- corFiml(bfi[1:10])     #FIML correlations require lavaan package
#lowerMat(fiml)  #to get pretty output
f3 <- fa(Thurstone,3)
f3r <- reflect(f3,2)  #reflect the second factor
#find the complexity of the response patterns of the iqitems.
round(shannon(iqitems),2) 
#test.all('BinNor')  #Does the BinNor package work when we are using other packages
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{multivariate}

