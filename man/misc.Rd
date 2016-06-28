\name{psych.misc}
\alias{psych.misc}
\alias{misc}
\alias{tableF}
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
progressBar prints out ...  as a calling routine (e.g., \code{\link{tetrachoric}}) works through a tedious calculation.  shannon finds the Shannon index (H) of diversity or of information. test.all tests all the examples in a package.  best.items sorts a factor matrix for absolute values and displays the expanded items names. fa.lookup returns sorted factor analysis output with item labels.
}
\usage{
psych.misc()
lowerCor(x,digits=2,use="pairwise",method="pearson")
lowerMat(R, digits = 2)
tableF(x,y)
reflect(f,flip=NULL)
progressBar(value,max,label=NULL) 
shannon(x,correct=FALSE,base=2)
test.all(pl,package="psych",dependencies 
      = c("Depends", "Imports", "LinkingTo"),find=FALSE,skip=NULL) 

}
\arguments{
  \item{R}{A rectangular matrix or data frame (probably a correlation matrix)}
  \item{x}{A data matrix or data frame or a vector depending upon the function.}
  \item{y}{A data matrix or data frame or a vector}
  \item{f}{The object returned from either a factor analysis (fa) or a principal components analysis (principal) }
  \item{digits}{round to digits}
  \item{use}{Should pairwise deletion be done, or one of the other options to cor}
  \item{method}{"pearson", "kendall", "spearman"} 
  \item{value}{the current value of some looping variable}
  \item{max}{The maximum value the loop will achieve}
  \item{label}{what function is looping}
  \item{flip}{The factor or components to be reversed keyed (by factor number)}
  \item{correct}{Correct for the maximum possible information in this item}
  \item{base}{What is the base for the log function (default=2, e implies base = exp(1))}
  \item{pl}{The name of a package (or list of packages) to be activated and then have all the examples tested.}
  \item{package}{Find the dependencies for this package, e.g., psych}
  \item{dependencies}{Which type of dependency to examine?}
  \item{find}{Look up the dependencies, and then test all of their examples}
  \item{skip}{Do not test these dependencies}
  

}
\value{ 

\code{\link{tableF}} is fast alternative to the table function for creating two way tables of numeric variables.  It does not have any of the elegant checks of the table function and thus is much faster.  Used in the \code{\link{tetrachoric}} and \code{\link{polychoric}} functions to maximize speed.

The lower triangle of a matrix, rounded to digits with titles abbreviated to digits + 3 (lowerMat)
or a series of dots (progressBar).

\code{\link{lowerCor}} prints the lower diagonal correlation matrix but returns (invisibly) the full correlation matrix found with the use and method parameters. The default values are for pairwise deletion of variables, and to print to 2 decimal places. 

\code{\link{tableF}} (for tableFast) is a cut down version of table that does no error checking, nor returns pretty output, but is significantly faster than table.  It will just work on two integer vectors.  This is used in polychoric an tetrachoric for about a 50\% speed improvement for large problems.  

\code{\link{shannon}} finds Shannon's H index of information.  Used for estimating the complexity or diversity of the distribution of responses in a vector or matrix. \deqn{H = -\sum{p_i log(p_i) }}

\code{\link{test.all}} allows one to test all the examples in specified package.  This allows us to make sure that those examples work when other packages (e.g., psych) are also loaded.  This is used when developing revisions to the psych package to make sure the the other packages work.  Some packages will not work and/or crash the system (e.g., DeducerPlugInScaling requires Java and even with Java, crashes when loaded, even if psych is not there!).  Alternatively, if testing a long list of dependencies, you can skip the first part by specifying them by name.

}

\details{

\code{\link{lowerCor}} prints out the lower off diagonal matrix rounded to digits with column names abbreviated to digits + 3 characters, but also returns the full and unrounded matrix.  By default, it uses pairwise deletion of variables.  It in turn calls

\code{\link{lowerMat}} which does the pretty printing.  

It is important to remember to not call \code{\link{lowerCor}} when all you need is \code{\link{lowerMat}}!  


}

\seealso{\code{\link{corr.test}} to find correlations, count the pairwise occurrences, and to give significance tests for each correlation.  \code{\link{r.test}} for a number of tests of correlations, including tests of the difference between correlations.  \code{\link{lowerUpper}} will display the differences between two matrices.}
\examples{
lowerMat(Thurstone)
lb <- lowerCor(bfi[1:10])  #finds and prints the lower correlation matrix, 
  # returns the square matrix.
#fiml <- corFiml(bfi[1:10])     #FIML correlations require lavaan package
#lowerMat(fiml)  #to get pretty output
f3 <- fa(Thurstone,3)
f3r <- reflect(f3,2)  #reflect the second factor
#find the complexity of the response patterns of the iqitems.
round(shannon(iqitems),2) 
#test.all('BinNor')  #Does the BinNor package work when we are using other packages
bestItems(lb,3,cut=.1)
#to make this a latex table 
#df2latex(bestItems(lb,2,cut=.2))
#
data(bfi.dictionary)
f2 <- fa(bfi[1:10],2)
fa.lookup(f2,bfi.dictionary)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{multivariate}

