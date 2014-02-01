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
\alias{best.items}
\alias{lookup}
\alias{fa.lookup} 
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
test.all(p) 
best.items(x,c1=1,cut=.3, abs=TRUE, contents=NULL,digits=2) 
lookup(x,y,c1=NULL)
fa.lookup(f,dictionary,digits=2)
}
\arguments{
  \item{R}{A rectangular matrix or data frame (probably a correlation matrix)}
  \item{x}{A data matrix or data frame or a vector depending upon the function.}
  \item{y}{A data matrix or data frame or a vector}
  \item{digits}{round to digits}
  \item{use}{Should pairwise deletion be done, or one of the other options to cor}
  \item{method}{"pearson", "kendall", "spearman"} 
  \item{value}{the current value of some looping variable}
  \item{max}{The maximum value the loop will achieve}
  \item{label}{what function is looping}
  \item{f}{The object returned from either a factor analysis (fa) or a principal components analysis (principal) }
  \item{flip}{The factor or components to be reversed keyed (by factor number)}
  \item{correct}{Correct for the maximum possible information in this item}
  \item{base}{What is the base for the log function (default=2, e implies base = exp(1))}
  \item{p}{The name of a package to be activated and then have all the examples tested.}
  \item{c1}{The column number or name of the factor loading to sort on for best.items or to find in lookup}
  \item{cut}{Return all values in abs(x[,c1]) > cut.}
  \item{abs}{if TRUE, sort by absolute value in best.items}
  \item{contents}{a two column data.frame of item names and item content.  Allows for long values for item output, even though the column names are short.}
  \item{dictionary}{a data.frame with rownames corresponding to rownames in the f$loadings matrix, and entries (may be multiple columns) of item content.}

  
}
\value{ 

The lower triangle of a matrix, rounded to digits with titles abbreviated to digits + 3 (lowerMat)
or a series of dots (progressBar).

\code{\link{lowerCor}} prints the lower diagonal correlation matrix but returns (invisibly) the full correlation matrix found with the use and method parameters. The default values are for pairwise deletion of variables, and to print to 2 decimal places. 

\code{\link{tableF}} (for tableFast) is a cut down version of table that does no error checking, nor returns pretty output, but is significantly faster than table.  It will just work on two integer vectors.  This is used in polychoric an tetrachoric for about a 50\% speed improvement for large problems.  

\code{\link{shannon}} finds Shannon's H index of information.  Used for estimating the complexity or diversity of the distribution of responses in a vector or matrix. \deqn{H = -\sum{p_i log(p_i) }}

\code{link{test.all}} allows one to test all the examples in specified package.  This allows us to make sure that those examples work when other packages are also loaded. 

\code{\link{best.items}} returns a sorted list of factor loadings or correlations with the labels as provided in the contents field.

\code{\link{lookup}} is a very simple implementation of the match function. 

\code{\link{fa.lookup}} takes a factor/cluster analysis object (or just a keys like matrix), sorts it using \code{\link{fa.sort}} and then matches by row.name to the corresponding dictionary entries.

}

\details{ \code{\link{corFiml}} is a very short wrapper function that uses the power of the lavaan package.  corFiml calls three lavaan functions: getMissingPatterns, getMissingPatternStats and estimate.moments.fiml.  corFiml is included here to allow for Full Information Maximum Likelihood (FIML) correlation and covariance matrices without digging into the lavaan package.  lavaan in turn requires a number of other packages, including boot, MASS, mnormt, pbivnorm, and quadprog.  

\code{\link{lowerCor}} prints out the lower off diagonal matrix rounded to digits with column names abbreviated to digits + 3 characters, but also returns the full and unrounded matrix.  By default, it uses pairwise deletion of variables.  It in turn calls

\code{\link{lowerMat}} which does the pretty printing.  

\code{\link{best.items}} and \code{\link{lookup}} are simple helper functions to summarize correlation matrices or factor loading matrices.  \code{\link{best.items}} will sort the specified column (c1) of x on the basis of the (absolute) value of the column.  The return as a default is just the rowname of the variable with those absolute values > cut.   If there is a dictionary of item content and item names, then include the contents as a two column matrix with item name and item content.

\code{\link{lookup}} is used by \code{\link{best.items}} and will find values in c1 of y that match those in x.  It returns those rows of y of that match x.
Suppose that you have a "dictionary" of the many variables in a study but you want to consider a small subset of them in a data set x.  Then, you can find the entries in the dictionary corresponding to x by lookup(rownames(x),y)  If the column is not specified, then it will match by rownames(y). 

\code{\link{fa.lookup}} is used when examining the output of a factor analysis and one wants the corresponding variable names and contents. The returned object may then be printed in LaTex by using the \code{\link{df2latex}} function with the char option set to TRUE.

Similarly, given a correlation matrix, r, of the x variables, if you want to find the items that most correlate with another item or scale, and then show the contents of that item from the dictionary, best.items(r,c1=column number or name of x, contents = y)


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
best.items(lb,3,cut=.1)
#to make this a latex table 
#df2latex(best.items(lb,2,cut=.2))
#
data(bfi.dictionary)
f2 <- fa(bfi[1:10],2)
fa.lookup(f2,bfi.dictionary)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{multivariate}

