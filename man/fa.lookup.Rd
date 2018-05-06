\name{fa.lookup}
\alias{lookup}
\alias{fa.lookup}
\alias{item.lookup}
\alias{keys.lookup}
\title{A set of functions for factorial and empirical scale construction} 

\description{ When constructing scales through rational, factorial, or empirical means, it is useful to examine the content of the items that relate most highly to each other (e.g., the factor loadings of \code{\link{fa.lookup}} of a set of items) , or to some specific set of criteria  (e.g., \code{\link{bestScales}}). Given a dictionary of item content, these routines will sort by factor loading or criteria correlations and display the item content. 
}
\usage{
lookup(x,y,criteria=NULL)
fa.lookup(f,dictionary,digits=2)
item.lookup(f,m, dictionary,cut=.3, digits = 2) 
keys.lookup(keys.list,dictionary)
}
\arguments{
\item{x}{A data matrix or data frame depending upon the function.}
\item{y}{A data matrix or data frame or a vector}
\item{criteria}{Which variables (by name or location) should be the empirical target for bestScales and bestItems.  May be a separate object. }
\item{f}{The object returned from either a factor analysis (fa) or a principal components analysis (principal) }
\item{keys.list}{A list of scoring keys suitable to use for make.keys}
\item{cut}{Return all values in abs(x[,c1]) > cut.}
\item{dictionary}{a data.frame with rownames corresponding to rownames in the f$loadings matrix or colnames of the data matrix or correlation matrix, and entries (may be multiple columns) of item content.}
\item{m}{A data frame of item means}
\item{digits}{round to digits}

  
}
\details{

\code{\link{fa.lookup}} and \code{\link{lookup}} are simple helper functions to summarize correlation matrices or factor loading matrices.  \code{\link{bestItems}} will sort the specified column (criteria) of x on the basis of the (absolute) value of the column.  The return as a default is just the rowname of the variable with those absolute values > cut.   If there is a dictionary of item content and item names, then include the contents as a two column matrix with rownames corresponding to the item name and then as many fields as desired for item content. (See the example dictionary \code{\link{bfi.dictionary}}).

\code{\link{lookup}} is used by \code{\link{bestItems}} and will find values in c1 of y that match those in x.  It returns those rows of y of that match x.
Suppose that you have a "dictionary" of the many variables in a study but you want to consider a small subset of them in a data set x.  Then, you can find the entries in the dictionary corresponding to x by lookup(rownames(x),y)  If the column is not specified, then it will match by rownames(y). 

\code{\link{fa.lookup}} is used when examining the output of a factor analysis and one wants the corresponding variable names and contents. The returned object may then be printed in LaTex by using the \code{\link{df2latex}} function with the char option set to TRUE.

Similarly, given a correlation matrix, r, of the x variables, if you want to find the items that most correlate with another item or scale, and then show the contents of that item from the dictionary, bestItems(r,c1=column number or name of x, contents = y)

\code{\link{item.lookup}} combines the output from a factor analysis \code{\link{fa}} with simple descriptive statistics (a data frame of means) with a dictionary.  Items are grouped by factor loadings > cut, and then sorted by item mean.  This allows a better understanding of how a scale works, in terms of the meaning of the item endorsements.  
}
\value{ 
\code{\link{bestItems}} returns a sorted list of factor loadings or correlations with the labels as provided in the dictionary.

\code{\link{lookup}} is a very simple implementation of the match function. 

\code{\link{fa.lookup}} takes a factor/cluster analysis object (or just a keys like matrix), sorts it using \code{\link{fa.sort}} and then matches by row.name to the corresponding dictionary entries.
}

\note{

To create a dictionary, create an object with row names as the item numbers, and the columns as the item content.  See the \code{link{bfi.dictionary}} as an example.
}
\references{
Revelle, W.  (in preparation) An introduction to psychometric theory with applications in {R}. Springer.  (Available online at \url{https://personality-project.org/r/book}). 

}
\author{William Revelle}
\note{Although empirical scale construction is appealing, it has the basic problem of capitalizing on chance.  Thus, be careful of over interpreting the results unless working with large samples.  Iteration and bootstrapping aggregation (bagging) gives information on the stability of the solutions.
}

\seealso{
\code{\link{fa}}, \code{\link{iclust}},\code{\link{principal}}, \code{\link{bestScales}} and \code{\link{bestItems}}
}
\examples{

f5 <- fa(bfi,5)
m <- colMeans(bfi,na.rm=TRUE)
item.lookup(f5,m,dictionary=bfi.dictionary[2])
fa.lookup(f5,dictionary=bfi.dictionary[2])  #just show the item content, not the source of the items

}

\keyword{ models }
\keyword{multivariate }

