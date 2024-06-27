\name{fa.lookup}
\alias{lookup}
\alias{lookupItems}
\alias{fa.lookup}
\alias{item.lookup}
\alias{itemSort}
\alias{keys.lookup}
\alias{lookupFromKeys}
\alias{lmCorLookup}
\title{A set of functions for factorial and empirical scale construction} 

\description{ When constructing scales through rational, factorial, or empirical means, it is useful to examine the content of the items that relate most highly to each other (e.g., the factor loadings of \code{\link{fa.lookup}} of a set of items), or to some specific set of criteria  (e.g., \code{\link{bestScales}}). Given a dictionary of item content, these routines will sort by factor loading, item means,  or criteria correlations and display the item content. 
}
\usage{
lookup(x,y,criteria=NULL)
lookupItems(content=NULL,dictionary=NULL,search=c("Item","Content","item"))
fa.lookup(f,dictionary=NULL,digits=2,cut=.0,n=NULL,sort=TRUE)
item.lookup(f,m, dictionary,cut=.3, digits = 2) 
itemSort(m, dictionary, ascending=TRUE, digits = 2)
keys.lookup(keys.list,dictionary)
lookupFromKeys(keys.list,dictionary,n=20,cors=NULL,sort=TRUE,suppress.names=FALSE,
      digits=2)
lmCorLookup(x,dictionary=NULL,cut=0,digits=2,p=.05) 
}
\arguments{
\item{x}{A data matrix or data frame depending upon the function.}
\item{y}{A data matrix or data frame or a vector}
\item{criteria}{Which variables (by name or location) should be the empirical target for bestScales and bestItems.  May be a separate object. }
\item{f}{The object returned from either a factor analysis (fa) or a principal components analysis (principal) }
\item{content}{The word(s) to search for from a dictionary}
\item{keys.list}{A list of scoring keys suitable to use for make.keys}
\item{cut}{Return all values in abs(x[,c1]) > cut.}
\item{n}{Return the n best items per factor (as long as they have their highest loading on that factor)}
\item{cors}{If provided (e.g. from scoreItems) will be added to the lookupFromKeys output}
\item{dictionary}{a data.frame with rownames corresponding to rownames in the f$loadings matrix or colnames of the data matrix or correlation matrix, and entries (may be multiple columns) of item content. See Notes for how to construct a dictionary.}
\item{search}{Column names of dictionary to search, defaults to "Item" or "Content" (dictionaries have different labels for this column), can search any column specified by search.}
\item{m}{A data frame of item means}
\item{digits}{round to digits}
\item{sort}{Should the factors be sorted first?}
\item{suppress.names}{In lookupFromKeys, should we suppress the column labels}
\item{p}{Show lmCor regressions with probability < p}
\item{ascending}{order to sort the means in itemSort -- see dfOrder}
  
}
\details{

\code{\link{fa.lookup}} and \code{\link{lookup}} are simple helper functions to summarize correlation matrices or factor loading matrices.  \code{\link{bestItems}} will sort the specified column (criteria) of x on the basis of the (absolute) value of the column.  The return as a default is just the rowname of the variable with those absolute values > cut.   If there is a dictionary of item content and item names, then include the contents as a two column (or more) matrix with rownames corresponding to the item name and then as many fields as desired for item content. (See the example dictionary \code{\link[psychTools]{bfi.dictionary}}).

\code{\link{lookup}} is used by \code{\link{bestItems}} and will find values in c1 of y that match those in x.  It returns those rows of y of that match x.
Suppose that you have a "dictionary" of the many variables in a study but you want to consider a small subset of them in a data set x.  Then, you can find the entries in the dictionary corresponding to x by lookup(rownames(x),y)  If the column is not specified, then it will match by rownames(y). 

\code{\link{fa.lookup}} is used when examining the output of a factor analysis and one wants the corresponding variable names and contents. The returned object may then be printed in LaTex by using the \code{\link[psychTools]{df2latex}} function with the char option set to TRUE.

\code{\link{fa.lookup}} will work with output from \code{\link{fa}}, \code{\link{pca}} or \code{\link{omega}}.  For omega output, the items are sorted by the non-general factor loadings.

Similarly, given a correlation matrix, r, of the x variables, if you want to find the items that most correlate with another item or scale, and then show the contents of that item from the dictionary, bestItems(r,c1=column number or name of x, contents = y)

\code{\link{item.lookup}} combines the output from a factor analysis \code{\link{fa}} with simple descriptive statistics (a data frame of means) with a dictionary.  Items are grouped by factor loadings > cut, and then sorted by item mean.  This allows a better understanding of how a scale works, in terms of the meaning of the item endorsements.  Note the means must be a one column matrix (with row names), not a vector (without rownames.) 

\code{\link{itemSort}} Combine item means and item content and then sort them by the item means.

\code{\link{lookupItems}} searches a dictionary for all items that have a certain content.  The rownames of the returned object are the item numbers which can then be used in other functions to find statistics (e.g. omega) of a scale with those items.   If an scales by items correlation matrix is given, then the item correlation with that scale are also shown. 
}
\value{ 
\code{\link{bestItems}} returns a sorted list of factor loadings or correlations with the labels as provided in the dictionary.

\code{\link{lookup}} is a very simple implementation of the match function. 

\code{\link{fa.lookup}} takes a factor/cluster analysis object (or just a keys like matrix), sorts it using \code{\link{fa.sort}} and then matches by row.name to the corresponding dictionary entries.
}

\references{
Revelle, W.  (in preparation) An introduction to psychometric theory with applications in {R}. Springer.  (Available online at \url{https://personality-project.org/r/book/}). 

}
\author{William Revelle}
\note{Although empirical scale construction is appealing, it has the basic problem of capitalizing on chance.  Thus, be careful of over interpreting the results unless working with large samples.  Iteration and bootstrapping aggregation (bagging) gives information on the stability of the solutions. See \code{\link{bestScales}}

To create a dictionary, create an object with row names as the item numbers, and the columns as the item content.  See the bfi.dictionary in the psychTools package as an example.  The bfi.dictionary was constructed from a spreadsheet with multiple columns, the first of which was the column names of the bfi. See the first (not run) example. 

}

\seealso{
\code{\link{fa}}, \code{\link{iclust}},\code{\link{principal}}, \code{\link{bestScales}} and \code{\link{bestItems}}
}
\examples{
#Tne following shows how to create a dictionary
#first, copy the spreadsheet to the clipboard
#the spreadsheet should have multiple columns
#col 1         col 2     col 3
#item	       content   label 
#A1           Am indifferent to the feelings of others.     (q_146)
#A2 	      Inquire about others' well-being. (q_1162
#...

# bfi.dictionary <- read.clipboard.tab()  #read from the clipboard
# rownames(bfi.dictionary) <- bfi.dictionary[1] #the first column had the names
# bfi.dictionary <- bfi.dictionary[-1]  #these are redundant, drop them

f5 <- fa(bfi,5)
m <- colMeans(bfi,na.rm=TRUE)
item.lookup(f5,m,dictionary=bfi.dictionary[2,drop=FALSE])
 #just show the item content, not the source of the items
fa.lookup(f5,dictionary=bfi.dictionary[2]) 

#just show the means and the items
#use the m vector we found above
itemSort(as.matrix(m),dictionary=bfi.dictionary[,2:3,drop=FALSE])


#show how to use lookupFromKeys 
bfi.keys <- 
list(agree=c("-A1","A2","A3","A4","A5"),conscientiousness=c("C1","C2","C3","-C4","-C5"),
extraversion=c("-E1","-E2","E3","E4","E5"),neuroticism=c("N1","N2","N3","N4","N5"),
openness = c("O1","-O2","O3","O4","-O5")) 
bfi.over <- scoreOverlap(bfi.keys,bfi) #returns the corrected for overlap values
lookupFromKeys(bfi.keys,bfi.dictionary,n=5, cors=bfi.over$item.cor)
 #show the keying information
 if(require(psychTools)){
lookupItems("life",psychTools::spi.dictionary) #find those items with "life" in the item
}

}

\keyword{ models }
\keyword{multivariate }

