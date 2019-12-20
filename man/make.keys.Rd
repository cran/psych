\name{make.keys}
\alias{make.keys}
\alias{keys2list}
\alias{selectFromKeys}
\title{ Create a keys matrix for use by score.items or cluster.cor}
\description{ When scoring items by forming composite scales either from the raw data using \code{\link{scoreItems}} or from the correlation matrix using \code{\link{cluster.cor}}, it used to be  necessary to create a keys matrix. This is no longer necessary as most of the scoring functions will directly use a keys list. \code{\link{make.keys}}  is just a short cut for creating a keys matrix.  The keys matrix is a nvar x nscales matrix of -1,0, 1 and defines the membership for each scale. Items can be specified by location or by name.
}
\usage{
make.keys(nvars, keys.list, item.labels = NULL, key.labels = NULL)
keys2list(keys,sign=TRUE)
selectFromKeys(keys.list)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nvars}{Number of variables items to be scored, or the name of the data.frame/matrix to be scored}
  \item{keys.list}{ A list of the scoring keys, one element for each scale}
  \item{item.labels}{ Typically, just the colnames of the items data matrix. }
  \item{key.labels}{ Labels for the scales can be specified here, or in the key.list }
  \item{keys}{A keys matrix returned from make.keys}
  \item{sign}{if TRUE, prefix negatively keyed items with - (e.g., ``-E2")}
}
\details{The easiest way to prepare keys for \code{\link{scoreItems}}, \code{\link{scoreOverlap}}, \code{\link{scoreIrt.1pl}}, or \code{\link{scoreIrt.2pl}}  is to specify a keys.list.  This is just a list specifying the name of the scales to be scores and the direction of the items to be used.

In earlier versions (prior to 1.6.9) keys were formed as a matrix of -1, 0, and 1s for all the items using make.keys.  This is no longer necessary, but make.keys is kept for compatibility with earlier versions.


There are three ways to create keys for the \code{\link{scoreItems}}, \code{\link{scoreOverlap}}, \code{\link{scoreIrt.1pl}}, or \code{\link{scoreIrt.2pl}} functions. One is to laboriously do it in a spreadsheet and then copy them into R.  The other is to just specify them by item number in a list. \code{\link{make.keys}} allows one to specify items by name or by location or a mixture of both.

\code{\link{keys2list}} reverses the \code{\link{make.keys}} process and returns a list of scoring keys with the item names for each  item to be keyed.  If sign=FALSE, this is just a list of the items to be scored. (Useful for \code{\link{scoreIrt.2pl}}

\code{\link{selectFromKeys}} will strip the signs from a keys.list and create a vector of item names (deleting duplicates) associated with those keys.  This is useful if using a keys.list to define scales and then just selecting those items that are in subset of the keys.list.  This is now done in the scoring functions in the interest of speed. 

Since these scoring functions \code{\link{scoreItems}}, \code{\link{scoreOverlap}}, \code{\link{scoreIrt.1pl}}, or \code{\link{scoreIrt.2pl}} can now (> version 1.6.9)  just take a keys.list as input, make.keys is not as important, but is kept for documentation purposes.

To address items by name it is necessary to specify item names, either by using the item.labels value, or by putting the name of the data file or the colnames of the data file to be scored into the first (nvars) position.

If specifying by number (location), then nvars is the total number of items in the object to be scored, not just the number of items used.

See the examples for the various options.

Note that make.keys was revised in Sept, 2013 to allow for keying by name.

It is also possible to do several make.keys operations and then combine them using \code{\link{superMatrix}}.


}
\value{
  
  \item{keys }{a nvars x nkeys matrix of -1, 0, or 1s describing how to score each scale. nkeys is the length of the keys.list}
 
  
}
\seealso{ \code{\link{scoreItems}}, \code{\link{scoreOverlap}}, \code{\link{cluster.cor}} \code{\link{superMatrix}} }
\examples{
data(attitude)  #specify the items by location
 key.list <- list(all=c(1,2,3,4,-5,6,7),
                  first=c(1,2,3),
                  last=c(4,5,6,7))
 keys <- make.keys(7,key.list,item.labels = colnames(attitude))
 keys
 #now, undo this 
new.keys.list <- keys2list(keys)  #note, these are now given the variable names

select <- selectFromKeys(key.list)

 
 #scores <- score.items(keys,attitude)
 #scores
 
# data(psychTools::bfi)
 #first create the keys by location (the conventional way)
 keys.list <- list(agree=c(-1,2:5),conscientious=c(6:8,-9,-10),
 extraversion=c(-11,-12,13:15),neuroticism=c(16:20),openness = c(21,-22,23,24,-25))   
 keys <- make.keys(25,keys.list,item.labels=colnames(psychTools::bfi)[1:25])
 new.keys.list <- keys2list(keys)  #these will be in the form of variable names
 
 #alternatively, create by a mixture of names and locations 
 keys.list <- list(agree=c("-A1","A2","A3","A4","A5"),
conscientious=c("C1","C2","C2","-C4","-C5"),extraversion=c("-E1","-E2","E3","E4","E5"),
neuroticism=c(16:20),openness = c(21,-22,23,24,-25)) 
keys <- make.keys(psychTools::bfi, keys.list) #specify the data file to be scored (bfi)
#or
keys <- make.keys(colnames(psychTools::bfi),keys.list) #specify the names of the variables 
#to be used
#or
#specify the number of variables to be scored and their names in all cases
keys <- make.keys(28,keys.list,colnames(psychTools::bfi)) 


 scores <- scoreItems(keys,psychTools::bfi)
 summary(scores)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ multivariate }
\keyword{models}