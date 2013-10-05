\name{make.keys}
\alias{make.keys}
\title{ Create a keys matrix for use by score.items or cluster.cor}
\description{ When scoring items by forming composite scales either from the raw data using \code{\link{score.items}} or from the correlation matrix using \code{\link{cluster.cor}}, it is necessary to create a keys matrix.  This is just a short cut for doing so.  The keys matrix is a nvar x nscales matrix of -1,0, 1 and defines the membership for each scale. Items can be specified by location or by name.
}
\usage{
make.keys(nvars, keys.list, item.labels = NULL, key.labels = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nvars}{Number of variables items to be scored}
  \item{keys.list}{ A list of the scoring keys,one element for each scale}
  \item{item.labels}{ Typically, just the colnames of the items data matrix. }
  \item{key.labels}{ Labels for the scales can be specified here, or in the key.list }
  
}
\details{There are two ways to create keys for the \code{\link{score.items}} function. One is to laboriously do it in a spreadsheet and then copy them into R.  The other is to just specify them by item number in a list. Make keys allows one to specify items by name or by location or a mixture of both.

To address items by name it is necessary to specify item names, either by using the item.labels value, or by putting the name of the data file or the colnames of the data file to be scored into the first (nvars) position.

If specifying by number, then nvars is the total number of items in the object to be scored, not just the number of items used.

See the examples for the various options.

Note that make.keys was revised in Sept, 2013 to allow for keying by name.
}
\value{
  
  \item{keys }{a nvars x nkeys matrix of -1, 0, or 1s describing how to score each scale. nkeys is the length of the keys.list}
 
  
}
\seealso{ \code{\link{score.items}},  \code{\link{cluster.cor}} }
\examples{
data(attitude)  #specify the items by location
 key.list <- list(all=c(1,2,3,4,-5,6,7),
                  first=c(1,2,3),
                  last=c(4,5,6,7))
 keys <- make.keys(7,key.list,item.labels = colnames(attitude))
 keys
 
 #scores <- score.items(keys,attitude)
 #scores
 
 data(bfi)
 #first create the keys by location (the conventional way)
 keys.list <- list(agree=c(-1,2:5),conscientious=c(6:8,-9,-10),
 extraversion=c(-11,-12,13:15),neuroticism=c(16:20),openness = c(21,-22,23,24,-25))   
 keys <- make.keys(25,keys.list,item.labels=colnames(bfi)[1:25])
 
 #alternatively, create by a mixture of names and locations 
 keys.list <- list(agree=c("-A1","A2","A3","A4","A5"),
conscientious=c("C1","C2","C2","-C4","-C5"),extraversion=c("-E1","-E2","E3","E4","E5"),
neuroticism=c(16:20),openness = c(21,-22,23,24,-25)) 
keys <- make.keys(bfi,keys.list) #specify the data file to be scored (bfi)
#or
keys <- make.keys(colnames(bfi),keys.list) #specify the names of the variables to be used
#or
#specify the number of variables to be scored and their names in all cases
keys <- make.keys(28,keys.list,colnames(bfi)) 


 scores <- score.items(keys,bfi)
 summary(scores)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ multivariate }
\keyword{models}