\name{make.keys}
\alias{make.keys}
\title{ Create a keys matrix for use by score.items or cluster.cor}
\description{ When scoring items by forming composite scales either from the raw data using \code{\link{score.items}} or from the correlatio matrix using \code{\link{cluster.cor}}, it is necessary to create a keys matrix.  This is just a short cut for doing so.
}
\usage{
make.keys(nvars, keys.list, key.labels = NULL, item.labels = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nvars}{Number of variables items to be scored}
  \item{keys.list}{ A list of the scoring keys,one element for each scale}
  \item{key.labels}{ Labels for the scales can be specified here, or in the key.list }
  \item{item.labels}{ Typically, just the colnames of the items data matrix. }
}
\details{
}
\value{
  
  \item{keys }{a nvars x nkeys matrix of -1, 0, or 1s describing how to score each scale. nkeys is the length of the keys.list}
 
  
}
\seealso{ \code{\link{score.items}}, \code{\link{cluster.cor}}, ~~~ }
\examples{
data(attitude)
 key.list <- list(all=c(1,2,3,4,-5,6,7),
                  first=c(1,2,3),
                  last=c(4,5,6,7))
 keys <- make.keys(7,key.list,item.labels = colnames(attitude))
 keys
 
 scores <- score.items(keys,attitude,short=TRUE)
 scores
 
 data(bfi)
 keys.list <- list(agree=c(-1,2:5),conscientious=c(6:8,-9,-10),extraversion=c(-11,-12,13:15),neuroticism=c(16:20),openness = c(21,-22,23,24,-25))
 keys <- make.keys(25,keys.list,item.labels=colnames(bfi))
 scores <- score.items(keys,bfi,short=TRUE)
 scores

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ multivariate }
\keyword{models}