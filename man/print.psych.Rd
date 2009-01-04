\name{print.psych}
\alias{print.psych}
\alias{summary.psych}


\title{ Print and summary functions for the psych class }
\description{Give limited output (print) or somewhat more detailed (summary) for the  \code{\link{factor.pa}}, \code{\link{omega}}, \code{\link{ICLUST}}, \code{\link{score.items}}, \code{\link{cluster.cor}}, or \code{\link{cluster.loadings}}   functions.
In addition, will supply the factor correlations for output from a promax rotation applied to a factanal output.
}
\usage{
print.psych(x,digits=2,all=FALSE,cutoff=NULL,sort=FALSE,...)
summary.psych(object,digits=2,items=FALSE,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ Output from a psych function (e.g., factor.pa, omega,ICLUST, score.items, cluster.cor}
  \item{object}{Output from a psych function}
  \item{items}{items=TRUE (default) does not print the item whole correlations}
  \item{digits}{Number of digits to use in printing}
  \item{all}{if all=TRUE, then the object is declassed and all output from the function is printed}
  \item{cutoff}{Cluster loadings < cutoff will not be printed.  For factor.pa, cutoff defaults to .3, for omega to .2.}
  \item{sort}{Cluster loadings are in sorted order}
  \item{...}{More options to pass to summary and print}
}
\details{Most of the psych functions produce too much output.  print.psych and summary.psych use generic methods for printing just the highlights.

To get complete output, unclass(theobject) and then print it.

As an added feature, if the promax function is applied to a factanal loadings matrix, the normal output just provides the rotation matrix.  print.psych will provide the factor correlations. (Following a suggestion by John Fox and Uli Keller to the R-help list). 

}
\value{ Various psych functions produce copious output.  This is a way to summarize the most important parts of the output of the score.items, cluster.scores, and ICLUST functions.  See those ( \code{\link{score.items}}, \code{\link{cluster.cor}}, \code{\link{cluster.loadings}}, or \code{\link{ICLUST}}) for details on what is produced.
  
}

\author{ William Revelle}
\note{ See \code{\link{score.items}}, \code{\link{cluster.cor}}, \code{\link{cluster.loadings}}, or \code{\link{ICLUST}}for details on what is printed.


}

\examples{
data(bfi)
 keys.list <- list(agree=c(-1,2:5),conscientious=c(6:8,-9,-10),extraversion=c(-11,-12,13:15),neuroticism=c(16:20),openness = c(21,-22,23,24,-25))
 keys <- make.keys(25,keys.list,item.labels=colnames(bfi))
 scores <- score.items(keys,bfi,short=TRUE)
 scores
 summary(scores)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ multivariate }