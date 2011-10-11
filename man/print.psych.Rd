\name{print.psych}
\alias{print.psych}
\alias{summary.psych}


\title{ Print and summary functions for the psych class }
\description{Give limited output (print) or somewhat more detailed (summary) for most of the functions in psych. 
}
\usage{
\method{print}{psych}(x,digits=2,all=FALSE,cut=NULL,sort=FALSE,short=TRUE,...)
\method{summary}{psych}(object,digits=2,items=FALSE,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ Output from a psych function (e.g., factor.pa, omega,ICLUST, score.items, cluster.cor}
  \item{object}{Output from a psych function}
  \item{items}{items=TRUE (default) does not print the item whole correlations}
  \item{digits}{Number of digits to use in printing}
  \item{all}{if all=TRUE, then the object is declassed and all output from the function is printed}
  \item{cut}{Cluster loadings < cut will not be printed.  For the factor analysis functions (fa and factor.pa etc.), cut defaults to 0, for ICLUST  to .3, for omega to .2.}
  \item{sort}{Cluster loadings are in sorted order}
  \item{short}{Controls how much to print} 
  \item{...}{More options to pass to summary and print}
}
\details{Most of the psych functions produce too much output.  print.psych and summary.psych use generic methods for printing just the highlights.  To see what else is available,  ask for the structure of the particular object: (str(theobject).

Alternatively, to get complete output, unclass(theobject) and then print it. This may be done by using the all=TRUE option.

As an added feature, if the promax function is applied to a factanal loadings matrix, the normal output just provides the rotation matrix.  print.psych will provide the factor correlations. (Following a suggestion by John Fox and Uli Keller to the R-help list).  The alternative is to just use the Promax function directly on the factanal object.

}
\value{ Various psych functions produce copious output.  This is a way to summarize the most important parts of the output of the score.items, cluster.scores, and ICLUST functions.  See those ( \code{\link{score.items}}, \code{\link{cluster.cor}}, \code{\link{cluster.loadings}}, or \code{\link{ICLUST}}) for details on what is produced.

  
}

\author{ William Revelle}
\note{ See \code{\link{score.items}}, \code{\link{cluster.cor}}, \code{\link{cluster.loadings}}, or \code{\link{ICLUST}}for details on what is printed.


}

\examples{
data(bfi)
 keys.list <- list(agree=c(-1,2:5),conscientious=c(6:8,-9,-10),extraversion=c(-11,-12,13:15),neuroticism=c(16:20),openness = c(21,-22,23,24,-25))
 keys <- make.keys(25,keys.list,item.labels=colnames(bfi[1:25]))
 scores <- score.items(keys,bfi[1:25])
 scores
 summary(scores)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ multivariate }
