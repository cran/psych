\name{biplot.psych}
\alias{biplot.psych}
\title{Draw biplots of factor or component scores by factor or component loadings}
\description{Applies the biplot function to the output of \code{\link{fa}}  or \code{\link{principal}}. Will plot factor scores and factor loadings in the same graph.  If the number of factors > 2, then all pairs of factors are plotted. Factor score histograms are plotted on the diagonal. The input is the resulting object from \code{\link{fa}}  or \code{\link{principal}} with the scores=TRUE option. 
}
\usage{
biplot.psych(x, labels=NULL,cex=c(.75,1),main="Biplot",hist.col="cyan",xlim=c(-3,3),ylim=c(-3,3),...)
}

\arguments{
  \item{x}{The output from \code{\link{fa}}  or \code{\link{principal}}  with the scores=TRUE option}
  \item{labels}{if NULL, draw the points with small o.  To identify the data points, specify labels= 1:n  where n is the number of observations, or labels =rownames(data) where data was the data set analyzed by the factor analysis.}
  \item{cex}{plot size of the data points and of the factor labels}
  \item{main}{A main title for a two factor biplot}
  \item{hist.col}{If plotting more than two factors, the color of the histogram of the factor scores}
  \item{xlim}{Defaults to plus/minus three sigma}
  \item{ylim}{Defaults to plus/minus three sigma}
  \item{\dots}{more options for graphics}
}
\details{
Uses the generic biplot function to take the output of a factor analysis \code{\link{fa}} or principal components analysis \code{\link{principal}} and plot the factor/component scores along with the factor/component loadings.  
}

\author{William Revelle}


\seealso{\code{\link{fa}}, \code{\link{principal}}, \code{\link{factor.plot}}, \code{\link{pairs.panels}}

}
\examples{
data(USArrests)
fa2 <- fa(USArrests,2,scores=TRUE)
biplot(fa2,labels=rownames(USArrests))
data(bfi)
fa3 <- fa(bfi[1:200,1:15],3,scores=TRUE)
biplot(fa3)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{multivariate }
\keyword{hplot }% __ONLY ONE__ keyword per line
