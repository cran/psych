\name{manhattan}
\alias{manhattan}
\title{"Manhattan" plots of correlations with a set of criteria.}
\description{
A useful way of showing the strength of many correlations with a particular criterion is the Manhattan plot.  This is just a plot of correlations ordered by some keying variable.  Useful to understand the basis of items used in \code{\link{bestScales}}.
}
\usage{
manhattan(x, criteria = NULL, keys = NULL, abs = TRUE, ylab = NULL, labels = NULL, 
log.p = FALSE,ci=.05, pch = 21, main = "Manhattan Plot of", adjust="holm",
ylim = NULL, ...)
}
\arguments{
\item{x}{A matrix or data.frame of items or a correlation matrix.}
\item{criteria}{What column names should be predicted.  If a separate file, what are the variables to predict.}
\item{keys}{a keys.list similar to that used in \code{\link{scoreItems}}  }
\item{abs}{Should we show the absolute value of the correlations.}
\item{ylab}{If NULL, will label as either correlations or log (10) of correlations}
\item{labels}{if NULL, will use the names of the keys}
\item{log.p}{Should we show the correlations (log.p = FALSE) or the log of the probabilities of the correlations (TRUE)}
\item{ci}{The probability for the upper and lower confidence intervals -- bonferroni adjusted}
\item{pch}{The default plot chararcter is a filled circle}
\item{main}{The title for each criterion}
\item{adjust}{Which adjustment for multiple correlations should be applied ("holm", "bonferroni", "none")}
\item{ylim}{If NULL will be the min and max of the data}
\item{\dots}{Other graphic parameters}

}
\details{When exploring the correlations of many items with a few criteria, it is useful to form scales from the most correlated items (see \code{\link{bestScales}}.  To get a feeling of the distribution of items across various measures, we can display their correlations (or the log of the probabilities) grouped by some set of scale keys.
}
\value{
The corelations or the log p values are returned (invisibily)
}

\author{
William Revelle
}


\seealso{ \code{\link{bestScales}}, \code{\link{error.dots}}
}
\examples{
op <- par(mfrow=(c(2,3))) #we want to compare two different sets of plots
manhattan(bfi[1:25],bfi[26:28],labels=colnames(bfi)[1:25])
manhattan(bfi[1:25],bfi[26:28],log.p=TRUE)

#Do it again, but now show items by the keys.list
bfi.keys <-
  list(agree=c("-A1","A2","A3","A4","A5"),conscientious=c("C1","C2","C3","-C4","-C5"),
 extraversion=c("-E1","-E2","E3","E4","E5"),neuroticism=c("N1","N2","N3","N4","N5"),
 openness = c("O1","-O2","O3","O4","-O5"))
 manhattan(bfi[1:25],bfi[26:28],keys=bfi.keys)
manhattan(bfi[1:25],bfi[26:28],keys=bfi.keys,log.p=TRUE) 
par <- op



}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{multivariate }
\keyword{hplot }