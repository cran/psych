\name{biplot.psych}
\alias{biplot.psych}
\title{Draw biplots of factor or component scores by factor or component loadings}
\description{Extends the biplot function to the output of \code{\link{fa}}, \code{\link{fa.poly}}  or \code{\link{principal}}. Will plot factor scores and factor loadings in the same graph.  If the number of factors > 2, then all pairs of factors are plotted. Factor score histograms are plotted on the diagonal. The input is the resulting object from \code{\link{fa}}, \code{\link{principal}}, or \}{code\{link{fa.poly}} with the scores=TRUE option. Points may be colored according to other criteria.
}
\usage{
\method{biplot}{psych}(x, labels=NULL,cex=c(.75,1),main="Biplot from fa",
hist.col="cyan",xlim.s=c(-3,3),ylim.s=c(-3,3),xlim.f=c(-1,1),ylim.f=c(-1,1),
maxpoints=100,adjust=1.2,col,pos, arrow.len = 0.1,pch=16,choose=NULL,
cuts=1,cutl=.0,group=NULL,smoother=FALSE,vars=TRUE,...)
}

\arguments{
  \item{x}{The output from \code{\link{fa}}, \code{\link{fa.poly}} or \code{\link{principal}}  with the scores=TRUE option}
  \item{labels}{if NULL, draw the points with the plot character (pch) specified. 
   To identify the data points, specify labels= 1:n  where n is the number of 
   observations, or labels =rownames(data) where data was the data set analyzed
    by the factor analysis.}
  \item{cex}{A vector of plot sizes of the data labels and of the factor labels}
  \item{main}{A main title for a two factor biplot}
  \item{hist.col}{If plotting more than two factors, the color of the histogram of the factor scores}
  \item{xlim.s}{x limits of the scores. Defaults to plus/minus three sigma}
  \item{ylim.s}{y limits of the scores.Defaults to plus/minus three sigma}
  \item{xlim.f}{x limits of the factor loadings.Defaults to plus/minus 1.0}
  \item{ylim.f}{y limits of the factor loadings.Defaults to plus/minus 1.0}
  \item{maxpoints}{When plotting 3 (or more) dimensions, at what size should we switch 
  from plotting "o" to plotting "."}
  \item{adjust}{an adjustment factor in the histogram}
  \item{col}{a vector of colors for the data points and for the factor loading labels}
  \item{pos}{If plotting labels, what position should they be in? 1=below, 2=left,
   3 top, 4 right.
  If missing, then the assumption is that labels should be printed instead of data points.}
  \item{arrow.len}{ the length of the arrow head}
  \item{pch}{The plotting character to use.  pch=16 gives reasonable size dots.
    pch="." gives tiny points.  If adding colors, use pch between 21 and 25. 
     (see examples).}
  \item{choose}{Plot just the specified factors}
  \item{cuts}{Do not label cases with abs(factor scores) < cuts) (Actually,
   the  distance of the x and y scores from 0) }
  \item{cutl}{Do not label variables with communalities in the two space < cutl}
  \item{group}{A vector of a grouping variable for the scores.  Show a different color 
  and symbol for each group.}
  \item{smoother}{If TRUE then do a smooth scatter plot (which shows the density rather than the data points). Only useful for large data sets.}
  \item{vars}{If TRUE, draw arrows for the variables, and plot the scores.  If FALSE, then draw arrows for the scores and plot the variables.}
  \item{\dots}{more options for graphics}
 
}
\details{
Uses the generic biplot function to take the output of a factor analysis \code{\link{fa}}, \code{\link{fa.poly}} or principal components analysis \code{\link{principal}} and plot the factor/component scores along with the factor/component loadings.

This is an extension of the generic biplot function to allow more control over plotting points in a two space and also to plot three or more factors (two at time).  

This will work for objects produced by \code{\link{fa}}, \code{\link{fa.poly}} or \code{\link{principal}}  if they applied to the original data matrix.  If however, one has a correlation matrix (e.g., based upon the output from \code{\link{tetrachoric}} or \code{\link{polychoric}}), and has done either \code{\link{fa}} or \code{\link{principal}} on the correlations, then obviously, we can not do a biplot.  

However, both of those functions produce a weights matrix, which, in combination with the original data can be used to find the scores by using \code{\link{factor.scores}}.  Since biplot.psych is looking for two elements of the x object: x$loadings and x$scores, you can create the appropriate object to plot, or add it to the factor object See the third and 	fourth examples.

In order to just plot the loadings, use \code{\link{fa.plot}}

 
}

\author{William Revelle}


\seealso{\code{\link{fa}}, \code{\link{fa.poly}},  \code{\link{principal}}, \code{\link{fa.plot}}, \code{\link{pairs.panels}}

}
\examples{
#the standard example
data(USArrests)
fa2 <- fa(USArrests,2,scores=TRUE)
biplot(fa2,labels=rownames(USArrests))

# plot the 3 factor solution
#data(bfi)
fa3 <- fa(psychTools::bfi[1:200,1:15],3,scores=TRUE)
biplot(fa3)
#just plot factors 1 and 3 from that solution
biplot(fa3,choose=c(1,3))

#
fa2 <- fa(psychTools::bfi[16:25],2)  #factor analysis
fa2$scores <- fa2$scores[1:100,]  #just take the first 100
#now plot with different colors and shapes for males and females
biplot(fa2,pch=c(24,21)[psychTools::bfi[1:100,"gender"]],
group =psychTools::bfi[1:100,"gender"],
   main="Biplot of Openness and Neuroticism by gender")

## An example from the correlation matrix
r <- cor(psychTools::bfi[1:200,1:10], use="pairwise") #find the correlations
f2 <- fa(r,2) 
#biplot(f2)  #this throws an error    (not run)

#f2 does not have scores, but we can find them
f2$scores <- factor.scores(psychTools::bfi[1:200,1:10],f2)

biplot(f2,main="biplot from correlation matrix and factor scores")


#or create a new object with the scores
 #find the correlations for all subjects
r <- cor(psychTools::bfi[1:10], use="pairwise") 
f2 <- fa(r,2) 
x <- list() 
#find the scores for just the first 200 subjects
x$scores <- factor.scores(psychTools::bfi[1:200,1:10],f2)
x$loadings <- f2$loadings
class(x) <- c('psych','fa')
biplot(x,main="biplot from correlation matrix combined with factor scores")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{multivariate }
\keyword{hplot }% __ONLY ONE__ keyword per line
