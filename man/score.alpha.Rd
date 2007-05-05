\name{score.alpha}
\alias{score.alpha}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Score scales and find Cronbach's alpha as well as associated statistics}
\description{Given a matrix or data.frame of k keys for m items (-1, 0, 1), and a matrix or data.frame of items scores for m items and n people, find the sum scores or average scores for each person and each scale.  In addition, report Cronbach's alpha, the average r, the scale intercorrelations, and the item by scale correlations.
}
\usage{
score.alpha(keys, items, labels = NULL, totals=TRUE,digits = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{keys}{ A matrix or dataframe of -1, 0, or 1 weights for each item on each scale }
  \item{items}{Data frame or matrix of raw item scores }
  \item{labels}{column names for the resulting scales}
  \item{totals}{Find sum scores (default) or average score}
  \item{digits}{Number of digits for answer (default =2) }
}
\details{The process of finding sum or average scores for a set of scales given a larger set of items is a typical problem in psychometric research.  Although the structure of scales can be determined from the item intercorrelations, to find scale means, variances, and do further analyses, it is typical to find the sum or the average scale score.  

Various estimates of scale reliability include ``Cronbach's alpha", and the average interitem correlation.  For k = number of items in a scale, and av.r = average correlation between items in the scale, alpha = k * av.r/(1+ (k-1)*av.r).  Thus, alpha is an increasing function of test length as well as the test homeogeneity.  

Alpha is a poor estimate of the general factor saturation of a test (see Zinbarg et al., 2005) for it can seriously overestimate the size of a general factor, and a better but not perfect estimate of total test reliability because it underestimates total reliability. None the less, it is a useful statistic to report. 

}
\value{
  \item{scores }{Sum or average scores for each subject on the k scales}
  \item{alpha }{Cronbach's coefficient alpha.  A simple (but non-optimal) measure of the internal consistency of a test. See also beta and omega.}
  \item{av.r}{The average correlation within a scale, also known as alpha 1 is a useful index of the internal consistency of a domain.}
  \item{n.items}{Number of items on each scale}
  \item{cor}{The intercorrelation of all the scales}
  \item{item.cor}{The correlation of each item with each scale.  Because this is not corrected for item overlap, it will overestimate the amount that an item correlates with the other items in a scale.}
}
\references{An introduction to psychometric theory with applications in R (in preparation).  \url{http://personality-project.org/r/book}}
\author{ William Revelle }

\seealso{ \code{\link{alpha.scale}}, \code{\link{correct.cor}}, \code{\link{alpha.scale}}, \code{\link{cluster.loadings}}, \code{\link{omega}} }
\examples{

y <- attitude     #from the datasets package
keys <- matrix(c(rep(1,7),rep(1,4),rep(0,7),rep(-1,3)),ncol=3)
labels <- c("first","second","third")
x <- score.alpha(keys,y,labels)

## The function is currently defined as
function (keys,items,labels=NULL,totals=TRUE,digits=2) {
    keys <- as.matrix(keys)   #just in case they were not matrices to start with
    items <- as.matrix(items)
    scores <- items \%*\%  keys  #this actually does all the work
    if (length(labels)>0) {colnames(scores) <- labels} #add labels
    abskeys <- abs(keys)
    item.var <- diag(var(items))  #find the item variances
    cov.scales  <- cov(scores)    #and total scale variance
    var.scales <- diag(cov.scales)
    cor.scales <- cor(scores)     #could do this as matrix operation, but why bother
    sum.item.var <- item.var \%*\% abskeys  
    num.item <- diag(t(abskeys) \%*\% abskeys) #how many items in each scale
    alpha.scale <- (var.scales - sum.item.var)*num.item/((num.item-1)*var.scales)
   if (length(labels)>0) {colnames(alpha.scale) <- labels}
   av.r <- alpha.scale/(num.item - alpha.scale*(num.item-1))  #alpha 1 = average r
   item.cor <- cor(items,scores)
    results <- list(scores=scores,alpha=round(alpha.scale,digits), av.r=round(av.r,digits), n.items = num.item, cor = round(cor.scales,digits), item.cor = round(item.cor,digits)) 
 }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{multivariate }
\keyword{models}
