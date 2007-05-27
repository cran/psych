\name{score.items}
\alias{score.items}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Score item composite scales and find Cronbach's alpha as well as associated statistics }
\description{Given a matrix or data.frame of k keys for m items (-1, 0, 1), and a matrix or data.frame of items scores for m items and n people, find the sum scores or average scores for each person and each scale.  In addition, report Cronbach's alpha, the average r, the scale intercorrelations, and the item by scale correlations.  Replace missing values with the item mean if desired. Will adjust scores for reverse scored items.
}
\usage{
score.items(keys, items, totals = TRUE, ilabels = NULL, missing = TRUE, min = NULL, max = NULL, digits = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{keys}{A matrix or dataframe of -1, 0, or 1 weights for each item on each scale }
  \item{items}{ Matrix or dataframe of raw item scores }
  \item{totals}{ if TRUE (default) find total scores, if FALSE, find average scores }
  \item{ilabels}{ a vector of item labels  }
  \item{missing}{TRUE: Replace missing values with the corresponding item mean. FALSE: do not score that subject  }
  \item{min}{May be specified as minimum item score allowed, else will be calculated from data}
  \item{max}{May be specified as maximum item score allowed, else will be calculated from data }
  \item{digits}{ Number of digits to report }
}
\details{The process of finding sum or average scores for a set of scales given a larger set of items is a typical problem in psychometric research.  Although the structure of scales can be determined from the item intercorrelations, to find scale means, variances, and do further analyses, it is typical to find scores based upon the sum or the average item score.  

Various estimates of scale reliability include ``Cronbach's alpha", and the average interitem correlation.  For k = number of items in a scale, and av.r = average correlation between items in the scale, alpha = k * av.r/(1+ (k-1)*av.r).  Thus, alpha is an increasing function of test length as well as the test homeogeneity.  

Alpha is a poor estimate of the general factor saturation of a test (see Zinbarg et al., 2005) for it can seriously overestimate the size of a general factor, and a better but not perfect estimate of total test reliability because it underestimates total reliability. None the less, it is a useful statistic to report. 

Correlations between scales are attenuated by a lack of reliability.  Correcting correlations for reliability (by dividing by the square roots of the reliabilities of each scale) sometimes help show structure.
}
\value{ 
 \item{scores }{Sum or average scores for each subject on the k scales}
  \item{alpha }{Cronbach's coefficient alpha.  A simple (but non-optimal) measure of the internal consistency of a test. See also beta and omega. Set to 1 for scales of length 1. }
  \item{av.r}{The average correlation within a scale, also known as alpha 1 is a useful index of the internal consistency of a domain. Set to 1 for scales with 1 item.}
  \item{n.items}{Number of items on each scale}
   \item{item.cor}{The correlation of each item with each scale.  Because this is not corrected for item overlap, it will overestimate the amount that an item correlates with the other items in a scale.}
  \item{cor}{The intercorrelation of all the scales}
  \item{corrected}{The correlations of all scales (below the diagonal), alpha on the diagonal, and the unattenuated correlations (above the diagonal)}
 
}
\references{An introduction to psychometric theory with applications in R (in preparation).  \url{http://personality-project.org/r/book}}
\author{ William Revelle }

\seealso{ \code{\link{alpha.scale}}, \code{\link{correct.cor}}, \code{\link{cluster.cor}} , \code{\link{cluster.loadings}}, \code{\link{omega}} }
\examples{
y <- attitude     #from the datasets package
keys <- matrix(c(rep(1,7),rep(1,4),rep(0,7),rep(-1,3)),ncol=3)
colnames(keys) <- c("first","second","third")
x <- score.items(keys,y)
#x    #to see the output 
## The function is currently defined as
function (keys,items,totals=TRUE,ilabels=NULL, missing=TRUE, min=NULL,max=NULL,digits=2) {
    keys <- as.matrix(keys)   #just in case they were not matrices to start with
    items <- as.matrix(items)
    item.means <- colMeans(items,na.rm=TRUE)
    if (is.null(min)) {min <- min(items,na.rm=TRUE)}
    if (is.null(max)) {max <- max(items,na.rm=TRUE)}
    if(missing) { miss.rep <- rowSums(is.na(items))
       item.means <- colMeans(items,na.rm=TRUE)
            miss <- which(is.na(items),arr.ind=TRUE)
           items[miss]<- item.means[miss[,2]]}
    n.keys <- dim(keys)[2]
    n.items <- dim(keys)[1]
    n.subjects <- dim(items)[1]
     scores<- items \%*\%  keys  #this actually does all the work but doesn't handle missing
    
    scores<- items \%*\%  keys  #this actually does all the work
    slabels <- colnames(keys)
    if (is.null(slabels)) {slabels<- paste("S",1:n.keys,sep="")}
    colnames(scores) <- slabels  
    
    abskeys <- abs(keys)
    item.var <- diag(var(items,use="pairwise"))  #find the item variances
    cov.scales  <- cov(scores,use="pairwise")    #and total scale variance
    var.scales <- diag(cov.scales)
    cor.scales <- cor(scores,use="pairwise")     #could do this as matrix operation, but why bother
    sum.item.var <- item.var \%*\% abskeys  
    num.item <- diag(t(abskeys) \%*\% abskeys) #how many items in each scale
    alpha.scale <- (var.scales - sum.item.var)*num.item/((num.item-1)*var.scales)
     alpha.scale[is.nan(alpha.scale)] <- 1      #set to 1 for one item scales
    colnames(alpha.scale) <- slabels
   av.r <- alpha.scale/(num.item - alpha.scale*(num.item-1))  #alpha 1 = average r
   item.cor <- cor(items,scores,use="pairwise")
   if(is.null(ilabels)) {ilabels <-  paste("I",1:n.items,sep="")}
   rownames(item.cor) <- ilabels
   correction <- (colSums(abs(keys)-(keys))/2)*(max+min)
   scores <- scores  + matrix(rep(correction,n.subjects),byrow=TRUE,nrow=n.subjects)
   
   if (!totals) {scores <- scores \%*\% diag(1/num.item)   #find averages
                 colnames(scores) <- paste("A",1:n.keys,sep="") }
 scale.cor <- correct.cor(cor.scales,t(alpha.scale))
   if (missing){results <-list(scores=scores,missing = miss.rep,alpha=round(alpha.scale,digits), av.r=round(av.r,digits), n.items = num.item,  item.cor = round(item.cor,digits),cor = round(cor.scales,digits) ,corrected = round(scale.cor,digits))} else{  
    results <- list(scores=scores,alpha=round(alpha.scale,digits), av.r=round(av.r,digits), n.items = num.item,  item.cor = round(item.cor,digits), cor = round(cor.scales,digits),corrected = round(scale.cor,digits))} 
 }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ multivariate }
\keyword{models}

