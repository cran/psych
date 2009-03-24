\name{score.items}
\alias{score.items}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Score item composite scales and find Cronbach's alpha, Guttman lambda 6 and item whole correlations }
\description{Given a matrix or data.frame of k keys for m items (-1, 0, 1), and a matrix or data.frame of items scores for m items and n people, find the sum scores or average scores for each person and each scale.  In addition, report Cronbach's alpha, Guttman's Lambda 6,  the average r, the scale intercorrelations, and the item by scale correlations (raw and corrected for item overlap).  Replace missing values with the item median or mean if desired. Will adjust scores for reverse scored items.  See \code{\link{make.keys}} for a convenient way to make the keys file.  If the input is a square matrix, then it is assumed that the input is a covariance or correlation matix and scores are not found, but the item statistics are reported.
}
\usage{
score.items(keys, items, totals = FALSE, ilabels = NULL, missing = TRUE,impute="median",  min = NULL, max = NULL, digits = 2,short=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{keys}{A matrix or dataframe of -1, 0, or 1 weights for each item on each scale. May be created by hand, or by using \code{\link{make.keys}} }
  \item{items}{ Matrix or dataframe of raw item scores }
  \item{totals}{ if TRUE  find total scores, if FALSE (default), find average scores }
  \item{ilabels}{ a vector of item labels. }
  \item{missing}{TRUE: Replace missing values with the corresponding item median or mean. FALSE: do not score that subject  }
  \item{impute}{impute="median" replaces missing values with the item median, impute = "mean" replaces values with the mean response.}
  \item{min}{May be specified as minimum item score allowed, else will be calculated from data}
  \item{max}{May be specified as maximum item score allowed, else will be calculated from data }
  \item{digits}{ Number of digits to report }
  \item{short}{if short is TRUE, then just give the item and scale statistics and do not report the scores}
}
\details{The process of finding sum or average scores for a set of scales given a larger set of items is a typical problem in psychometric research.  Although the structure of scales can be determined from the item intercorrelations, to find scale means, variances, and do further analyses, it is typical to find scores based upon the sum or the average item score.  For some strange reason, personality scale scores are typically given as totals, but attitude scores as averages.  The default for score.items is the average.  

Various estimates of scale reliability include ``Cronbach's alpha", Guttman's Lambda 6, and the average interitem correlation.  For k = number of items in a scale, and av.r = average correlation between items in the scale, alpha = k * av.r/(1+ (k-1)*av.r).  Thus, alpha is an increasing function of test length as well as the test homeogeneity.  

Alpha is a poor estimate of the general factor saturation of a test (see Zinbarg et al., 2005) for it can seriously overestimate the size of a general factor, and a better but not perfect estimate of total test reliability because it underestimates total reliability. None the less, it is a useful statistic to report.  To estimate the omega cofficient, use the \code{\link{omega}} function. 

Correlations between scales are attenuated by a lack of reliability.  Correcting correlations for reliability (by dividing by the square roots of the reliabilities of each scale) sometimes help show structure.  

By default, missing values are replaced with the corresponding median value for that item.  Means can be used instead (impute="mean"), or subjects with missing data can just be dropped (missing = FALSE).  
}
\value{ 
 \item{scores }{Sum or average scores for each subject on the k scales}
  \item{alpha }{Cronbach's coefficient alpha.  A simple (but non-optimal) measure of the internal consistency of a test. See also beta and omega. Set to 1 for scales of length 1. }
  \item{av.r}{The average correlation within a scale, also known as alpha 1 is a useful index of the internal consistency of a domain. Set to 1 for scales with 1 item.}
  \item{n.items}{Number of items on each scale}
   \item{item.cor}{The correlation of each item with each scale.  Because this is not corrected for item overlap, it will overestimate the amount that an item correlates with the other items in a scale.}
  \item{cor}{The intercorrelation of all the scales}
  \item{corrected}{The correlations of all scales (below the diagonal), alpha on the diagonal, and the unattenuated correlations (above the diagonal)}
  \item{item.corrected}{The item by scale correlations for each item, corrected for item overlap by replacing the item variance with the smc for that item}
 
}
\references{An introduction to psychometric theory with applications in R (in preparation).  \url{http://personality-project.org/r/book}}
\author{ William Revelle }

\seealso{\code{\link{make.keys}} for a convenient way to create the keys file,  \code{\link{score.multiple.choice}} for multiple choice items, \cr \code{\link{alpha.scale}}, \code{\link{correct.cor}}, \code{\link{cluster.cor}} , \code{\link{cluster.loadings}}, \code{\link{omega}}  for item/scale analysis}
\examples{

#see  the example including the bfi data set
data(bfi)
 keys.list <- list(agree=c(-1,2:5),conscientious=c(6:8,-9,-10),extraversion=c(-11,-12,13:15),neuroticism=c(16:20),openness = c(21,-22,23,24,-25))
 keys <- make.keys(25,keys.list,item.labels=colnames(bfi))
 scores <- score.items(keys,bfi)
 scores


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ multivariate }
\keyword{models}

