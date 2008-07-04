\name{score.multiple.choice}
\alias{score.multiple.choice}

\title{Score multiple choice items and provide basic test statistics }
\description{Ability tests are typically multiple choice with one right answer.  score.multiple.choice takes a scoring key and a data matrix (or data.frame) and finds total or average number right for each participant.  Basic test statistics (alpha, average r, item means, item-whole correlations) are also reported. 
}
\usage{
score.multiple.choice(key, data, score = TRUE, totals = FALSE, ilabels = NULL, missing = TRUE, impute = "median", digits = 2,short=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{key}{ A vector of the correct item alternatives}
  \item{data}{a matrix or data frame of items to be scored.}
  \item{score}{score=FALSE, just convert to right (1) or wrong (0).\cr
               score=TRUE, find the totals or average scores and do item analysis}
  \item{totals}{total=FALSE: find the average number correct \cr
                total=TRUE: find the total number correct}
  \item{ilabels}{item labels }
  \item{missing}{missing=TRUE: missing values are replaced with means or medians \cr   mising=FALSE  missing values are not scored }
  \item{impute}{impute="median", replace missing items with the median score
               \cr impute="mean": replace missing values with the item mean}
  \item{digits}{ How many digits of output }
  \item{short}{short=TRUE, just report the item statistics, \cr
               short=FALSE, report item statistics and subject scores as well}
}
\details{Basically combines \code{\link{score.items}} with a conversion from multiple choice to right/wrong.

The item-whole correlation is inflated because of item overlap.}

\value{

  \item{scores }{Subject scores on one scale}
  \item{missing }{Number of missing items for each subject}
  \item{item.stats}{scoring key, response frequencies, item whole correlations, n subjects scored, mean, sd, skew, kurtosis and se for each item}
  \item{alpha}{Cronbach's coefficient alpha}
  \item{av.r}{Average interitem correlation}
}

\author{William Revelle}

\seealso{  \code{\link{score.items}}, \code{\link{omega}}}
\examples{
data(iqitems)
iq.keys <- c(4,4,3,1,4,3,2,3,1,4,1,3,4,3)
score.multiple.choice(iq.keys,iqitems)
#just convert the items to true or false 
iq.tf <- score.multiple.choice(iq.keys,iqitems,score=FALSE)
describe(iq.tf)  #compare to previous results

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{multivariate }
\keyword{models}
