\name{irt.item.diff.rasch}
\alias{irt.item.diff.rasch}
\alias{irt.discrim}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Simple function to estimate item difficulties using IRT concepts}
\description{Steps toward a  very crude and preliminary IRT program. These two functions  estimate item difficulty and discrimination parameters.
}
\usage{
irt.item.diff.rasch(items)
irt.discrim(item.diff,theta,items)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{items}{ a matrix of items }
  \item{item.diff}{a vector of item difficulties (found by irt.item.diff)}
  \item{theta}{ability estimate from irt.person.theta}
}
\details{Item Response Theory (aka "The new psychometrics") models individual responses to items with a logistic function and an individual (theta) and item difficulty (diff) parameter.

irt.item.diff.rasch finds item difficulties with the assumption of theta=0 for all subjects and that all items are equally discriminating.

irt.discrim takes those  difficulties and theta estimates from \code{\link{irt.person.rasch}} to find item discrimination (beta) parameters.

A far better package with these features is the ltm package.  The IRT functions in the psych-package are for pedagogical rather than production purposes.  They are believed to be accurate, but are not guaranteed. They do seem to be slightly more robust to missing data structures associated with SAPA data sets than the ltm package. 
}
\value{a vector of item difficulties or item discriminations.

}
\references{ }
\author{William Revelle }
\note{ Under development. Not recommended for public consumption.}
\seealso{ \code{\link{irt.person.rasch}} }
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function(items) {
 ncases <- nrow(items)
 item.mean <- colMeans(items,na.rm=TRUE)
 item.mean[item.mean<(1/ncases)] <- 1/ncases
 irt.item.diff.rasch <- log((1/item.mean)- 1) }
 
 
 "irt.discrim" <- 
function(item.diff,theta,items) {
#find the item discrimination parameter  (beta)
#find the item discrimination parameter -- optimized by item.discrim
 irt.item.discrim <- function(x,diff,theta,scores) {
  fit <- -1*(log(scores/(1+exp(x*(diff-theta))) + (1-scores)/(1+exp(x*(theta-diff)))))
  mean(fit,na.rm=TRUE)
  }
 nitems <- length(item.diff)
 discrim <- matrix(NA,nitems,2)
 for (i in 1:nitems) {
    item.fit <- optimize(irt.item.discrim,c(-5,5),diff=item.diff[i],theta=theta,scores = items[,i])
   discrim[i,1] <- item.fit$minimum
   discrim[i,2] <- item.fit$objective}
  irt.discrim <- discrim
 }	
}
\keyword{ multivariate}% at least one, from doc/KEYWORDS
\keyword{ models}% __ONLY ONE__ keyword per line
