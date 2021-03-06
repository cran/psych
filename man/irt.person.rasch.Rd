\name{irt.1p}
\alias{irt.0p}
\alias{irt.1p}
\alias{irt.2p}
\alias{irt.person.rasch}
\title{Item Response Theory estimate of theta (ability) using a Rasch (like) model}
\description{Item Response Theory models individual responses to items by estimating individual ability (theta) and item difficulty (diff) parameters.  
This is an early and crude attempt to capture this modeling procedure. A better procedure is to use  \code{\link{irt.fa}}. 
}
\usage{
irt.person.rasch(diff, items)
irt.0p(items)
irt.1p(delta,items)
irt.2p(delta,beta,items) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{diff}{ A vector of item difficulties  --probably taken from irt.item.diff.rasch}
  \item{items}{A matrix of 0,1 items  nrows = number of subjects, ncols = number of items}
  \item{delta}{delta is the same as diff and is the item difficulty parameter}
  \item{beta}{beta is the item discrimination parameter  found in \code{\link{irt.discrim}} } 
}
\details{A very preliminary IRT estimation procedure.
Given scores xij for ith individual on jth item \cr
Classical Test Theory ignores item difficulty and defines ability as expected score : abilityi = theta(i) = x(i.)
A zero parameter model rescales these mean scores from 0 to 1 to a quasi logistic scale ranging from - 4 to 4
This is merely a non-linear transform of the raw data to reflect a logistic mapping.

Basic 1 parameter (Rasch) model considers item difficulties (delta j):  
	p(correct on item j for the ith subject |theta i, deltaj) = 1/(1+exp(deltaj - thetai))
	If we have estimates of item difficulty (delta), then we can find theta i by optimization 
	
Two parameter model adds item sensitivity (beta j):
	p(correct on item j for subject i |thetai, deltaj, betaj) = 1/(1+exp(betaj *(deltaj- theta i)))
Estimate delta, beta, and theta  to maximize fit of model to data.

The procedure used here is to first find the item difficulties assuming theta = 0
Then find theta given those deltas
Then find beta given delta and theta.

This is not an "official" way to do IRT, but is useful for basic item development. See \code{\link{irt.fa}} and \code{\link{score.irt}}  for far better options. 

}
\value{ a data.frame with estimated ability (theta) and quality of fit. (for irt.person.rasch)
   \cr
   a data.frame with the raw means, theta0, and the number of items completed}
  

\author{ William Revelle}
\note{ Not recommended for serious use.  This code is under development. Much better functions are in the ltm and eRm packages. Similar analyses can be done using \code{\link{irt.fa}} and \code{\link{score.irt}}. 
}
\seealso{ \code{\link{sim.irt}}, \code{\link{sim.rasch}},  \code{\link{logistic}}, \code{\link{irt.fa}}, \code{\link{tetrachoric}},  \code{\link{irt.item.diff.rasch}}   }

\keyword{multivariate}% at least one, from doc/KEYWORDS
\keyword{models}% __ONLY ONE__ keyword per line
