\name{sim.item}
\alias{item.sim}
\alias{sim.item}
\alias{sim.dichot}
\alias{item.dichot}
\alias{sim.circ}
\alias{circ.sim}
\title{Generate simulated data structures for circumplex or simple structure }
\description{Rotations of factor analysis and principal components analysis solutions typically try to represent correlation matrices as simple structured.  An alternative structure, appealing to some, is a circumplex structure where the variables are uniformly spaced on the perimeter of a circle in a two dimensional space.  Generating simple structure and circumplex data is straightforward, and is useful for exploring alternative solutions to affect and personality structure.
}
\usage{
sim.item(nvar = 72, nsub = 500, circum = FALSE, xloading = 0.6, yloading = 0.6, gloading = 0, xbias = 0, ybias = 0, categorical = FALSE, low = -3, high = 3, truncate = FALSE, cutpoint = 0)
sim.circ(nvar = 72, nsub = 500, circum = TRUE, xloading = 0.6, yloading = 0.6, gloading = 0, xbias = 0, ybias = 0, categorical = FALSE, low = -3, high = 3, truncate = FALSE, cutpoint = 0)
sim.dichot(nvar = 72, nsub = 500, circum = FALSE, xloading = 0.6, yloading = 0.6, gloading = 0, xbias = 0, ybias = 0, low = 0, high = 0) 
item.dichot(nvar = 72, nsub = 500, circum = FALSE, xloading = 0.6, yloading = 0.6, gloading = 0, xbias = 0, ybias = 0, low = 0, high = 0) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nvar}{ Number of variables to simulate }
  \item{nsub}{Number of subjects to simulate }
  \item{circum}{ circum=TRUE is circumplex structure, FALSE is simple structure}
  \item{xloading}{the average loading on the first dimension }
  \item{yloading}{Average loading on the second dimension }
  \item{gloading}{Average loading on a general factor (default=0)}
  \item{xbias}{To introduce skew, how far off center is the first dimension }
  \item{ybias}{To introduce skew on the second dimension}
  \item{categorical}{ continuous or categorical variables.  }
  \item{low}{ values less than low are forced to low (or 0 in item.dichot)}
  \item{high}{ values greater than high are forced to high (or 1 in item.dichot) }
  \item{truncate}{Change all values less than cutpoint to cutpoint. }
  \item{cutpoint}{What is the cutpoint }
}
\details{This simulation was originally developed to compare the effect of skew on the measurement of affect (see Rafaeli and Revelle, 2005).  It has been extended to allow for a general simulation of affect or personality items with either a simple structure or a circumplex structure.  Items can be continuous normally distributed, or broken down into n categories (e.g, -2, -1, 0, 1, 2).  Items can be distorted by limiting them to these ranges, even though the items have a mean of (e.g., 1).  

The addition of item.dichot allows for testing structures with dichotomous items of different difficulty (endorsement) levels.  Two factor data with either simple structure or circumplex structure are generated for two sets of items, one giving a score of 1 for all items greater than the low (easy) value, one giving a 1 for all items greater than the high (hard) value. The default values for low and high are 0.  That is, all items are assumed to have a 50 percent endorsement rate.  To examine the effect of item difficulty, low could be   -1, high  1. This will lead to item endorsements of .84 for the easy and .16 for the hard.  Within each set of difficulties, the first 1/4 are assigned to the first factor factor, the second to the second factor, the third to the first factor (but with negative loadings) and the fourth to the second factor (but with negative loadings). 

It is useful to compare the results of sim.item with sim.hierarchical.  sim.item will produce a general factor that runs through all the items as well as two orthogonal factors.  This produces a data set that is hard to represent with standard rotation techniques.  Extracting 3 factors without rotation and then rotating the 2nd and 3rd factors reproduces the correct solution.  But simple oblique rotation of 3 factors, or an \code{\link{omega}} analysis do not capture the underlying structure.  See the last example.
}
\value{
  A data matrix of (nsub) subjects by (nvar) variables.  
}
\references{ Variations of a routine used in Rafaeli and Revelle, 2006;
Rafaeli, E. & Revelle, W. (2006). A premature consensus: Are happiness and sadness truly opposite affects? Motivation and Emotion. 

Acton, G. S. and Revelle, W. (2004) Evaluation of Ten Psychometric Criteria for Circumplex Structure. Methods of Psychological Research Online, Vol. 9, No. 1 \url{http://www.dgps.de/fachgruppen/methoden/mpr-online/issue22/mpr110_10.pdf} }
\author{ William Revelle }

\seealso{  See Also the implementation in this to generate numerous simulations.  \code{\link{simulation.circ}},  \code{\link{circ.tests}} as well as other simulations ( \code{\link{sim.structural}}  \code{\link{sim.hierarchical}})}
\examples{

round(cor(circ.sim(nvar=8,nsub=200)),2)
plot(fa(circ.sim(16,500),2)$loadings,main="Circumplex Structure") #circumplex structure
#
#
plot(fa(item.sim(16,500),2)$loadings,main="Simple Structure") #simple structure
#
cluster.plot(fa(item.dichot(16,low=0,high=1),2))

#not run
#x12 <- sim.item(12,gloading=.6)
#f3 <- fa(x12,3,rotate="none")
#f3  #observe the general factor
#oblimin(f3$loadings[,2:3])  #show the 2nd and 3 factors.
#f3 <- fa(x12,3)   #now do it with oblimin rotation
#f3  # not what one naively expect.

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ multivariate}
\keyword{datagen}

