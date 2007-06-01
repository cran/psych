\name{item.sim}
\alias{item.sim}
\alias{circ.sim}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Generate simulated data structures for circumplex or simple structure }
\description{Rotations of factor analysis and principal components analysis solutions typically try to represent correlation matrices as simple structured.  An alternative structure, appealing to some, is a circumplex structure where the variables are uniformly spaced on the perimeter of a circle in a two dimensional space.  Generating these data is straightforward, and is useful for exploring alternative solutions to affect and personality structure.
}
\usage{
item.sim(nvar = 72, nsub = 500, circum = FALSE, xloading = 0.6, yloading = 0.6, gloading = 0, xbias = 0, ybias = 0, categorical = FALSE, low = -3, high = 3, truncate = FALSE, cutpoint = 0)
circ.sim(nvar = 72, nsub = 500, circum = TRUE, xloading = 0.6, yloading = 0.6, gloading = 0, xbias = 0, ybias = 0, categorical = FALSE, low = -3, high = 3, truncate = FALSE, cutpoint = 0)
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
  \item{low}{ values less than low are forced to low }
  \item{high}{ values greater than high are forced to high }
  \item{truncate}{Change all values less than cutpoint to cutpoint. }
  \item{cutpoint}{What is the cutpoint }
}
\details{This simulation was originally developed to compare the effect of skew on the measurement of affect (see Rafaeli and Revelle, 2005).  It has been extended to allow for a general simulation of affect or personality items with either a simple structure or a circumplex structure.  Items can be continuous normally distributed, or broken down into n categories (e.g, -2, -1, 0, 1, 2).  Items can be distorted by limiting them to these ranges, even though the items have a mean of (e.g., 1).  
}
\value{
  A data matrix of (nsub) subjects by (nvar) variables.  
}
\references{ Variations of a routine used in Rafaeli and Revelle, 2006;
Rafaeli, E. & Revelle, W. (2006). A premature consensus: Are happiness and sadness truly opposite affects? Motivation and Emotion. 

Acton, G. S. and Revelle, W. (2004) Evaluation of Ten Psychometric Criteria for Circumplex Structure. Methods of Psychological Research Online, Vol. 9, No. 1 \url{http://www.dgps.de/fachgruppen/methoden/mpr-online/issue22/mpr110_10.pdf} }
\author{ William Revelle }

\seealso{  See Also the implementation in this to generate numerous simulations.  \code{\link{circ.simulation}},  \code{\link{circ.tests}} }
\examples{

round(cor(circ.sim(nvar=8,nsub=200)),2)
plot(factor.pa(circ.sim(16,500),2)$loadings) #circumplex structure
#
#
plot(factor.pa(item.sim(16,500),2)$loadings) #simple structure
#

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ multivariate}

