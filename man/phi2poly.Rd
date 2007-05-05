\name{phi2poly}
\alias{phi2poly}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Convert a phi coefficient to a polychoric correlation }
\description{Given a phi coefficient (a Pearson r calculated on two dichotomous variables), and the marginal frequencies, what is the corresponding estimate of the polychoric correlation?
}
\usage{
phi2poly(ph, cp, cc)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ph}{phi }
  \item{cp}{ probability of the predictor -- the so called selection ratio }
  \item{cc}{probability of the criterion -- the so called success rate. }
}
\details{Uses John Fox's polycor function.
}
\value{a polychoric correlation
}
\author{ William Revelle}
\seealso{  \code{\link{polychor.matrix}} }
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function(ph,cp,cc) {
	require(polycor)
     #ph is the phi coefficient
     #cp is the selection ratio of the predictor
     #cc is the success rate of the criterion
     r.marg<-rep(0,2)
     c.marg<- rep(0,2)
     p<-array(rep(0,4),dim=c(2,2))
     r.marg[1]<- cp
     r.marg[2]<- 1 -cp 
     c.marg[1]<- cc
     c.marg[2]<- 1-cc
     
     p[1,1]<- r.marg[1]*c.marg[1]+ ph*sqrt(prod(r.marg,c.marg))
     p[2,2]<- r.marg[2]*c.marg[2]+ ph*sqrt(prod(r.marg,c.marg))
     p[1,2]<- r.marg[1]*c.marg[2]- ph*sqrt(prod(r.marg,c.marg))
     p[2,1]<- r.marg[2]*c.marg[1]- ph*sqrt(prod(r.marg,c.marg))
     
     result<-polychor(p ) 
     return(result)}
}
\keyword{ models }% at least one, from doc/KEYWORDS
\keyword{ models }% __ONLY ONE__ keyword per line
