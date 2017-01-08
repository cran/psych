\name{test.irt}
\alias{test.irt}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{A simple demonstration (and test) of various IRT scoring algorthims.
}
\description{ Item Response Theory provides a number of alternative ways of estimating latent scores.  Here we compare 6 different ways to estimate the latent variable associated with a pattern of responses. Originally developed as a test for scoreIrt, but perhaps useful for demonstration purposes.  Items are simulated using \code{\link{sim.irt}} and then scored using factor scores from \code{\link{factor.scores}} using statistics found using \code{\link{irt.fa}}, simple weighted models for 1 and 2 PL and 2 PN. Results show almost perfect agreement with estimates from MIRT and ltm for the dichotomous case and with MIRT for the polytomous case.  (Results from ltm are unstable for the polytomous case, sometimes agreeing with \code{\link{scoreIrt}} and MIRT, sometimes being much worse.)  
}
\usage{
test.irt(nvar = 9, n.obs = 1000, mod = "logistic",type="tetra", low = -3, high = 3,
 seed = NULL)
}
\arguments{
  \item{nvar}{Number of variables to create (simulate) and score}
  \item{n.obs}{Number of simulated subjects}
  \item{mod}{"logistic" or "normal" theory data are generated}
  \item{type}{"tetra" for dichotomous, "poly" for polytomous}
  \item{low}{items range from low to high}
  \item{high}{items range from low to high}
  \item{seed}{Set the random number seed using some non-nul value.  Otherwise, use the existing sequence of random numbers}
}
\details{n.obs observations (0/1)  on nvar variables are simulated using either a logistic or normal theory model.  Then, a number of different scoring algorithms are applied and shown graphically.  Requires the ltm package to be installed to compare ltm scores.

}
\value{
A dataframe of scores as well as the generating theta true score. A graphic display of the correlations is also shown.}

\author{William Revelle}



\seealso{
 \code{\link{scoreIrt}},\code{\link{irt.fa}}
}
\examples{
#not run
#test.irt(9,1000)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{multivariate }% use one of  RShowDoc("KEYWORDS")
\keyword{models }% __ONLY ONE__ keyword per line
