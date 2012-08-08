\name{sim.multilevel}
\alias{sim.multilevel}
\title{Simulate multilevel data with specified within group and between group correlations}
\description{Multilevel data occur when observations are nested within groups. This can produce correlational structures that are sometimes difficult to understand. This simulation allows for demonstrations that correlations within groups do not imply, nor are implied by, correlations between group means. The correlations of aggregated data is sometimes called an 'ecological correlation'. That group level and individual level correlations are independent makes such inferences problematic.   
}
\usage{
sim.multilevel(nvar = 9, ngroups = 4, ncases = 16, rwg, rbg, eta)
}
\arguments{
  \item{nvar}{Number of variables to simulate}
  \item{ngroups}{The number of groups to simulate}
  \item{ncases}{The number of simulated cases}
  \item{rwg}{The within group correlational structure}
\item{rbg}{The between group correlational structure}
\item{eta}{The correlation of the data with the within data}
}
\details{The basic concepts of the independence of within group and between group correlations is discussed very clearly by Pedhazur (1997) as well as by Bliese (2009).  This function merely simulates pooled correlations (mixtures of between group and within group correlations) to allow for a better understanding of the problems inherent in multi-level modeling.  

Data (wg) are created with a particular within group structure (rwg).  Independent data (bg) are also created with a between group structure (rbg). Note that although there are ncases rows to this data matrix, there are only ngroups independent cases.  That is, every ngroups case is a repeat. The resulting data frame (xy) is a weighted sum of the wg and bg.  This is the inverse procedure for estimating estimating rwg and rbg from an observed rxy which is done by the \code{\link{statsBy}} function.  

 

}
\value{
\item{wg }{A matrix (ncases * nvar) of simulated within group scores}
\item{bg}{A matrix (ncases * nvar) of simulated between group scores}
\item{xy}{A matrix ncases * (nvar +1) of pooled data}
}


\references{
P. D. Bliese. Multilevel modeling in R (2.3) a brief introduction to R, the multilevel package and the nlme package, 2009.

Pedhazur, E.J. (1997) Multiple regression in behavioral research: explanation and prediction.  Harcourt Brace.

Revelle, W. An introduction to psychometric theory with applications in R (in prep) Springer. Draft chapters available at  \url{http://personality-project.org/r/book/}

}
\author{William Revelle}




\seealso{\code{\link{statsBy}}  for the decomposition of multi level data and \code{\link{withinBetween}} for an example data set.
}
\examples{
#get some parameters to simulate
data(withinBetween)
wb.stats <- statsBy(withinBetween,"Group")
rwg <- wb.stats$rwg
rbg <- wb.stats$rbg
eta <- rep(.5,9)

#simulate them.  Try this again to see how it changes
XY <- sim.multilevel(ncases=100,ngroups=10,rwg=rwg,rbg=rbg,eta=eta)
lowerCor(XY$wg)  #based upon 89 df
lowerCor(XY$bg)  #based upon 9 df   -- 
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{multivariate}
\keyword{models}
