\name{phi.demo}
\alias{phi.demo}
\title{ A simple demonstration of the Pearson, phi, and polychoric corelation} 
\description{A not very interesting demo of what happens if bivariate continuous data are dichotomized.  Bascially a demo of r, phi, and polychor.  }
\usage{
phi.demo(n=1000,r=.6, cuts=c(-2,-1,0,1,2))
}
\arguments{
  \item{n}{number of cases to simulate}
  \item{r}{ correlation between latent and observed }
  \item{cuts}{form dichotomized variables at the value of cuts}
  }
 
%- maybe also 'usage' for other objects documented here.
\details{A demonstration of the problem of different base rates on the phi correlation, and how these are partially solved by using the polychoric correlation. Not one of my more interesting demonstrations.  See 
\url{http://personality-project.org/r/simulating-personality.html} and 
\url{http://personality-project.org/r/r.datageneration.html} for better demonstrations of data generation.
}
\value{a matrix of correlations and a graphic plot.   The items above the diagonal are the tetrachoric correlations, below the diagonal are raw correlations.
  
}
\references{\url{http://personality-project.org/r/simulating-personality.html} and 
\url{http://personality-project.org/r/r.datageneration.html} for better demonstrations of data generation. }
\author{ William Revelle}

\seealso{\code{\link{VSS.simulate}},\code{\link{item.sim}}}
\examples{
if(require(polycor)) {demo <- phi.demo() #compare the phi (lower off diagonal and polychoric correlations (upper off diagonal)
#show the result from poly.mat
round(demo$tetrachoric,2)
#show the result from phi2poly
#tetrachorics above the diagonal, phi below the diagonal 
round(demo$phis,2) }
}
\keyword{ multivariate }% at least one, from doc/KEYWORDS
\keyword{ models }% __ONLY ONE__ keyword per line
