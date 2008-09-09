\name{phi.demo}
\alias{phi.demo}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Create demo data for psychometrics }
\description{A not very interesting demo of what happens if bivariate continuous data are dichotomized.  Bascially a demo of r, phi, and polychor. }
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
\value{a matrix of correlations and a graphic plot)
  
}
\references{\url{http://personality-project.org/r/simulating-personality.html} and 
\url{http://personality-project.org/r/r.datageneration.html} for better demonstrations of data generation. }
\author{ William Revelle}

\seealso{\code{\link{VSS.simulate}},\code{\link{item.sim}}}
\examples{
round(phi.demo() ,2)  #compare the phi (lower off diagonal and polychoric correlations (upper off diagonal)
}
\keyword{ multivariate }% at least one, from doc/KEYWORDS
\keyword{ models }% __ONLY ONE__ keyword per line
