\name{fa.poly}
\alias{factor.pa} 
\alias{factor.minres} 
\alias{factor.wls}  
\alias{fa.poly}


\title{Deprecated Exploratory Factor analysis functions.  Please use fa}
\description{After 6  years, it is time to stop using these deprecated functions!  Please see \code{\link{fa}} which includes all of the functionality of these older functions.  
}

\usage{

fa.poly(x,nfactors=1,n.obs = NA, n.iter=1, rotate="oblimin", SMC=TRUE,  missing=FALSE,
 impute="median", min.err = .001, max.iter=50, symmetric=TRUE, warnings=TRUE,
 fm="minres",alpha=.1, p =.05,scores="regression", oblique.scores=TRUE,
        weight=NULL,global=TRUE,...)  #deprecated
        
factor.minres(r, nfactors=1, residuals = FALSE, rotate = "varimax",n.obs = NA,
scores = FALSE,SMC=TRUE, missing=FALSE,impute="median",min.err = 0.001, digits = 2,
 max.iter = 50,symmetric=TRUE,warnings=TRUE,fm="minres")    #deprecated

factor.wls(r,nfactors=1,residuals=FALSE,rotate="varimax",n.obs = NA,
scores=FALSE,SMC=TRUE,missing=FALSE,impute="median", min.err = .001,
 digits=2,max.iter=50,symmetric=TRUE,warnings=TRUE,fm="wls")  #deprecated
}

\arguments{
  \item{r}{deprecated.}
  \item{x}{deprecated}
  \item{nfactors}{ deprecated }
  \item{n.obs}{deprecated }
   \item{rotate}{deprecated}
 \item{n.iter}{deprecated}
  \item{residuals}{deprecated }
  \item{scores}{deprecated}
  \item{SMC}{deprecated }

  \item{missing}{deprecated}
  \item{impute}{deprecated}
  \item{max.iter}{deprecated}
  \item{symmetric}{deprecated}
  \item{warnings}{deprecated}
  \item{fm}{deprecated}
  \item{alpha}{deprecated}
  \item{p}{deprecated}
  \item{oblique.scores}{deprecated}
  \item{weight}{deprecated}
  
  \item{global}{deprecated}
  \item{digits}{deprecated}
  \item{min.err}{deprecated}
  \item{...}{deprecated}
  
 } 
 
\details{Please see the writeup for \code{\link{fa}} for all of the functionality in these older functions.
}
\value{
Pleases see the writeup for \code{\link{fa}}
  
 }



\references{
Revelle, William. (in prep) An introduction to psychometric theory with applications in R. Springer.  Working draft available at \url{https://personality-project.org/r/book/} 


}
\author{ William Revelle }

\note{These functions have been deprecated for 8 years.  Don't use them.




}
\examples{
#none, you should see fa
#using the Harman 24 mental tests, compare a principal factor with a principal components solution

}

\keyword{ multivariate }
\keyword{ models}