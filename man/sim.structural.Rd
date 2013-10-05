\name{sim.structure}
\alias{sim.structure}
\alias{sim.structural}

\title{Create correlation matrices or data matrices with a particular measurement and structural model }
\description{Structural Equation Models decompose correlation or correlation matrices into a measurement (factor) model and a structural (regression) model.  sim.structural creates data sets with known measurement and structural properties. Population or sample correlation matrices with known properties are generated. Optionally raw data are produced. 

It is also possible to specify a measurement model for a set of x variables separately from a set of y variables.  They are then combined into one model with the correlation structure between the two sets.
}
\usage{
sim.structure(fx=NULL,Phi=NULL, fy=NULL, f=NULL, n=0, uniq=NULL, raw=TRUE, 
  items = FALSE, low=-2,high=2,d=NULL,cat=5, mu=0)
sim.structural(fx=NULL, Phi=NULL, fy=NULL, f=NULL, n=0, uniq=NULL, raw=TRUE,
      items = FALSE, low=-2,high=2,d=NULL,cat=5, mu=0)  #deprecated

}
\arguments{
  \item{fx}{The measurement model for x}
  \item{Phi}{The structure matrix of the latent variables}
  \item{fy}{The measurement model for y}
  \item{f}{ The measurement model}
  \item{n}{ Number of cases to simulate.  If n=0, the population matrix is returned.}
  \item{uniq}{The uniquenesses if creating a covariance matrix}
  \item{raw}{if raw=TRUE, raw data are returned as well for n > 0.}
  \item{items}{TRUE if simulating items, FALSE if simulating scales}
  \item{low}{Restrict the item difficulties to range from low to high}
  \item{high}{Restrict the item difficulties to range from low to high}
  \item{d}{A vector of item difficulties, if NULL will range uniformly from low to high}
  \item{cat}{Number of categories when creating binary (2) or polytomous items}
  \item{mu}{A vector of means, defaults to 0}
 }
\details{Given the measurement model, fx and the structure model Phi, the model is  f \%*\% Phi \%*\%  t(f).   Reliability is f \%*\% t(f). \eqn{f \phi f'} and the reliability for each test is the items communality or just the diag of the model. 

If creating a correlation matrix, (uniq=NULL) then the diagonal is set to 1, otherwise the diagonal is diag(model) + uniq and the resulting structure is a covariance matrix.

Given the model, raw data are generated using the mvnorm function.

A special case of a structural model are one factor models such as parallel tests, tau equivalent tests, and congeneric tests.  These may be created by letting the structure matrix = 1 and then defining a vector of factor loadings. Alternatively, make.congeneric will do the same. 

      
}
\value{
  \item{model }{The implied population correlation or covariance matrix} 
  \item{reliability }{The population reliability values}
  \item{r}{The sample correlation or covariance matrix} 
  \item{observed}{If raw=TRUE, a sample data matrix}
}
\references{Revelle, W. (in preparation) An Introduction to Psychometric Theory with applications in R. Springer. at \url{http://personality-project.org/r/book/}    }
 
\author{ William Revelle }


\seealso{ \code{\link{make.hierarchical}} for another structural model and \code{\link{make.congeneric}} for the one factor case. \code{\link{structure.list}} and \code{\link{structure.list}} for making symbolic structures.
} 
\examples{
fx <-matrix(c( .9,.8,.6,rep(0,4),.6,.8,-.7),ncol=2)  
fy <- matrix(c(.6,.5,.4),ncol=1)
rownames(fx) <- c("V","Q","A","nach","Anx")
rownames(fy)<- c("gpa","Pre","MA")
Phi <-matrix( c(1,0,.7,.0,1,.7,.7,.7,1),ncol=3)
gre.gpa <- sim.structural(fx,Phi,fy)
print(gre.gpa,2)  
#correct for attenuation to see structure

round(correct.cor(gre.gpa$model,gre.gpa$reliability),2)  
congeneric <- sim.structure(f=c(.9,.8,.7,.6)) # a congeneric model 
congeneric 

}

\keyword{multivariate }
\keyword{datagen}
