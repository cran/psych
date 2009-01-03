\name{sim.structural}
\alias{make.structural}
\alias{sim.structural}

\title{Create correlation matrices or data matrices with a particular measurement and structural model }
\description{Structural Equation Models decompose correlation or correlation matrices into a measurement (factor) model and a structural (regression) model.  sim.structural creates data sets with known measurement and structural properties. Population or sample correlation matrices with known properties are generated. Optionally raw data are produced. 

It is also possible to specify a measurement model for a set of x variables separately from a set of y variables.  They are then combined into one model.
}
\usage{
make.structural(s = NULL, f = NULL, n = 0, raw = FALSE)
sim.structural(s = NULL, f = NULL, n = 0,mx= NULL, my=NULL, raw = FALSE)
}
\arguments{
  \item{s}{ A correlation  (structure) matrix of latent variables. }
  \item{f}{ The measurement model}
  \item{mx}{The measurement model for the x variables}
  \item{my}{The measurement model for the y variables}
  \item{n}{ Number of cases to simulate.  If n=0, the population matrix is returned.}
  \item{raw}{if raw=TRUE, raw data are returned.}
 }
\details{Given the measurement model, f and the structure model f, the model is  f %*% s %*%  t(f).   Reliability is f %*% t(f).

Given the model, raw data are generated using the mvnorm function.

A special case of a structural model are one factor models such as parallel tests, tau equivalent tests, and congneneric tests.  These may be created by letting the structure matrix = 1 and then defining a vector of factor loadings. Alternatively, make.congeneric will do the same. 

If all parameters are NULL, a default example is produced.
}
\value{
  \item{model }{The correlation matrix} 
  \item{reliability }{The population reliability values}
  \item{data }{If raw=TRUE, a sample data matrix}
}
\references{Revelle, W. (in preparation) An Introduction to Psychometric Theory with applications in R. Springer. at \url{http://personality-project.org/r/book/}    }
 
\author{ William Revelle }


\seealso{ \code{\link{make.hierarchical}} for another structural model and \code{\link{make.congeneric}} for the one factor case.}
\examples{
gre.gpa <- sim.structural()
print(gre.gpa,2)  #the raw correlation matrix
round(correct.cor(gre.gpa$model,gre.gpa$reliability),2)  #correct for attenuation to see structure
congeneric <- sim.structural(1,c(.9,.8,.7,.6)) # a congeneric model 
congeneric 
}

\keyword{multivariate }
