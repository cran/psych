\name{sim.structure}
\alias{sim.structure}
\alias{sim.structural}
\alias{sim.correlation}
\alias{simCor}

\title{Create correlation matrices or data matrices with a particular measurement and structural model }
\description{Structural Equation Models decompose correlation or correlation matrices into a measurement (factor) model and a structural (regression) model.  sim.structural creates data sets with known measurement and structural properties. Population or sample correlation matrices with known properties are generated. Optionally raw data are produced. 

It is also possible to specify a measurement model for a set of x variables separately from a set of y variables.  They are then combined into one model with the correlation structure between the two sets.

Finally, the general case is given a population correlation matrix, generate data that will reproduce (with sampling variability) that correlation matrix.  \code{\link{sim.correlation}}.
}
\usage{
sim.structure(fx=NULL,Phi=NULL, fy=NULL, f=NULL, n=0, uniq=NULL, raw=TRUE, 
  items = FALSE, low=-2,high=2,d=NULL,cat=5, mu=0)
sim.structural(fx=NULL, Phi=NULL, fy=NULL, f=NULL, n=0, uniq=NULL, raw=TRUE,
      items = FALSE, low=-2,high=2,d=NULL,cat=5, mu=0)  #deprecated
simCor(R,n=1000,data=FALSE,scale=TRUE, skew=c("none","log","lognormal",
   "sqrt","abs"),vars=NULL,latent=FALSE,quant=NULL)
sim.correlation(R,n=1000,data=FALSE,scale=TRUE, skew=c("none","log","lognormal",
    "sqrt","abs"),vars=NULL,latent=FALSE,quant=NULL)

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
  \item{R}{The correlation matrix to reproduce}
  \item{data}{if TRUE, return the raw data, otherwise return the sample correlation matrix.}
  \item{scale}{standardize the simulated data?}
  \item{skew}{Defaults to none (the multivariate normal case.  Alternatives take the log, the squareroot, or the absolute value of latent or observed data )}
  \item{vars}{Apply the skewing or cuts to just these variables.  If NULL, to all the variables/}
  \item{latent}{Should the skewing transforms be applied to the latent variables, or the observed variables?}
  \item{quant}{Either a single number or a vector length nvar.  The data will be dichotomized at quant.}
 }
\details{Given the measurement model, fx and the structure model Phi, the model is  f \%*\% Phi \%*\%  t(f).   Reliability is f \%*\% t(f). \eqn{f \phi f'} and the reliability for each test is the items communality or just the diag of the model. 

If creating a correlation matrix, (uniq=NULL) then the diagonal is set to 1, otherwise the diagonal is diag(model) + uniq and the resulting structure is a covariance matrix.

%Given the model, raw data are generated using the mvnorm function.

A special case of a structural model are one factor models such as parallel tests, tau equivalent tests, and congeneric tests.  These may be created by letting the structure matrix = 1 and then defining a vector of factor loadings. Alternatively, \code{\link{sim.congeneric}} will do the same. 

The general case is to use \code{\link{simCor}} aka \code{\link{sim.correlation}} which will create data sampled from a specified correlation matrix for a particular sample size. If desired, it will just return the sample correlation matrix.  With data=TRUE, it will return the sample data as well.  It uses an eigen value decomposition of the original matrix times a matrix of random normal deviates (code adapted from the mvnorm function of Brian Ripley's MASS package).  These resulting scores may be transformed using a number of transforms (see the skew option) or made into dichotomous variables (see quant option) for all or a select set (vars option) of the variables.  

      
}
\value{
  \item{model }{The implied population correlation or covariance matrix} 
  \item{reliability }{The population reliability values}
  \item{r}{The sample correlation or covariance matrix} 
  \item{observed}{If raw=TRUE, a sample data matrix}
}
\references{Revelle, W. (in preparation) An Introduction to Psychometric Theory with applications in R. Springer. at \url{https://personality-project.org/r/book/}    }
 
\author{ William Revelle }


\seealso{ \code{\link{make.hierarchical}} for another structural model and \code{\link{make.congeneric}} for the one factor case. \code{\link{structure.list}} and \code{\link{structure.list}} for making symbolic structures.
} 
\examples{
#First, create a sem like model with a factor model of x and ys with correlation Phi
fx <-matrix(c( .9,.8,.6,rep(0,4),.6,.8,-.7),ncol=2)  
fy <- matrix(c(.6,.5,.4),ncol=1)
rownames(fx) <- c("V","Q","A","nach","Anx")
rownames(fy)<- c("gpa","Pre","MA")
Phi <-matrix( c(1,0,.7,.0,1,.7,.7,.7,1),ncol=3)
#now create this structure
gre.gpa <- sim.structural(fx,Phi,fy)
print(gre.gpa,2)  
#correct for attenuation to see structure
#the raw correlations are below the diagonal, the adjusted above
round(correct.cor(gre.gpa$model,gre.gpa$reliability),2) 

#These are the population values,
# we can also create a correlation matrix sampled from this population
GRE.GPA  <- sim.structural(fx,Phi,fy,n=250,raw=FALSE)
lowerMat(GRE.GPA$r)

#or we can show data sampled from such a population
GRE.GPA  <- sim.structural(fx,Phi,fy,n=250,raw=TRUE)
lowerCor(GRE.GPA$observed)


 
congeneric <- sim.structure(f=c(.9,.8,.7,.6)) # a congeneric model 
congeneric 

#now take this correlation matrix as a population value and create samples from it
example.congeneric <- sim.correlation(congeneric$model,n=200) #create a sample matrix
lowerMat(example.congeneric ) #show the correlation matrix
#or create another sample and show the data
example.congeneric.data <- simCor(congeneric$model,n=200,data=TRUE) 
describe(example.congeneric.data)
lowerCor(example.congeneric.data )
example.skewed <- simCor(congeneric$model,n=200,vars=c(1,2),data=TRUE,skew="log") 
describe(example.skewed)

}

\keyword{multivariate }
\keyword{datagen}
