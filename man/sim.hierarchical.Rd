\name{sim.hierarchical}
\alias{sim.hierarchical}
\alias{make.hierarchical}
\title{Create a population or sample correlation matrix, perhaps with hierarchical structure. }
\description{Create a population orthogonal or hierarchical correlation matrix from a set of factor loadings and factor intercorrelations. Samples of size n may be then be drawn from this population.  Return either the sample data, sample correlations, or population correlations.  This is used to create sample data sets for instruction and demonstration.
}
\usage{
sim.hierarchical(gload=NULL, fload=NULL, n = 0, raw = FALSE,mu = NULL)
make.hierarchical(gload=NULL, fload=NULL, n = 0, raw = FALSE)  #deprecated
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gload}{ Loadings of group factors on a general factor }
  \item{fload}{ Loadings of items on the group factors }
  \item{n}{ Number of subjects to generate: N=0 => population values }
  \item{raw}{ raw=TRUE, report the raw data, raw=FALSE, report the sample  correlation matrix. }
  \item{mu}{means for the individual variables}
}
\details{Many personality and cognitive tests have a hierarchical factor structure.  For demonstration purposes, it is useful to be able to create such matrices, either with population values, or sample values. 

Given a matrix of item factor loadings (fload) and of loadings of these factors on a general factor (gload), we create a population correlation matrix by using the general factor law (R = F' theta F where theta = g'g).  

To create sample values, we use  code adapted from the \code{\link{mvrnorm}} function in MASS.

The default is to return population correlation matrices. Sample correlation matrices are generated if n >0.  Raw data are returned if raw = TRUE.

The default values for gload and fload create a data matrix discussed by Jensen and Weng, 1994.

Although written to create hierarchical structures, if the gload matrix is all 0, then a non-hierarchical structure will be generated.
}
\value{ a matrix of correlations or a data matrix
}
\references{ \url{http://personality-project.org/r/r.omega.html }
\cr Jensen, A.R., Weng, L.J. (1994) What is a Good g? Intelligence, 18, 231-258.}
\author{ William Revelle }

\seealso{ \code{\link{omega}},   \code{\link{schmid}},  \code{\link{ICLUST}},  \code{\link{VSS}} for ways of analyzing these data.  Also see \code{\link{sim.structure}} to simulate a variety of structural models (e.g., multiple correlated factor models).  The simulation uses the  \code{\link{mvrnorm}} function from the MASS package.}
\examples{

gload <-  gload<-matrix(c(.9,.8,.7),nrow=3)    # a higher order factor matrix
fload <-matrix(c(                    #a lower order (oblique) factor matrix
           .8,0,0,
           .7,0,.0,
           .6,0,.0,
            0,.7,.0,
            0,.6,.0,
            0,.5,0,
            0,0,.6,
            0,0,.5,
            0,0,.4),   ncol=3,byrow=TRUE)
            
jensen <- sim.hierarchical(gload,fload)    #the test set used by omega
round(jensen,2)     

#simulate a non-hierarchical structure
fload <- matrix(c(c(c(.9,.8,.7,.6),rep(0,20)),c(c(.9,.8,.7,.6),rep(0,20)),
    c(c(.9,.8,.7,.6),rep(0,20)),c(c(c(.9,.8,.7,.6),rep(0,20)),c(.9,.8,.7,.6))),ncol=5)
gload <- matrix(rep(0,5))
five.factor <- sim.hierarchical(gload,fload,500,TRUE) #create sample data set
#do it again with a hierachical structure
gload <- matrix(rep(.7,5)  )
five.factor.g <- sim.hierarchical(gload,fload,500,TRUE) #create sample data set
#compare these two with omega
#not run
#om.5 <- omega(five.factor$observed,5)
#om.5g <- omega(five.factor.g$observed,5)
}
\keyword{ multivariate}% at least one, from doc/KEYWORDS
\keyword{models }% __ONLY ONE__ keyword per line
\keyword{datagen}
