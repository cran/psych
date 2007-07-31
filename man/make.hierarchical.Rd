\name{make.hierarchical}
\alias{make.hierarchical}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Create a population or sample correlation matrix with hierachical structure. }
\description{Create a population hierarchical correlation matrix from a set of factor loadings and factor intercorrelations. Samples of size n may be then be drawn from this population.  Return either the sample data, sample correlations, or population correlations.  This is used to create sample data sets for instruction and demonstration.
}
\usage{
make.hierarchical(gload=NULL, fload=NULL, n = 0, raw = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gload}{ Loadings of group factors on a general factor }
  \item{fload}{ Loadings of items on the group factor }
  \item{n}{ Number of subjects to generate: N=0 => population values }
  \item{raw}{ raw=TRUE, report the raw data, raw=FALSE, report the sample  correlation matrix. }
}
\details{Many personality and cognitive tests have a hierarchical factor structure.  For demonstration purposes, it is useful to be able to create such matrices, either with population values, or sample values. 

Given a matrix of item factor loadings (fload) and of loadings of these factors on a general factor (gload), we create a population correlation matrix by using the general factor law (R = F' theta F where theta = g'g).  

To create sample values, we use the \code{\link{mvrnorm}} function from MASS.

The default is to return population correlation matrices. Sample correlation matrices are generated if n >0.  Raw data are returned if raw = TRUE.

The default values for gload and fload create a data matrix discussed by Jensen and Weng, 1994.
}
\value{ a matrix of correlations or a data matrix
}
\references{ \url{http://personality-project.org/r/r.omega.html }
\cr Jensen, A.R., Weng, L.J. (1994) What is a Good g? Intelligence, 18, 231-258.}
\author{ William Revelle }

\seealso{ \code{\link{omega}},   \code{\link{schmid}},  \code{\link{ICLUST}},  \code{\link{VSS}},  \code{\link{mvrnorm}} }
\examples{
\dontrun{
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
            
jensen <- make.hierarchical(gload,fload)    #the test set used by omega
round(jensen,2)     
}


}
\keyword{ multivariate}% at least one, from doc/KEYWORDS
\keyword{models }% __ONLY ONE__ keyword per line
