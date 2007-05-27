\name{fa.parallel}
\alias{fa.parallel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Scree plots of data or correlation matrix compared to random ``parallel" matrix }
\description{One way to determine the number of factors or components in a data matrix or a correlation matrix is to examine the ``scree" plot of the successive eigenvalues.  Sharp breaks in the plot suggest the appropriate number of components or factors to extract.  ``Parallel" analyis is an alternative technique that compares the scree of the observed data with that of a random data matrix of the same size as the original.  
}
\usage{
fa.parallel(x, ncases = 0, main = "Parallel Analysis Scree Plots")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ A data.frame or data matrix of scores.  If the matrix is square, it is assumed to be a correlation matrix.  Otherwise, correlations (with pairwise deletion) will be found }
  \item{ncases}{ncases=0 implies a data matrix/data.frame.  Otherwise, how many cases were used to find the correlations. }
  \item{main}{ a title for the analysis }
}
\details{Cattell's ``scree" test is one of most simple tests for the number of factors problem.  Humphreys and Montanelli's ``parallel" analysis is an equally compelling procedure.  Other procedures for determining the most optimal number of factors include finding the Very Simple Structure (VSS) criterion (\code{\link{VSS}}).
}
\value{
  A plot of the eigen values for the original data, a resampling of the  original data, and of a equivalent size matrix of random normal deviates.
}
\references{ }
\author{ William Revelle }


\seealso{ \code{\link{VSS}},\code{\link{VSS.plot}}, \code{\link{VSS.parallel}} }
\examples{

#not run
#test.data <- Harman74.cor$cov 
#fa.parallel(test.data,ncases=200)
#
#fa.parallel(attitude) 
#


## The function is currently defined as
function(x,ncases=0,main="Parallel Analysis Scree Plots")  {     
 nsub <- dim(x)[1]
 nvariables <- dim(x)[2]
 if (ncases >0) { nsub <- ncases 
  rx <- x} else {
  rx <- cor(x,use="pairwise")
  }
 if((ncases ==0) && (nsub==nvariables)) {warning("It seems as if you are using a correlation matrix, but have not specified the number of cases")}
  if(ncases == 0) {sampledata <- matrix(sample(unlist(x),size=nsub*nvariables,replace=TRUE),nrow=nsub,ncol=nvariables)
  values.samp <- eigen(cor(sampledata,use="pairwise"))$values}
  simdata=matrix(rnorm(nsub*nvariables),nrow=nsub,ncol=nvariables)    #make up simulated data
  values.sim <- eigen(cor(simdata))$values

 valuesx  <- eigen(rx)$values 

plot(valuesx,type="b", main = main ) 
points(values.sim,type ="b",lty="dotted")
if(ncases ==0) {points(values.samp,type ="b",lty="dashed")}
legend("topright", c("   Actual Data", "   Simulated Data", "  Resampled Data"), col = c(3,4,5),
       text.col = "green4", lty = c("solid","dotted", "dashed"),
       merge = TRUE, bg = 'gray90')
abline(h=1) }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{multivariate}

