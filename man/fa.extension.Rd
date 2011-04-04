\name{fa.extension}
\alias{fa.extension}
\title{Apply Dwyer's factor extension to find factor loadings for extended variables}

\description{Dwyer (1937) introduced a method for finding factor loadings for variables not included in the original analysis.  This is basically finding the unattenuated correlation of the extension variables with the factor scores.  An alternative, which does not correct for factor reliability was proposed by Gorsuch (1997). Both options are an application of exploratory factor analysis with extensions to new variables.  
}
\usage{
fa.extension(Roe,fo,correct=TRUE)
}

\arguments{
  \item{Roe}{The correlations of the original variables with the extended variables}
  \item{fo}{The output from the \code{\link{fa}} or \code{\link{omega}} functions applied to the original variables.}
  \item{correct}{correct=TRUE produces Dwyer's solution, correct=FALSE produces Gorsuch's solution}

}
\details{It is sometimes the case that factors are derived from a set of variables (the Fo factor loadings) and we want to see what the loadings of an extended set of variables (Fe) would be. Given the original correlation matrix Ro and the correlation of these original variables with the extension variables of Roe, it is a straight forward calculation to find the loadings Fe of the extended variables on the original factors.  This technique was developed by Dwyer (1937) for the case of adding new variables to a factor analysis without doing all the work over again. But, as discussed by Horn (1973) factor extension is also appropriate when one does not want to include the extension variables in the original factor analysis, but does want to see what the loadings would be anyway.

This could be done by estimating the factor scores and then finding the covariances of the extension variables with the factor scores. But if the original data are not available, but just the covariance or correlation matrix is, then the use of \code{\link{fa.extension}} is most appropriate. 

The factor analysis results from either \code{\link{fa}} or \code{\link{omega}} functions applied to the original correlation matrix is extended to the extended variables given the correlations (Re) of the extended variables with the original variables.

\code{\link{fa.extension}} assumes that the original factor solution was found by the \code{\link{fa}} function.

For a very nice discussion of the relationship between factor scores, correlation matrices, and the factor loadings in a factor extension, see Horn (1973).
}
\value{Factor Loadings of the exended variables on the original factors}

\references{Paul S. Dwyer (1937), The determination of the factor loadings of a given test from the known factor loadings of other tests. Psychometrika, 3, 173-178

 Gorsuch, Richard L. (1997) New procedure for extension analysis in exploratory factor analysis,  Educational and Psychological Measurement, 57, 725-740
 
 Horn, John L. (1973) On extension analysis and ite relation to correlations between variables and factor scores.  Multivariate Behavioral Research, 8, (4), 477-489.
}
\author{William Revelle
}




\seealso{ See Also as \code{\link{fa}}, \code{\link{principal}}, \code{\link{Dwyer}}
}
\examples{ #The Dwyer Example
Ro <- Dwyer[1:7,1:7]
Roe <- Dwyer[1:7,8]
fo <- fa(Ro,2,rotate="none")
fa.extension(Roe,fo)

#an example from simulated data
set.seed(42) 
 d <- sim.item(12)    #two orthogonal factors 
 R <- cor(d)
 Ro <- R[c(1,2,4,5,7,8,10,11),c(1,2,4,5,7,8,10,11)]
 Roe <- R[c(1,2,4,5,7,8,10,11),c(3,6,9,12)]
 fo <- fa(Ro,2)
 fa.extension(Roe,fo)
 
 #create two correlated factors
 fx <- matrix(c(.9,.8,.7,.85,.75,.65,rep(0,12),.9,.8,.7,.85,.75,.65),ncol=2)
 Phi <- matrix(c(1,.6,.6,1),2)
 sim.data <- sim.structure(fx,Phi,n=1000,raw=TRUE)
 R <- cor(sim.data$observed)
 Ro <- R[c(1,2,4,5,7,8,10,11),c(1,2,4,5,7,8,10,11)]
 Roe <- R[c(1,2,4,5,7,8,10,11),c(3,6,9,12)]
 fo <- fa(Ro,2)
 fa.extension(Roe,fo)
 
 #an example of extending an omega analysis
 
 
fload <- matrix(c(c(c(.9,.8,.7,.6),rep(0,20)),c(c(.9,.8,.7,.6),rep(0,20)),c(c(.9,.8,.7,.6),rep(0,20)),c(c(c(.9,.8,.7,.6),rep(0,20)),c(.9,.8,.7,.6))),ncol=5)
 gload <- matrix(rep(.7,5))
 five.factor <- sim.hierarchical(gload,fload,500,TRUE) #create sample data set
 ss <- c(1,2,3,5,6,7,9,10,11,13,14,15,17,18,19)
 Ro <- cor(five.factor$observed[,ss])
 Re <- cor(five.factor$observed[,ss],five.factor$observed[,-ss])
 om5 <-omega(Ro,5)   #the omega analysis
 fa.extension(Re,om5) #the extension analysis
 
}
\keyword{ multivariate }% at least one, from doc/KEYWORDS