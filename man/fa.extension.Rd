\name{fa.extension}
\alias{fa.extension}
\title{Apply Dwyer's factor extension to find factor loadings for extended variables}

\description{Dwyer (1937) introduced a method for finding factor loadings for variables not included in the original analysis.  This is basically finding the unattenuated correlation of the extension variables with the factor scores.  An alternative, which does not correct for factor reliability was proposed by Gorsuch (1997). This is an application of exploratory factor analysis with extensions to new variables.  
}
\usage{
fa.extension(Ro,Roe,fo,correct=TRUE)
}

\arguments{
  \item{Ro}{The correlation matrix that was factored}
  \item{Roe}{The correlations of the original variables with the extended variables
}
  \item{fo}{The output from the \code{\link{fa}} function.}
  \item{correct}{correct=TRUE produces Dwyer's solution, correct=FALSE produces Gorsuch's solution}

}
\details{It is assumed that the original factor solution was found by the \code{\link{fa}} function.
}
\value{Factor Loadings of the exended variables on the original factors}

\references{Paul S. Dwyer (1937), The determination of the factor loadings of a given test from the known factor loadings of other tests. Psychometrika, 3, 173-178

 Gorsuch, Richard L. (1997) New Procedure for Extension Analysis in Exploratory Factor Analysis  Educational and Psychological Measurement, 57, 725-740
}
\author{William Revelle
}




\seealso{ See Also as \code{\link{fa}}, \code{\link{principal}}, \code{\link{Dwyer}}
}
\examples{ #The Dwyer Example
Ro <- Dwyer[1:7,1:7]
Roe <- Dwyer[1:7,8]
fo <- fa(Ro,2,rotate="none")
fa.extension(Ro,Roe,fo)

#an example from simulated data
 d <- sim.item(12)
 R <- cor(d)
 Ro <- R[c(1,2,4,5,7,8,10,11),c(1,2,4,5,7,8,10,11)]
 Roe <- R[c(1,2,4,5,7,8,10,11),c(3,6,9,12)]
 fo <- fa(Ro,2)
 fa.extension(Ro,Roe,fo)
}
\keyword{ multivariate }% at least one, from doc/KEYWORDS