\name{bifactor}
\alias{bifactor}
\alias{Holzinger}
\alias{Reise}
\docType{data}
\title{Two data sets showing a bifactor solution.}
\description{Holzinger-Swineford (1937) introduced the bifactor model of a general factor and uncorrelated group factors. The Holzinger correlation matrix is a 14 * 14 matrix from their paper. The Reise data set is 16 * 16 correlation matrix of mental health items.  
}
\usage{data(bifactor)}

\details{Holzinger and Swineford introduced the bifactor model (one general factor and several group factors) for mental abilities.  This is a nice demonstration data set of a hierarchical factor structure that can be analyzed using the \code{\link{omega}} function or using sem. The bifactor model is typically used in measures of cognitive ability.

The 14 variables are ordered to reflect 3 spatial tests, 3 mental speed tests, 4 motor speed tests, and 4 verbal tests. 

More recent applications are to the measurement of psychological status. The Reise data set is a correlation matrix based upon >35,000 observations to the Consumer Assessment of Health Care Provideers and Systems survey instrument. Reise, Morizot, and Hays (2007) describe a bifactor solution based upon 1,000 cases.    

The five factors from Reise et al. reflect Getting care quickly (1-3), Doctor communicates well (4-7), Courteous and helpful staff (8,9), Getting needed care (10-13), and and Health plan customer service (14-16). 

}
\source{Holzinger:  Holzinger and Swineford (1937) \cr
Reise: Steve Reise (personal communication)
}
\references{ Holzinger, Karl and Swineford, Frances (1937) The Bi-factor method.  Psychometrika, 2, 41-54

Reise, Steven and Morizot, Julien and Hays, Ron (2007) The role of the bifactor model in resolving dimensionality issues in health outcomes measures. Quality of Life Research. 16, 19-31.

}
\examples{
data(bifactor)
holz <- omega(Holzinger,4, title = "14 ability tests from Holzinger-Swineford")
bf <- omega(Reise,5,title="16 health items from Reise") 
omega(Reise,5,labels=colnames(Reise),title="16 health items from Reise")

}
\keyword{datasets}
